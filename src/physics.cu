#define GLM_ENABLE_EXPERIMENTAL

#include "physics.cuh"
#include <cmath>
#include <glm/gtx/norm.hpp>
#include <sstream>
#include "helper_cuda.h"

__device__ float cudaCalcWpoly6(glm::vec3 i, glm::vec3 j) {    
    float r2 = glm::distance2(i, j);
    return (float)(r2 <= P_H2 ? POLY6_COEFF * std::pow(P_H2 - r2, 3) : 0);
}

__device__ glm::vec3 cudaCalcWspiky(glm::vec3 i, glm::vec3 j) {
    float r = glm::distance(i, j);

    if ((EPS <= r) && (r <= P_H)) {
        float coeff = (float)(-SPIKY_COEFF * std::pow(P_H - r, 2));
        return coeff * glm::normalize(i - j);
    } else {
        return glm::vec3(0.0f);
    }
}

void cudaApplyBodyForces(ParticleSystem *psystem) {
    for (size_t i = 0; i < psystem->num_particles; i++) {
        psystem->vel[i].z -= G * DT;
        psystem->pos[i] = psystem->prevpos[i] + DT * psystem->vel[i];
    }
}

void cudaCalcPartition(ParticleSystem *psystem) {
    Box *box = psystem->box;
    box->clear_partitions();

    for (size_t i = 0; i < psystem->num_particles; i++) {
        box->add_particle(i, psystem->pos[i]);
    }
}

__global__ void cudaKNearestNeighbors(size_t num_particles, glm::vec3 *pos,
                                      size_t *neighbors, size_t *num_neighbors_arr,
                                      size_t x_partitions, size_t y_partitions,
                                      size_t z_partitions, size_t *partitions,
                                      size_t *partition_sizes
                                      ) {
    size_t p = blockIdx.x * blockDim.x + threadIdx.x;
    while (p < num_particles) {
        glm::vec3 pi = pos[p];

        // Discretize particle positions into a 3D grid of size P_H
        int x = (int)((pi.x + EPS) / P_H);
        int y = (int)((pi.y + EPS) / P_H);
        int z = (int)((pi.z + EPS) / P_H);

        float dist, max;
        size_t max_idx = 0, num_neighbors = 0, part_sz;
        glm::vec3 pj;

        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                for (int k = z - 1; k <= z + 1; k++) {
                    if (i < 0 || j < 0 || k < 0 || (size_t)i >= x_partitions ||
                        (size_t)j >= y_partitions || (size_t)k >= z_partitions) continue;

                    part_sz = partition_sizes[i * y_partitions * z_partitions + j * z_partitions + k];
                    
                    for (size_t l = 0; l < part_sz; l++) {
                        size_t neighbor = partitions[i * y_partitions *
                                    z_partitions * num_particles +
                                    j * z_partitions * num_particles + 
                                    k * num_particles + l];
                        
                        if (neighbor == p) continue;
                        pj = pos[neighbor];

                        dist = glm::distance(pi, pj);
                        if (num_neighbors < MAX_NEIGHBORS) {
                            if (dist < P_H) {
                                neighbors[p + num_neighbors * num_particles] = neighbor;
                                num_neighbors++;
                            }
                        } else {
                            max = 0.0f;
                            for (size_t idx = 0; idx < num_neighbors; idx++) {
                                size_t neighbor_idx = neighbors[p + idx * num_particles];
                                float d = glm::distance(pi, pos[neighbor_idx]);
                                if (d > max) {
                                    max = d;
                                    max_idx = idx;
                                }
                            }
                            if (dist < max && dist < P_H) {
                                neighbors[p + max_idx * num_particles] = neighbor;
                            }
                        }
                    }
                }
            }
        }
        num_neighbors_arr[p] = num_neighbors;

        p += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcLambda(size_t num_particles, glm::vec3 *pos, size_t *neighbors,
                               size_t *num_neighbors, float *lambda) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        glm::vec3 pi = pos[i];

        glm::vec3 gradPjCi;  // Temporary store for calculated gradients

        // Accumulators for summing over neighbors j
        float rhoI = 0.0f;
        float sumGradPkCi2 = 0.0f;
        glm::vec3 sumGradPiCi = glm::vec3(0.0f);
        glm::vec3 pj;

        for (size_t j = 0; j < num_neighbors[i]; j++) {
            size_t neighbor_idx = neighbors[i + j * num_particles];
            pj = pos[neighbor_idx];

            rhoI += cudaCalcWpoly6(pi, pj);  // eq (1)
            gradPjCi = cudaCalcWspiky(pi, pj) * RHO_0_INV;  // eq (8)

            sumGradPiCi += gradPjCi;  // eq (9), denominator, k = j
            sumGradPkCi2 += glm::length2(gradPjCi);  // eq (9), k = i
        }
        sumGradPkCi2 += glm::length2(sumGradPiCi);   // eq (9), k = i

        float numerator = rhoI * RHO_0_INV - 1.0f;  // eq (11)
        float denominator = sumGradPkCi2 + RELAXATION_EPS;  // eq (11)

        lambda[i] = -numerator / denominator;  // eq (11)

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcDeltaPos(size_t num_particles, glm::vec3 *pos, float *lambda,
                                 size_t *neighbors, size_t *num_neighbors, glm::vec3 *deltapos) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        glm::vec3 pi = pos[i];
        float lambda_i = lambda[i];

        glm::vec3 dq = pi + DELTA_Q * glm::vec3(1.0f);
        float denom = 1.0f / cudaCalcWpoly6(pi, dq);  // eq (13)

        glm::vec3 dpi = glm::vec3(0.0f);
        glm::vec3 pj;
        float sCorrBase, sCorr, lambda_j;

        for (size_t j = 0; j < num_neighbors[i]; j++) {
            size_t neighbor_idx = neighbors[i + j * num_particles];
            pj = pos[neighbor_idx];
            lambda_j = lambda[neighbor_idx];

            sCorrBase = cudaCalcWpoly6(pi, pj) * denom;
            sCorr = -SCORR_K * SCORR_N(sCorrBase);

            dpi += (lambda_i + lambda_j + sCorr) * cudaCalcWspiky(pi, pj);  // eq (14)
        }
        deltapos[i] = dpi * RHO_0_INV;  // eq (14)

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaUpdatePos(size_t num_particles, glm::vec3 *pos, glm::vec3 *deltapos) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        pos[i] += deltapos[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaSavePrevPos(size_t num_particles, glm::vec3 *prevpos, glm::vec3 *pos) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        prevpos[i] = pos[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcVel(size_t num_particles, glm::vec3 *vel, glm::vec3 *pos,
                            glm::vec3 *prevpos, glm::vec3 *nextvel) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        vel[i] = (pos[i] - prevpos[i]) / DT;
        nextvel[i] = vel[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaUpdateVel(size_t num_particles, glm::vec3 *vel, glm::vec3 *nextvel) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        vel[i] = nextvel[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaApplyCollisionResponse(size_t num_particles, glm::vec3 *pos,
                                           float shake_t, glm::vec3 *vel,
                                           size_t xybound, size_t zbound) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        if (pos[i].x < SHAKE(shake_t)) {
            pos[i].x = 2*SHAKE(shake_t) -pos[i].x;
            vel[i].x *= -1;
        }
        if (pos[i].x > xybound + SHAKE(shake_t)) {
            pos[i].x = 2*(xybound + SHAKE(shake_t)) - pos[i].x;
            vel[i].x *= -1;
        }

        if (pos[i].y < 0) {
            pos[i].y *= -1;
            vel[i].y *= -1;
        }
        if (pos[i].y > xybound) {
            pos[i].y = 2*xybound - pos[i].y;
            vel[i].x *= -1;
        }

        if (pos[i].z < 0) {
            pos[i].z *= -1;
            vel[i].z *= -1;
        }
        if (pos[i].z > zbound) {
            pos[i].z = 2*zbound - pos[i].z;
            vel[i].z *= -1;
        }

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcVorticityViscosity(size_t num_particles, glm::vec3 *pos,
                                           glm::vec3 *vel, size_t *neighbors,
                                           size_t *num_neighbors, glm::vec3 *nextvel,
                                           glm::vec3 *vorticity) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        glm::vec3 pi = pos[i];
        glm::vec3 vi = vel[i];

        glm::vec3 vij, pj;
        glm::vec3 wi    = glm::vec3(0.0f);
        glm::vec3 vXSPH = glm::vec3(0.0f);

        for (size_t j = 0; j < num_neighbors[i]; j++) {
            size_t neighbor_idx = neighbors[i + j * num_particles];
            pj = pos[neighbor_idx];
            vij = vel[neighbor_idx] - vi;

            wi += glm::cross(vij, cudaCalcWspiky(pi, pj));  // eq (15)
            vXSPH += vij * cudaCalcWpoly6(pi, pj);  // eq (17)
        }
        nextvel[i] += XSPH_C * vXSPH; // eq (17)
        vorticity[i] = wi;

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaApplyVorticityCorrection(size_t num_particles, glm::vec3 *pos,
                                             size_t *neighbors, size_t *num_neighbors,
                                             glm::vec3 *vorticity, glm::vec3 *nextvel) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < num_particles) {
        glm::vec3 pi = pos[i];
        glm::vec3 gradwi = glm::vec3(0.0f);
        glm::vec3 pj, wj;

        for (size_t j = 0; j < num_neighbors[i]; j++) {
            size_t neighbor_idx = neighbors[i + j * num_particles];
            pj = pos[neighbor_idx];
            wj = vorticity[neighbor_idx];

            gradwi += glm::length(wj) * cudaCalcWspiky(pi, pj);
        }

        /* Avoid normalizing zero values to 1 */
        if (glm::length(gradwi) > EPS) gradwi = glm::normalize(gradwi);

        nextvel[i] += DT * VORTICITY_EPS * glm::cross(gradwi, vorticity[i]);

        i += blockDim.x * gridDim.x;
    }
}

void cudaUpdate(ParticleSystem *psystem, ParticleSystem *gpu_psystem, bool shake) {
    // performed on CPU
    cudaApplyBodyForces(psystem);

    cudaMemcpy(gpu_psystem->pos,           psystem->pos,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->deltapos,      psystem->deltapos,      psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->prevpos,       psystem->prevpos,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->vel,           psystem->vel,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->nextvel,       psystem->nextvel,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->vorticity,     psystem->vorticity,     psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->lambda,        psystem->lambda,        psystem->num_particles * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->neighbors,     psystem->neighbors,     psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->num_neighbors, psystem->num_neighbors, psystem->num_particles * sizeof(size_t), cudaMemcpyHostToDevice);

    // performed on CPU to avoid race conditions
    cudaCalcPartition(psystem);

    cudaMemcpy(gpu_psystem->box->partitions, psystem->box->partitions,
        psystem->box->total_partitions * psystem->num_particles * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->box->partition_sizes, psystem->box->partition_sizes,
        psystem->box->total_partitions * sizeof(size_t), cudaMemcpyHostToDevice);

    cudaKNearestNeighbors<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                         gpu_psystem->pos,
                                                         gpu_psystem->neighbors,
                                                         gpu_psystem->num_neighbors,
                                                         psystem->box->x_partitions,
                                                         psystem->box->y_partitions,
                                                         psystem->box->z_partitions,
                                                         gpu_psystem->box->partitions,
                                                         gpu_psystem->box->partition_sizes
                                                         );

    for (size_t i = 0; i < SOLVER_ITERATIONS; i++) {
        cudaCalcLambda<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                      gpu_psystem->pos,
                                                      gpu_psystem->neighbors,
                                                      gpu_psystem->num_neighbors,
                                                      gpu_psystem->lambda);

        cudaCalcDeltaPos<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                        gpu_psystem->pos,
                                                        gpu_psystem->lambda,
                                                        gpu_psystem->neighbors,
                                                        gpu_psystem->num_neighbors,
                                                        gpu_psystem->deltapos);

        cudaUpdatePos<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                     gpu_psystem->pos,
                                                     gpu_psystem->deltapos);

        cudaApplyCollisionResponse<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                                  gpu_psystem->pos,
                                                                  psystem->shake_t,
                                                                  gpu_psystem->vel,
                                                                  psystem->box->xybound,
                                                                  psystem->box->zbound);
    }

    cudaCalcVel<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                               gpu_psystem->vel,
                                               gpu_psystem->pos,
                                               gpu_psystem->prevpos,
                                               gpu_psystem->nextvel);

    cudaCalcVorticityViscosity<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                              gpu_psystem->pos,
                                                              gpu_psystem->vel,
                                                              gpu_psystem->neighbors,
                                                              gpu_psystem->num_neighbors,
                                                              gpu_psystem->nextvel,
                                                              gpu_psystem->vorticity);
    
    cudaApplyVorticityCorrection<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                                gpu_psystem->pos,
                                                                gpu_psystem->neighbors,
                                                                gpu_psystem->num_neighbors,
                                                                gpu_psystem->vorticity,
                                                                gpu_psystem->nextvel);
    
    cudaUpdateVel<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                 gpu_psystem->vel,
                                                 gpu_psystem->nextvel);

    cudaSavePrevPos<<<BLOCKS, THREADS_PER_BLOCK>>>(psystem->num_particles,
                                                   gpu_psystem->prevpos,
                                                   gpu_psystem->pos);
    

    cudaMemcpy(psystem->pos,           gpu_psystem->pos,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->deltapos,      gpu_psystem->deltapos,      psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->prevpos,       gpu_psystem->prevpos,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->vel,           gpu_psystem->vel,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->nextvel,       gpu_psystem->nextvel,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->vorticity,     gpu_psystem->vorticity,     psystem->num_particles * sizeof(glm::vec3), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->lambda,        gpu_psystem->lambda,        psystem->num_particles * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->neighbors,     gpu_psystem->neighbors,     psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(psystem->num_neighbors, gpu_psystem->num_neighbors, psystem->num_particles * sizeof(size_t), cudaMemcpyDeviceToHost);


    psystem->t += DT;
    if (shake) psystem->shake_t += DT;
}


void cudaMallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem) {
    glm::vec3 *pos;
    CUDA_CALL(cudaMalloc(&pos, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->pos = pos;

    glm::vec3 *deltapos;
    CUDA_CALL(cudaMalloc(&deltapos, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->deltapos = deltapos;

    glm::vec3 *prevpos;
    CUDA_CALL(cudaMalloc(&prevpos, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->prevpos = prevpos;

    glm::vec3 *vel;
    CUDA_CALL(cudaMalloc(&vel, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->vel = vel;

    glm::vec3 *nextvel;
    CUDA_CALL(cudaMalloc(&nextvel, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->nextvel = nextvel;

    glm::vec3 *vorticity;
    CUDA_CALL(cudaMalloc(&vorticity, psystem->num_particles * sizeof(glm::vec3)));
    gpu_psystem->vorticity = vorticity;

    float *lambda;
    CUDA_CALL(cudaMalloc(&lambda, psystem->num_particles * sizeof(float)));
    gpu_psystem->lambda = lambda;
    
    size_t *neighbors;
    CUDA_CALL(cudaMalloc(&neighbors, psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t)));
    gpu_psystem->neighbors = neighbors;

    size_t *num_neighbors;
    CUDA_CALL(cudaMalloc(&num_neighbors, psystem->num_particles * sizeof(size_t)));
    gpu_psystem->num_neighbors = num_neighbors;

    size_t *partitions;
    CUDA_CALL(cudaMalloc(&partitions, psystem->box->total_partitions * psystem->num_particles * sizeof(size_t)));
    gpu_psystem->box->partitions = partitions;

    size_t *partition_sizes;
    CUDA_CALL(cudaMalloc(&partition_sizes, psystem->box->total_partitions * sizeof(size_t)));
    gpu_psystem->box->partition_sizes = partition_sizes;
}

void cudaReallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem) {
    cudaFree(gpu_psystem->pos);
    cudaFree(gpu_psystem->deltapos);
    cudaFree(gpu_psystem->prevpos);
    cudaFree(gpu_psystem->vel);
    cudaFree(gpu_psystem->nextvel);
    cudaFree(gpu_psystem->vorticity);
    cudaFree(gpu_psystem->lambda);
    cudaFree(gpu_psystem->neighbors);
    cudaFree(gpu_psystem->num_neighbors);
    cudaFree(gpu_psystem->box->partitions);
    cudaFree(gpu_psystem->box->partition_sizes);

    cudaMallocPsystem(psystem, gpu_psystem);
}