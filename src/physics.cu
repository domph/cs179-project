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

__global__ void cudaKNearestNeighbors(ParticleSystem *psystem) {
    Box *box = psystem->box;

    size_t p = blockIdx.x * blockDim.x + threadIdx.x;
    while (p < psystem->num_particles) {
        glm::vec3 pi = psystem->pos[p];

        // Discretize particle positions into a 3D grid of size P_H
        int x = (int)((pi.x + EPS) / P_H);
        int y = (int)((pi.y + EPS) / P_H);
        int z = (int)((pi.z + EPS) / P_H);

        float dist, max;
        size_t max_idx = 0, num_neighbors = 0;
        glm::vec3 pj;

        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                for (int k = z - 1; k <= z + 1; k++) {
                    if (i < 0 || j < 0 || k < 0 || (size_t)i >= box->x_partitions ||
                        (size_t)j >= box->y_partitions || (size_t)k >= box->z_partitions) continue;
                    
                    for (size_t l = 0; l < box->get_part_sz(i, j, k); l++) {
                        size_t neighbor = box->get_id_at(i, j, k, l);
                        if (neighbor == p) continue;
                        pj = psystem->pos[neighbor];

                        dist = glm::distance(pi, pj);
                        if (num_neighbors < MAX_NEIGHBORS) {
                            if (dist < P_H) {
                                psystem->neighbors[p + num_neighbors * psystem->num_particles] = neighbor;
                                num_neighbors++;
                            }
                        } else {
                            max = 0.0f;
                            for (size_t idx = 0; idx < num_neighbors; idx++) {
                                size_t neighbor_idx = psystem->neighbors[p + idx * psystem->num_particles];
                                float d = glm::distance(pi, psystem->pos[neighbor_idx]);
                                if (d > max) {
                                    max = d;
                                    max_idx = idx;
                                }
                            }
                            if (dist < max && dist < P_H) {
                                psystem->neighbors[p + max_idx * psystem->num_particles] = neighbor;
                            }
                        }
                    }
                }
            }
        }
        psystem->num_neighbors[p] = num_neighbors;

        p += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcLambda(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        glm::vec3 pi = psystem->pos[i];

        glm::vec3 gradPjCi;  // Temporary store for calculated gradients

        // Accumulators for summing over neighbors j
        float rhoI = 0.0f;
        float sumGradPkCi2 = 0.0f;
        glm::vec3 sumGradPiCi = glm::vec3(0.0f);
        glm::vec3 pj;

        for (size_t j = 0; j < psystem->num_neighbors[i]; j++) {
            size_t neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];

            rhoI += cudaCalcWpoly6(pi, pj);  // eq (1)
            gradPjCi = cudaCalcWspiky(pi, pj) * RHO_0_INV;  // eq (8)

            sumGradPiCi += gradPjCi;  // eq (9), denominator, k = j
            sumGradPkCi2 += glm::length2(gradPjCi);  // eq (9), k = i
        }
        sumGradPkCi2 += glm::length2(sumGradPiCi);   // eq (9), k = i

        float numerator = rhoI * RHO_0_INV - 1.0f;  // eq (11)
        float denominator = sumGradPkCi2 + RELAXATION_EPS;  // eq (11)

        psystem->lambda[i] = -numerator / denominator;  // eq (11)

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcDeltaPos(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        glm::vec3 pi = psystem->pos[i];
        float lambda_i = psystem->lambda[i];

        glm::vec3 dq = pi + DELTA_Q * glm::vec3(1.0f);
        float denom = 1.0f / cudaCalcWpoly6(pi, dq);  // eq (13)

        glm::vec3 dpi = glm::vec3(0.0f);
        glm::vec3 pj;
        float sCorrBase, sCorr, lambda_j;

        for (size_t j = 0; j < psystem->num_neighbors[i]; j++) {
            size_t neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            lambda_j = psystem->lambda[neighbor_idx];

            sCorrBase = cudaCalcWpoly6(pi, pj) * denom;
            sCorr = -SCORR_K * SCORR_N(sCorrBase);

            dpi += (lambda_i + lambda_j + sCorr) * cudaCalcWspiky(pi, pj);  // eq (14)
        }
        psystem->deltapos[i] = dpi * RHO_0_INV;  // eq (14)

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaUpdatePos(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        psystem->pos[i] += psystem->deltapos[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaSavePrevPos(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        psystem->prevpos[i] = psystem->pos[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcVel(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        psystem->vel[i] = (psystem->pos[i] - psystem->prevpos[i]) / DT;
        psystem->nextvel[i] = psystem->vel[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaUpdateVel(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        psystem->vel[i] = psystem->nextvel[i];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaApplyCollisionResponse(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        if (psystem->pos[i].x < SHAKE(psystem->shake_t)) {
            psystem->pos[i].x = 2*SHAKE(psystem->shake_t) -psystem->pos[i].x;
            psystem->vel[i].x *= -1;
        }
        if (psystem->pos[i].x > psystem->box->xybound + SHAKE(psystem->shake_t)) {
            psystem->pos[i].x = 2*(psystem->box->xybound + SHAKE(psystem->shake_t)) - psystem->pos[i].x;
            psystem->vel[i].x *= -1;
        }

        if (psystem->pos[i].y < 0) {
            psystem->pos[i].y *= -1;
            psystem->vel[i].y *= -1;
        }
        if (psystem->pos[i].y > psystem->box->xybound) {
            psystem->pos[i].y = 2*psystem->box->xybound - psystem->pos[i].y;
            psystem->vel[i].x *= -1;
        }

        if (psystem->pos[i].z < 0) {
            psystem->pos[i].z *= -1;
            psystem->vel[i].z *= -1;
        }
        if (psystem->pos[i].z > psystem->box->zbound) {
            psystem->pos[i].z = 2*psystem->box->zbound - psystem->pos[i].z;
            psystem->vel[i].z *= -1;
        }

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaCalcVorticityViscosity(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        glm::vec3 pi = psystem->pos[i];
        glm::vec3 vi = psystem->vel[i];

        glm::vec3 vij, pj;
        glm::vec3 wi    = glm::vec3(0.0f);
        glm::vec3 vXSPH = glm::vec3(0.0f);

        for (size_t j = 0; j < psystem->num_neighbors[i]; j++) {
            size_t neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            vij = psystem->vel[neighbor_idx] - vi;

            wi += glm::cross(vij, cudaCalcWspiky(pi, pj));  // eq (15)
            vXSPH += vij * cudaCalcWpoly6(pi, pj);  // eq (17)
        }
        psystem->nextvel[i] += XSPH_C * vXSPH; // eq (17)
        psystem->vorticity[i] = wi;

        i += blockDim.x * gridDim.x;
    }
}

__global__ void cudaApplyVorticityCorrection(ParticleSystem *psystem) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    while (i < psystem->num_particles) {
        glm::vec3 pi = psystem->pos[i];
        glm::vec3 gradwi = glm::vec3(0.0f);
        glm::vec3 pj, wj;

        for (size_t j = 0; j < psystem->num_neighbors[i]; j++) {
            size_t neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            wj = psystem->vorticity[neighbor_idx];

            gradwi += glm::length(wj) * cudaCalcWspiky(pi, pj);
        }

        /* Avoid normalizing zero values to 1 */
        if (glm::length(gradwi) > EPS) gradwi = glm::normalize(gradwi);

        psystem->nextvel[i] += DT * VORTICITY_EPS * glm::cross(gradwi, psystem->vorticity[i]);

        i += blockDim.x * gridDim.x;
    }
}

void cudaUpdate(ParticleSystem *psystem, ParticleSystem *gpu_psystem, bool shake) {
    unsigned int blocks = 512;
    unsigned int threads_per_block = 512;

    // performed on CPU
    cudaApplyBodyForces(psystem);

    cudaMemcpy(gpu_psystem, psystem, offsetof(ParticleSystem, pos), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->pos,           psystem->pos,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->deltapos,      psystem->deltapos,      psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->prevpos,       psystem->prevpos,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->vel,           psystem->vel,           psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->nextvel,       psystem->nextvel,       psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->vorticity,     psystem->vorticity,     psystem->num_particles * sizeof(glm::vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->lambda,        psystem->lambda,        psystem->num_particles * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->neighbors,     psystem->neighbors,     psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->num_neighbors, psystem->num_neighbors, psystem->num_particles * sizeof(size_t), cudaMemcpyHostToDevice);
    
    // performed on CPU
    cudaCalcPartition(psystem);

    cudaMemcpy(gpu_psystem->box, psystem->box, offsetof(Box, partitions), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->box->partitions, psystem->box->partitions,
        psystem->box->total_partitions * psystem->num_particles * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_psystem->box->partition_sizes, psystem->box->partition_sizes,
        psystem->box->total_partitions * sizeof(size_t), cudaMemcpyHostToDevice);

    cudaKNearestNeighbors<<<blocks, threads_per_block>>>(gpu_psystem);

    for (size_t i = 0; i < SOLVER_ITERATIONS; i++) {
        cudaCalcLambda<<<blocks, threads_per_block>>>(gpu_psystem);
        cudaCalcDeltaPos<<<blocks, threads_per_block>>>(gpu_psystem);
        cudaUpdatePos<<<blocks, threads_per_block>>>(gpu_psystem);
        cudaApplyCollisionResponse<<<blocks, threads_per_block>>>(gpu_psystem);
    }

    cudaCalcVel<<<blocks, threads_per_block>>>(gpu_psystem);
    cudaCalcVorticityViscosity<<<blocks, threads_per_block>>>(gpu_psystem);
    cudaApplyVorticityCorrection<<<blocks, threads_per_block>>>(gpu_psystem);
    cudaUpdateVel<<<blocks, threads_per_block>>>(gpu_psystem);

    cudaSavePrevPos<<<blocks, threads_per_block>>>(gpu_psystem);


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



void cudaMallocPsystem(ParticleSystem *psystem, ParticleSystem **gpu_psystem) {
    CUDA_CALL(cudaMalloc(gpu_psystem, sizeof(ParticleSystem)));
    //CUDA_CALL(cudaMalloc(gpu_psystem, sizeof(ParticleSystem)));
    glm::vec3 *pos;
    CUDA_CALL(cudaMalloc(&pos,           psystem->num_particles * sizeof(glm::vec3)));
    CUDA_CALL(cudaMemcpy())
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->deltapos,      psystem->num_particles * sizeof(glm::vec3)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->prevpos,       psystem->num_particles * sizeof(glm::vec3)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->vel,           psystem->num_particles * sizeof(glm::vec3)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->nextvel,       psystem->num_particles * sizeof(glm::vec3)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->vorticity,     psystem->num_particles * sizeof(glm::vec3)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->lambda,        psystem->num_particles * sizeof(float)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->neighbors,     psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->num_neighbors, psystem->num_particles * sizeof(size_t)));

    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->box, sizeof(Box)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->box->partitions,
    //     psystem->box->total_partitions * psystem->num_particles * sizeof(size_t)));
    // CUDA_CALL(cudaMalloc(&(*gpu_psystem)->box->partition_sizes,
    //     psystem->box->total_partitions * sizeof(size_t)));
}

void cudaReallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem) {
    // cudaFree(gpu_psystem->pos);
    // cudaFree(gpu_psystem->deltapos);
    // cudaFree(gpu_psystem->prevpos);
    // cudaFree(gpu_psystem->vel);
    // cudaFree(gpu_psystem->nextvel);
    // cudaFree(gpu_psystem->vorticity);
    // cudaFree(gpu_psystem->lambda);
    // cudaFree(gpu_psystem->neighbors);
    // cudaFree(gpu_psystem->num_neighbors);

    // cudaMalloc(&gpu_psystem->pos,           psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->deltapos,      psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->prevpos,       psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->vel,           psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->nextvel,       psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->vorticity,     psystem->num_particles * sizeof(glm::vec3));
    // cudaMalloc(&gpu_psystem->lambda,        psystem->num_particles * sizeof(float));
    // cudaMalloc(&gpu_psystem->neighbors,     psystem->num_particles * MAX_NEIGHBORS * sizeof(size_t));
    // cudaMalloc(&gpu_psystem->num_neighbors, psystem->num_particles * sizeof(size_t));

    // cudaFree(gpu_psystem->box->partitions);
    
    // cudaMalloc(&gpu_psystem->box->partitions,
    //     psystem->box->total_partitions * psystem->num_particles * sizeof(size_t));
}