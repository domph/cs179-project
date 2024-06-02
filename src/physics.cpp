#define GLM_ENABLE_EXPERIMENTAL

#include <cmath>
#include "physics.h"

// TODO: consider using "glm/gtx/fast_exponential.hpp"

float calcWpoly6(glm::vec3 i, glm::vec3 j) {    
    float r2 = glm::distance2(i, j);
    return (EPS <= r2) && (r2 <= P_H2) ? POLY6_COEFF * std::pow(P_H2 - r2, 3) : 0;
}

glm::vec3 calcWspiky(glm::vec3 i, glm::vec3 j) {
    float r = glm::distance(i, j); // r = ||i - j||

    if ((EPS <= r) && (r <= P_H)) {
        float coeff = -SPIKY_COEFF * std::pow(P_H - r, 2);
        return coeff * glm::normalize(i - j);
    } else {
        return glm::vec3(0.0f);
    }
}

void calcPartition(ParticleSystem *psystem) {
    Box *box = psystem->box;
    box->clear_partitions();

    for (int i = 0; i < psystem->num_particles; i++) {
        box->add_particle(i, psystem->pos[i]);
    }
}

void kNearestNeighbors(ParticleSystem *psystem) {
    Box *box = psystem->box;
    for (int p = 0; p < psystem->num_particles; p++) {
        glm::vec3 pi = psystem->pos[p];

        // Discretize particle positions into a 3D grid of size P_H
        int x = (pi.x + EPS) / P_H;
        int y = (pi.y + EPS) / P_H;
        int z = (pi.z + EPS) / P_H;

        float dist, max;
        int max_idx, num_neighbors = 0;
        glm::vec3 pj;

        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                for (int k = z - 1; k <= z + 1; k++) {
                    if (x >= box->x_partitions || y >= box->y_partitions ||
                        z >= box->z_partitions || x < 0 || y < 0 || z < 0) continue;
                    
                    for (int neighbor: box->partitions[x][y][z]) {
                        if (neighbor == p) continue;  // should not happen
                        pj = psystem->pos[neighbor];

                        dist = glm::distance(pi, pj);
                        if (num_neighbors < MAX_NEIGHBORS) {
                            if (dist < P_H) {
                                psystem->neighbors[p + num_neighbors * psystem->num_particles] = neighbor;
                                num_neighbors++;
                            }
                        } else {
                            max = 0.0f;
                            for (int l = 0; l < MAX_NEIGHBORS; l++) {
                                int neighbor_idx = psystem->neighbors[p + l * psystem->num_particles];
                                float d = glm::distance(pi, psystem->pos[neighbor_idx]);
                                if (d > max) {
                                    max = d;
                                    max_idx = l;
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
    }
}

void applyBodyForces(ParticleSystem *psystem)
{
    for (int i = 0; i < psystem->num_particles; i++) {
        psystem->vel[i].z -= G * DT;
        psystem->pos[i] = psystem->prevpos[i] + DT * psystem->vel[i];
    }
}

void calcLambda(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        glm::vec3 pi = psystem->pos[i];

        glm::vec3 gradPjCi;  // Temporary store for calculated gradients

        // Accumulators for summing over neighbors j
        float rhoI = 0;
        float sumGradPkCi2 = 0;
        glm::vec3 sumGradPiCi = glm::vec3(0.0f);
        glm::vec3 pj;

        for (int j = 0; j < psystem->num_neighbors[i]; j++) {
            int neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];

            rhoI += calcWpoly6(pi, pj);  // eq (1)
            gradPjCi = calcWspiky(pi, pj) * RHO_0_INV;  // eq (8)

            sumGradPiCi += gradPjCi;  // eq (9), denominator, k = j
            sumGradPkCi2 += glm::length2(gradPjCi);  // eq (9), k = i
        }
        sumGradPkCi2 += glm::length2(sumGradPiCi);   // eq (9), k = i

        float numerator = rhoI * RHO_0_INV - 1;  // eq (11)
        float denominator = sumGradPkCi2 + RELAXATION_EPS;  // eq (11)

        psystem->lambda[i] = -numerator / denominator;  // eq (11)
    }
}

void calcDeltaPos(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        glm::vec3 pi = psystem->pos[i];
        float lambda_i = psystem->lambda[i];

        glm::vec3 dq = pi + DELTA_Q * glm::vec3(1.0f);
        float denom = 1.0f / calcWpoly6(pi, dq);  // eq (13)

        glm::vec3 dpi = glm::vec3(0.0f);
        glm::vec3 pj;
        float sCorrBase, sCorr, lambda_j;

        for (int j = 0; j < psystem->num_neighbors[i]; j++) {
            int neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            lambda_j = psystem->lambda[neighbor_idx];

            sCorrBase = calcWpoly6(pi, pj) * denom;
            sCorr = -SCORR_K * SCORR_N(sCorrBase);

            dpi += (lambda_i + lambda_j + sCorr) * calcWspiky(pi, pj);  // eq (14)
        }
        psystem->deltapos[i] = dpi * RHO_0_INV;  // eq (14)
    }
}

void updatePos(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        psystem->pos[i] += psystem->deltapos[i];
    }
}

void savePrevPos(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        psystem->prevpos[i] = psystem->pos[i];
    }
}

void calcVel(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        psystem->vel[i] = (psystem->pos[i] - psystem->prevpos[i]) / DT;
        psystem->nextvel[i] = psystem->vel[i];
    }
}

void updateVel(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        psystem->vel[i] = psystem->nextvel[i];
    }
}

void applyCollisionResponse(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {   
        if (psystem->pos[i].x < 0) {
            psystem->pos[i].x *= -1;
            psystem->vel[i].x *= -1;
        }
        if (psystem->pos[i].x > psystem->box->w) {
            psystem->pos[i].x = 2*psystem->box->w - psystem->pos[i].x;
            psystem->vel[i].x *= -1;
        }

        if (psystem->pos[i].y < 0) {
            psystem->pos[i].y *= -1;
            psystem->vel[i].y *= -1;
        }
        if (psystem->pos[i].y > psystem->box->w) {
            psystem->pos[i].y = 2*psystem->box->w - psystem->pos[i].y;
            psystem->vel[i].x *= -1;
        }

        if (psystem->pos[i].z < 0) {
            psystem->pos[i].z *= -1;
            psystem->vel[i].z *= -1;
        }
        if (psystem->pos[i].z > psystem->box->h) {
            psystem->pos[i].z = 2*psystem->box->h - psystem->pos[i].z;
            psystem->vel[i].z *= -1;
        }
    }
}

void calcVorticityViscosity(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        glm::vec3 pi = psystem->pos[i];
        glm::vec3 vi = psystem->vel[i];

        glm::vec3 vij, pj;
        glm::vec3 wi    = glm::vec3(0.0f);
        glm::vec3 vXSPH = glm::vec3(0.0f);

        for (int j = 0; j < psystem->num_neighbors[i]; j++) {
            int neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            vij = psystem->vel[neighbor_idx] - vi;

            wi += glm::cross(vij, calcWspiky(pi, pj));  // eq (15)
            vXSPH += vij * calcWpoly6(pi, pj);  // eq (17)
        }
        psystem->nextvel[i] += XSPH_C * vXSPH; // eq (17)
        psystem->vorticity[i] = wi;
    }
}

void applyVorticityCorrection(ParticleSystem *psystem) {
    for (int i = 0; i < psystem->num_particles; i++) {
        glm::vec3 pi = psystem->pos[i];
        glm::vec3 gradwi = glm::vec3(0.0f);
        glm::vec3 pj, wj;

        for (int j = 0; j < psystem->num_neighbors[i]; j++) {
            int neighbor_idx = psystem->neighbors[i + j * psystem->num_particles];
            pj = psystem->pos[neighbor_idx];
            wj = psystem->vorticity[neighbor_idx];

            gradwi += glm::length(wj) * calcWspiky(pi, pj);
        }

        /* Avoid normalizing zero values to 1 */
        if (glm::length(gradwi) > EPS) gradwi = glm::normalize(gradwi);

        psystem->nextvel[i] += DT * VORTICITY_EPS * glm::cross(gradwi, psystem->vorticity[i]);
    }
}

void update(ParticleSystem *psystem) {
    applyBodyForces(psystem);
    calcPartition(psystem);
    kNearestNeighbors(psystem);

    for (int i = 0; i < SOLVER_ITERATIONS; i++) {
        calcLambda(psystem);
        calcDeltaPos(psystem);
        updatePos(psystem);
        applyCollisionResponse(psystem);
    }

    calcVel(psystem);
    calcVorticityViscosity(psystem);
    applyVorticityCorrection(psystem);
    updateVel(psystem);

    savePrevPos(psystem);
}