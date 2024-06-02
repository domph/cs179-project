#define GLM_ENABLE_EXPERIMENTAL

#include <cmath>
#include "physics.h"

// TODO: consider using "glm/gtx/fast_exponential.hpp"

float calcWpoly6(glm::vec3 i, glm::vec3 j) {    
    float r2 = glm::distance2(i, j);
    return (EPS <= r2) && (r2 <= P_H2) ? POLY6_COEFF * std::pow(P_H2 - r2, 3) : 0;
}

glm::vec3 calcWspiky(glm::vec3 i, glm::vec3 j) {
    float r = glm::l2Norm(i, j); // r = ||i - j||

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

    for (Particle *particle : psystem->particles) box->add_particle(particle);
}

void kNearestNeighbors(ParticleSystem *psystem) {
    Box *box = psystem->box;
    for (Particle *particle : psystem->particles) {
        particle->neighbors.clear();
        glm::vec3 pi = particle->pos;

        // Discretize particle positions into a 3D grid of size P_H
        int x = (pi.x + EPS) / P_H;
        int y = (pi.y + EPS) / P_H;
        int z = (pi.z + EPS) / P_H;

        float dist, max;
        int max_idx;
        glm::vec3 pj;

        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                for (int k = z - 1; k <= z + 1; k++) {
                    if (x >= box->x_partitions || y >= box->y_partitions ||
                        z >= box->z_partitions || x < 0 || y < 0 || z < 0) continue;
                    
                    for (Particle *neighbor: box->partitions[x][y][z]) {
                        if (neighbor == particle) continue;
                        pj = neighbor->pos;

                        dist = glm::distance(pi, pj);
                        if (particle->neighbors.size() < MAX_NEIGHBORS) {
                            if (dist < P_H) {
                                particle->neighbors.push_back(neighbor);
                            }
                        } else {
                            // max = 0.0f;
                            // for (int l = 0; l < particle->neighbors.size(); l++) {
                            //     float d = glm::distance(pi, particle->neighbors[l]->pos);
                            //     if (d > max) {
                            //         max = d;
                            //         max_idx = l;
                            //     }
                            // }
                            // if (dist < max && dist < P_H) {
                            //     particle->neighbors[max_idx] = neighbor;
                            // }
                        }
                    }
                }
            }
        }
    }
}

void applyBodyForces(ParticleSystem *psystem)
{
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {
            particle->vel.z -= G * DT;
            particle->pos = particle->prevpos + DT * particle->vel;
        }
    }
}

void calcLambda(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {
            glm::vec3 pi = particle->pos;

            glm::vec3 gradPjCi;  // Temporary store for calculated gradients

            // Accumulators for summing over neighbors j
            float rhoI = 0;
            float sumGradPkCi2 = 0;
            glm::vec3 sumGradPiCi = glm::vec3(0.0f);
            glm::vec3 pj;

            for (Particle *neighbor : particle->neighbors) {
                pj = neighbor->pos;

                rhoI += calcWpoly6(pi, pj);  // eq (1)
                gradPjCi = calcWspiky(pi, pj) * RHO_0_INV;  // eq (8)

                sumGradPiCi += gradPjCi;  // eq (9), denominator, k = j
                sumGradPkCi2 += glm::length2(gradPjCi);  // eq (9), k = i
            }
            sumGradPkCi2 += glm::length2(sumGradPiCi);   // eq (9), k = i

            float numerator = rhoI * RHO_0_INV - 1;  // eq (11)
            float denominator = sumGradPkCi2 + RELAXATION_EPS;  // eq (11)

            particle->lambda = -numerator / denominator;  // eq (11)
        }
    }
}

void calcDeltaPos(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {
            glm::vec3 pi = particle->pos;
            float lambda_i = particle->lambda;

            glm::vec3 dq = pi + DELTA_Q * glm::vec3(1.0f);
            float denom = 1.0f / calcWpoly6(pi, dq);  // eq (13)

            glm::vec3 dpi = glm::vec3(0.0f);
            glm::vec3 pj;
            float sCorrBase, sCorr, lambda_j;

            for (Particle *neighbor : particle->neighbors) {
                pj = neighbor->pos;
                lambda_j = neighbor->lambda;

                sCorrBase = calcWpoly6(pi, pj) * denom;
                sCorr = -SCORR_K * SCORR_N(sCorrBase);

                dpi += (lambda_i + lambda_j + sCorr) * calcWspiky(pi, pj);  // eq (14)
            }
            particle->deltapos = dpi * RHO_0_INV;  // eq (14)
        }
    }
}

void updatePos(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        particle->pos += particle->deltapos;
    }
}

void savePrevPos(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        particle->prevpos = particle->pos;
    }
}

void calcVel(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        particle->vel = (particle->pos - particle->prevpos) / DT;
        particle->nextvel = particle->vel;
    }
}

void updateVel(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        particle->vel = particle->nextvel;
    }
}

void applyCollisionResponse(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {            
            if (particle->pos.x < 0) {
                particle->pos.x *= -1;
                particle->vel.x *= -1;
            }
            if (particle->pos.x > psystem->box->w) {
                particle->pos.x = 2*psystem->box->w - particle->pos.x;
                particle->vel.x *= -1;
            }

            if (particle->pos.y < 0) {
                particle->pos.y *= -1;
                particle->vel.y *= -1;
            }
            if (particle->pos.y > psystem->box->w) {
                particle->pos.y = 2*psystem->box->w - particle->pos.y;
                particle->vel.x *= -1;
            }

            if (particle->pos.z < 0) {
                particle->pos.z *= -1;
                particle->vel.z *= -1;
            }
            if (particle->pos.z > psystem->box->h) {
                particle->pos.z = 2*psystem->box->h - particle->pos.z;
                particle->vel.z *= -1;
            }
        }
    }
}

void calcVorticityViscosity(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {
            glm::vec3 pi = particle->pos;
            glm::vec3 vi = particle->vel;

            glm::vec3 vij, pj;
            glm::vec3 wi    = glm::vec3(0.0f);
            glm::vec3 vXSPH = glm::vec3(0.0f);

            for (Particle *neighbor : particle->neighbors) {
                pj = neighbor->pos;
                vij = neighbor->vel - vi;

                wi += glm::cross(vij, calcWspiky(pi, pj));  // eq (15)
                vXSPH += vij * calcWpoly6(pi, pj);  // eq (17)
            }
            particle->nextvel += XSPH_C * vXSPH; // eq (17)
            particle->vorticity = wi;
        }
    }
}

void applyVorticityCorrection(ParticleSystem *psystem) {
    for (Particle *particle : psystem->particles) {
        if (!particle->fixed) {
            glm::vec3 pi = particle->pos;
            glm::vec3 gradwi = glm::vec3(0.0f);
            glm::vec3 pj, wj;

            for (Particle *neighbor : particle->neighbors) {
                pj = neighbor->pos;
                wj = neighbor->vorticity;

                gradwi += glm::length(wj) * calcWspiky(pi, pj);
            }

            /* Avoid normalizing zero values to 1 */
            if (glm::length(gradwi) > EPS) gradwi = glm::normalize(gradwi);

            particle->nextvel += DT * VORTICITY_EPS * glm::cross(gradwi, particle->vorticity);
        }
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