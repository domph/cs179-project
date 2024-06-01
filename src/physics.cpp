#define GLM_ENABLE_EXPERIMENTAL

#include <cmath>
#include <glm/gtx/norm.hpp>
#include "particle.h"

// TODO: consider using "glm/gtx/fast_exponential.hpp"

/* Calculates the density smoothing kernel between positions i and j from
   Muller et al. (2003), Eq (20) */
float calcWpoly6(glm::vec3 i, glm::vec3 j) {    
    float r2 = glm::distance2(i, j);
    return (EPS <= r2) && (r2 <= H2) ? POLY6_COEFF * std::pow(H2 - r2, 3) : 0;
}

/* Calculates the pressure smoothing kernel between positions i and j from
   Muller et al. (2003), Eq (21) */
glm::vec3 calcWspiky(glm::vec3 i, glm::vec3 j) {
    float r = glm::l2Norm(i, j); // r = ||i - j||

    if ((EPS <= r) && (r <= H)) {
        float coeff = -SPIKY_COEFF * std::pow(H - r, 2);
        return coeff * glm::normalize(i - j);
    } else {
        return glm::vec3(0.0f);
    }
}

/* Clears partitions and partitions particles into cells. */
void calcPartition(particle_system_t *particles) {
    box_t *box = particles->box;
    box->clear_partitions();

    for (int i = 0; i < particles->num_particles; i++) {
        box->add_particle(&particles->particles[i]);
    }
}

/* k-nearest-neighbors algorithm for each particle in the particle system.
   For efficiency, only the particles in the 3x3x3 partition cells surrounding
   each particles are considered as candidate neighbors. A brute-force search
   is then run across these candidates to find the k nearest neighbors.*/
void kNearestNeighbors(particle_system_t *particles) {
    for (int i = 0; i  < particles->num_particles; i++) {
        int neighbors_idx = i * MAX_NEIGHBORS;
        particle_t *particle = &particles->particles[i];
        glm::vec3 pi = particle->pos;
        box_t *box = particles->box;

        // Discretize particle positions into a 3D grid of size H
        int x = (pi.x + EPS) / H;
        int y = (pi.y + EPS) / H;
        int z = (pi.z + EPS) / H;

        float dist, max;
        int maxIdx;
        int num_neighbors = 0;

        for (int i = x - 1; i <= x + 1; i++) {
            for (int j = y - 1; j <= y + 1; j++) {
                for (int k = z - 1; k <= z + 1; k++) {
                    if (x >= box->x_partitions || y >= box->y_partitions ||
                        z >= box->z_partitions || x < 0 || y < 0 || z < 0) continue;
                    
                    for (particle_t *neighbor: box->partitions[x][y][z]) {
                        if (neighbor == particle) continue;
                        glm::vec3 pj = neighbor->pos;

                        dist = glm::distance(pi, pj);
                        if (num_neighbors < MAX_NEIGHBORS && dist < H) {
                            particles->neighbors[neighbors_idx + num_neighbors] = neighbor;
                            num_neighbors++;
                        } else {
                            max = 0.0f;
                            for (int l = 0; l < MAX_NEIGHBORS; l++) {
                                float d = glm::distance(pi, particles->neighbors[neighbors_idx + l].pos);
                                if (d > max) {
                                    max = d;
                                    maxIdx = l;
                                }
                            }
                            if (dist < max && dist < H) {
                                particles->neighbors[neighbors_idx + maxIdx] = neighbor;
                            }
                        }
                    }
                }
            }
        }
        particles->num_neighbors[i] = num_neighbors;
    }
}

/* Applies body forces to each particle - Algorithm 1 Line 2*/
void applyBodyForces(particle_system_t *particles)
{
    for (int i = 0; i < particles->num_particles; i++) {
        particles->particles[i].vel += G * DT;
        particles->particles[i].pos = particles->particles[i].prevpos +
                                      DT * particles->particles[i].vel;
    }
}

/* Calculates the regularization parameter lambda for every particle from
   Macklin & Muller (2013), Eq (11) - Algorithm 1 Line 10 */
void calcLambda(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        if (!particles->particles[i].fixed) {
            glm::vec3 pi = particles->particles[i].pos;

            glm::vec3 gradPjCi;  // Temporary store for calculated gradients

            // Accumulators for summing over neighbors j
            float rhoI = 0;
            float sumGradPkCi2 = 0;
            glm::vec3 sumGradPiCi = glm::vec3(0.0f);
            
            for (int j = 0; j < particles->num_neighbors[i]; j++) {
                glm::vec3 pj = particles->neighbors[i * MAX_NEIGHBORS + j].pos;

                rhoI += calcWpoly6(pi, pj);  // eq (1)
                gradPjCi = calcWspiky(pi, pj) * RHO_0_INV;  // eq (8)

                sumGradPiCi += gradPjCi;  // eq (9), denominator, k = j
                sumGradPkCi2 += glm::length2(gradPjCi);  // eq (9), k = i
            }
            sumGradPkCi2 += glm::length2(sumGradPiCi);   // eq (9), k = i

            float numerator = rhoI * RHO_0_INV - 1;  // eq (11)
            float denominator = sumGradPkCi2 + RELAXATION_EPS;  // eq (11)

            particles->particles[i].lambda = -numerator / denominator;  // eq (11)
        }
    }
}

/* Calculates the position update from Macklin & Muller (2013),
   Eq (11) - Algorithm 1 Line 13 */
void calcDeltaPos(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        if (!particles->particles[i].fixed) {
            glm::vec3 pi = particles->particles[i].pos;
            float lambdaI = particles->particles[i].lambda;

            glm::vec3 dq = pi + DELTA_Q * glm::vec3(1.0f);
            float denom = 1.0f / calcWpoly6(pi, dq);  // eq (13)

            glm::vec3 dpi = glm::vec3(0.0f);
            float sCorrBase, sCorr, lambdaJ;

            for (int j = 0; j < particles->num_neighbors[i]; j++) {
                glm::vec3 pj = particles->neighbors[i * MAX_NEIGHBORS + j].pos;
                lambdaJ = particles->particles[j].lambda;

                sCorrBase = calcWpoly6(pi, pj) * denom;
                sCorr = -SCORR_K * SCORR_N(sCorrBase);

                dpi += (lambdaI + lambdaJ + sCorr) * calcWspiky(pi, pj);  // eq (14)
            }
            particles->particles[i].deltapos = dpi * RHO_0_INV;  // eq (14)
        }
    }
}

/* Updates the position of all particles based on a previously calculated
   position update - Algorithm 1 Line 17 */
void updatePos(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        particles->particles[i].pos += particles->particles[i].deltapos;
    }
}

/* Saves the previous position of all particles - Algorithm 1 Line 23 */
void savePrevPos(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        particles->particles[i].prevpos = particles->particles[i].pos;

    }
}

/* Calculates the velocity of all particles based on the previous position,
   the current position and the simulation timestep - Algorithm 1 Line 21 */
void calcVel(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        particle_t *particle = &particles->particles[i];
        particle->vel = (particle->pos - particle->prevpos) / DT;
    }
}

/* Updates the velocity of all particles after accounting for vorticity
   and viscosity - Algorithm 1 Line 22*/
void updateVel(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        particles->particles[i].vel = particles->particles[i].nextvel;
    }
}

/* Apply a collision response between every particle in the particle
   system using an elastic model - Algorithm 1 Line 14 */
void applyCollisionResponse(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        if (!particles->particles[i].fixed) {
            particle_t *particle = &particles->particles[i];
            
            if (particle->pos.x < 0) {
                particle->pos.x *= -1;
                particle->vel.x *= -1;
            }
            if (particle->pos.x > particles->box->w) {
                particle->pos.x = 2*particles->box->w - particle->pos.x;
                particle->vel.x *= -1;
            }

            if (particle->pos.y < 0) {
                particle->pos.y *= -1;
                particle->vel.y *= -1;
            }
            if (particle->pos.y > particles->box->w) {
                particle->pos.y = 2*particles->box->w - particle->pos.y;
                particle->vel.x *= -1;
            }

            if (particle->pos.z < 0) {
                particle->pos.z *= -1;
                particle->vel.z *= -1;
            }
            if (particle->pos.z > particles->box->h) {
                particle->pos.z = 2*particles->box->h - particle->pos.z;
                particle->vel.z *= -1;
            }
        }
    }
}

/* Calculate the vorticity corrective force from Macklin & Muller (2013).
   Calculate XSPH viscosity and update the next velocity for each particle.
   Algorithm 1 Line 22 */
void calcVorticityViscosity(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        if (!particles->particles[i].fixed) {
            glm::vec3 pi = particles->particles[i].pos;
            glm::vec3 vi = particles->particles[i].vel;

            glm::vec3 vij;
            glm::vec3 wi    = glm::vec3(0.0f);
            glm::vec3 vXSPH = glm::vec3(0.0f);

            for (int j = 0; j < particles->num_neighbors[i]; j++) {
                glm::vec3 pj = particles->neighbors[i * MAX_NEIGHBORS + j].pos;
                glm::vec3 vj = particles->neighbors[i * MAX_NEIGHBORS + j].vel;
                vij = vj - vi;

                wi += glm::cross(vij, calcWspiky(pi, pj));  // eq (15)
                vXSPH += vij * calcWpoly6(pi, pj);  // eq (17)
            }
            particles->particles[i].nextvel = vi + XSPH_C * vXSPH; // eq (17)
            particles->particles[i].vorticity = wi;
        }
    }
}

/* Apply the vorticity corrective force from Macklin & Muller (2013).
   Algorithm 1 Line 22 */
void applyVorticityCorrection(particle_system_t *particles) {
    for (int i = 0; i < particles->num_particles; i++) {
        if (!particles->particles[i].fixed) {
            glm::vec3 pi = particles->particles[i].pos;
            glm::vec3 gradwi = glm::vec3(0.0f);

            for (int j = 0; j < particles->num_neighbors[i]; j++) {
                glm::vec3 pj = particles->neighbors[i * MAX_NEIGHBORS + j].pos;
                glm::vec3 wj = particles->neighbors[i * MAX_NEIGHBORS + j].vorticity;

                gradwi += glm::length(wj) * calcWspiky(pi, pj);
            }

            /* Avoid normalizing zero values to 1 */
            if (glm::length(gradwi) > EPS) gradwi = glm::normalize(gradwi);

            particles->particles[i].nextvel = DT * VORTICITY_EPS * glm::cross(gradwi,
                                              particles->particles[i].vorticity);
        }
    }
}

/* Perform a physics-based update to the particle system over one timestep. */
void update(particle_system_t *particles) {
    applyBodyForces(particles);
    calcPartition(particles);
    kNearestNeighbors(particles);

    for (int i = 0; i < SOLVER_ITERATIONS; i++) {
        calcLambda(particles);
        calcDeltaPos(particles);
        updatePos(particles);
        applyCollisionResponse(particles);
    }

    calcVel(particles);
    calcVorticityViscosity(particles);
    applyVorticityCorrection(particles);
    updateVel(particles);

    savePrevPos(particles);
}