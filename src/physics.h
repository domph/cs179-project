#pragma once
#define GLM_ENABLE_EXPERIMENTAL

#include <glm/gtx/norm.hpp>
#include "particle.h"

// TODO: consider using "glm/gtx/fast_exponential.hpp"

/* Calculates the density smoothing kernel between positions i and j from
   Muller et al. (2003), Eq (20) */
float calcWpoly6(glm::vec3 i, glm::vec3 j);

/* Calculates the pressure smoothing kernel between positions i and j from
   Muller et al. (2003), Eq (21) */
glm::vec3 calcWspiky(glm::vec3 i, glm::vec3 j);

/* Clears partitions and partitions particles into cells. */
void calcPartition(ParticleSystem *particles);

/* k-nearest-neighbors algorithm for each particle in the particle system.
   For efficiency, only the particles in the 3x3x3 partition cells surrounding
   each particles are considered as candidate neighbors. A brute-force search
   is then run across these candidates to find the k nearest neighbors.*/
void kNearestNeighbors(ParticleSystem *particles);

/* Applies body forces to each particle - Algorithm 1 Line 2*/
void applyBodyForces(ParticleSystem *particles, float t);

/* Calculates the regularization parameter lambda for every particle from
   Macklin & Muller (2013), Eq (11) - Algorithm 1 Line 10 */
void calcLambda(ParticleSystem *particles);

/* Calculates the position update from Macklin & Muller (2013),
   Eq (11) - Algorithm 1 Line 13 */
void calcDeltaPos(ParticleSystem *particles);

/* Updates the position of all particles based on a previously calculated
   position update - Algorithm 1 Line 17 */
void updatePos(ParticleSystem *particles);

/* Saves the previous position of all particles - Algorithm 1 Line 23 */
void savePrevPos(ParticleSystem *particles);

/* Calculates the velocity of all particles based on the previous position,
   the current position and the simulation timestep - Algorithm 1 Line 21 */
void calcVel(ParticleSystem *particles);

/* Updates the velocity of all particles after accounting for vorticity
   and viscosity - Algorithm 1 Line 22*/
void updateVel(ParticleSystem *particles);

/* Apply a collision response between every particle in the particle
   system using an elastic model - Algorithm 1 Line 14 */
void applyCollisionResponse(ParticleSystem *particles);

/* Calculate the vorticity corrective force from Macklin & Muller (2013).
   Calculate XSPH viscosity and update the next velocity for each particle.
   Algorithm 1 Line 22 */
void calcVorticityViscosity(ParticleSystem *particles);

/* Apply the vorticity corrective force from Macklin & Muller (2013).
   Algorithm 1 Line 22 */
void applyVorticityCorrection(ParticleSystem *particles);

/* Perform a physics-based update to the particle system over one timestep. */
void update(ParticleSystem *particles, float t);