#pragma once

#include "particle.cuh"

/* Perform a physics-based update to the particle system over one timestep. */
void update(ParticleSystem *psystem, bool shake);