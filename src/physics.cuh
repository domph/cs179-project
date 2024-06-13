#pragma once

#include "particle.h"

/* Perform a physics-based update to the particle system over one timestep. */
void cudaUpdate(ParticleSystem *psystem, ParticleSystem *gpu_psystem, bool shake);