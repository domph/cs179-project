#pragma once

#include "particle.h"

/* Perform a physics-based update to the particle system over one timestep. */
void update(ParticleSystem *psystem, bool shake);