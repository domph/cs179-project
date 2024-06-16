#pragma once

#include "particle.h"

/* Perform a CUDA-accelerated, physics-based update to the particle system over one timestep. */
void cudaUpdate(ParticleSystem *psystem, ParticleSystem *gpu_psystem, bool shake);

/* Allocate memory to the CUDA psystem upon initial creation. */
void cudaMallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem);

/* Reallocate memory for the CUDA psystem upon parcel spawn. */
void cudaReallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem);