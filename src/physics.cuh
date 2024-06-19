#pragma once

#include "particle.h"

/* Perform a CUDA-accelerated, physics-based update to the particle system over one timestep. */
void cudaUpdate(ParticleSystem *psystem, ParticleSystem *gpu_psystem, bool shake, bool copy_to_device, bool copy_to_host);

/* Allocate memory to the CUDA psystem upon initial creation. */
void cudaMallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem);

/* Reallocate memory for the CUDA psystem upon parcel spawn. */
void cudaReallocPsystem(ParticleSystem *psystem, ParticleSystem *gpu_psystem);

/* Copy physics data from host to device */
void cudaCopyHostToDevice(ParticleSystem *psystem, ParticleSystem *gpu_psystem);

/* Copy physics data from device to host */
void cudaCopyDeviceToHost(ParticleSystem *psystem, ParticleSystem *gpu_psystem);