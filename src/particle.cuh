#pragma once

#include <cstring>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "constants.h"
#include <cuda_runtime.h>

/* Represents a box with a square base of length xybound, and height zbound */
struct Box {
    size_t xybound;
    size_t zbound;

    size_t x_partitions;
    size_t y_partitions;
    size_t z_partitions;
    size_t total_partitions;
    size_t num_particles;

    size_t *partitions;
    size_t *partition_sizes;

    Box(size_t xybound, size_t zbound, size_t num_particles);

    __host__ __device__ size_t part_idx(size_t x, size_t y, size_t z, size_t n);
    __host__ __device__ size_t get_id_at(size_t x, size_t y, size_t z, size_t n);

    __host__ __device__ size_t part_sz_idx(size_t x, size_t y, size_t z);
    __host__ __device__ size_t get_part_sz(size_t x, size_t y, size_t z);

    void add_particle(size_t id, glm::vec3 pos);
    void clear_partitions();
};

struct ParticleSystem {
    size_t num_particles;
    float t;
    float shake_t;

    glm::vec3 *pos;
    glm::vec3 *deltapos;
    glm::vec3 *prevpos;
    
    glm::vec3 *vel;
    glm::vec3 *nextvel;

    glm::vec3 *vorticity;
    float     *lambda;

    // Arranged as:
    // [p1n1, p2n1, ..., pNn1, p1n2, p2n2, ..., pNnM]
    size_t       *neighbors;
    
    size_t       *num_neighbors;

    Box *box;

    void init_particle(size_t id, glm::vec3 p, glm::vec3 v) {
        pos[id] = p;
        prevpos[id] = p;
        vel[id] = v;
    }

    void spawn_init() {
        size_t id = 0;
        PSYSTEM_INIT_SPAWN(init_particle(id++, glm::vec3(i, j, k*INIT_STEP), glm::vec3(0.0f)));
    }

    ParticleSystem() {
        t = 0;
        shake_t = 0;

        // Initial calculation of num_particles
        num_particles = 0;
        PSYSTEM_INIT_SPAWN(num_particles++);

        pos       = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        deltapos  = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        prevpos   = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        vel       = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        nextvel   = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        vorticity = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        lambda        = (float *) malloc(num_particles * sizeof(float));
        neighbors     = (size_t *)   malloc(num_particles * MAX_NEIGHBORS * sizeof(size_t));
        num_neighbors = (size_t *)   malloc(num_particles * sizeof(size_t));

        box = new Box((size_t)XYBOUND, (size_t)ZBOUND, num_particles);

        spawn_init();
    }

    void respawn() {
        // Recalculate num_particles
        num_particles = 0;
        PSYSTEM_INIT_SPAWN(num_particles++);

        free(box);
        box = new Box((size_t)XYBOUND, (size_t)ZBOUND, num_particles);
        
        t = 0;
        shake_t = 0;
        spawn_init();
    }

    void spawn_parcel(float x, float y, float z, float z_vel, float r) {
        size_t old_num_particles = num_particles;

        // Initial calculation of num_particles
        PSYSTEM_PARCEL_SPAWN(r, num_particles++);

        pos       = (glm::vec3 *) realloc(pos, num_particles * sizeof(glm::vec3));
        deltapos  = (glm::vec3 *) realloc(deltapos, num_particles * sizeof(glm::vec3));
        prevpos   = (glm::vec3 *) realloc(prevpos, num_particles * sizeof(glm::vec3));
        vel       = (glm::vec3 *) realloc(vel, num_particles * sizeof(glm::vec3));
        nextvel   = (glm::vec3 *) realloc(nextvel, num_particles * sizeof(glm::vec3));
        vorticity = (glm::vec3 *) realloc(vorticity, num_particles * sizeof(glm::vec3));
        lambda        = (float *) realloc(lambda, num_particles * sizeof(float));
        neighbors     = (size_t *)   realloc(neighbors, num_particles * MAX_NEIGHBORS * sizeof(size_t));
        num_neighbors = (size_t *)   realloc(num_neighbors, num_particles * sizeof(size_t));

        box->partitions = (size_t *) realloc(box->partitions, box->total_partitions * num_particles * sizeof(size_t));

        x += SHAKE(shake_t);
        glm::vec3 p = glm::vec3(x - r, y - r, z - r);
        glm::vec3 v = glm::vec3(0.0f, 0.0f, z_vel);
        PSYSTEM_PARCEL_SPAWN(r, init_particle(old_num_particles++, glm::vec3(i, j, k) + p, v));
    }
};