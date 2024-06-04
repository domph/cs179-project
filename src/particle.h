#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "constants.h"


/* Represents a box with a xybound x xybound square base and height zbound */
struct Box {
    size_t xybound;
    size_t zbound;
    size_t *partitions;
    size_t *partition_sizes;

    size_t x_partitions;
    size_t y_partitions;
    size_t z_partitions;
    size_t total_partitions;
    size_t num_particles;

    Box(size_t xybound, size_t zbound, size_t num_particles) : xybound(xybound), zbound(zbound), num_particles(num_particles) {
        x_partitions = (float)xybound / P_H + 1;
        y_partitions = (float)xybound / P_H + 1;
        z_partitions = (float)zbound / P_H + 1;
        total_partitions = x_partitions * y_partitions * z_partitions;

        partitions = (size_t *) malloc(total_partitions * num_particles * sizeof(size_t));
        partition_sizes = (size_t *) calloc(total_partitions, sizeof(size_t));
    }

    size_t part_idx(size_t x, size_t y, size_t z, size_t n) {
        return x * y_partitions * z_partitions * num_particles +
                                    y * z_partitions * num_particles + 
                                    z * num_particles + n;
    }

    size_t get_id_at(size_t x, size_t y, size_t z, size_t n) {
        return partitions[part_idx(x, y, z, n)];
    }

    size_t part_sz_idx(size_t x, size_t y, size_t z) {
        return x * y_partitions * z_partitions + y * z_partitions + z;
    }

    size_t get_part_sz(size_t x, size_t y, size_t z) {
        return partition_sizes[part_sz_idx(x, y, z)];
    }

    void add_particle(size_t id, glm::vec3 pos) {
        int x = (pos.x + EPS) / P_H;
        int y = (pos.y + EPS) / P_H;
        int z = (pos.z + EPS) / P_H;

        if (x >= 0 && (size_t)x < x_partitions &&
            y >= 0 && (size_t)y < y_partitions &&
            z >= 0 && (size_t)z < z_partitions) {
            size_t n = partition_sizes[part_sz_idx(x, y, z)]++;
            partitions[part_idx(x, y, z, n)] = id;
        } else {
            printf("Box error: particle out of bounds!\n");
            printf("loc:   x: %d, y: %d, z: %d\n", x, y, z);
            printf("bound: x: %zu, y: %zu, z: %zu\n", x_partitions, y_partitions, z_partitions);
        }
    }

    void clear_partitions() {
        for (size_t i = 0; i < total_partitions; i++) partition_sizes[i] = 0;
    }
};

struct ParticleSystem {
    size_t num_particles;
    float t;

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

    void init_particle(size_t id, glm::vec3 p) {
        pos[id] = p;
        prevpos[id] = p;
        vel[id] = glm::vec3(0.0f);
    }

    void spawn_init() {
        size_t id = 0;
        PSYSTEM_INIT_SPAWN(init_particle(id++, glm::vec3(i, j, k*INIT_STEP)));
    }

    ParticleSystem() {
        t = 0;

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

        box = new Box(XYBOUND, ZBOUND, num_particles);

        spawn_init();
    }

    void respawn() {
        // Recalculate num_particles
        num_particles = 0;
        PSYSTEM_INIT_SPAWN(num_particles++);

        free(box);
        box = new Box(XYBOUND, ZBOUND, num_particles);
        
        spawn_init();
    }

    void spawn_parcel(float x, float y, float z) {
        size_t old_num_particles = num_particles;

        // Initial calculation of num_particles
        PSYSTEM_PARCEL_SPAWN(num_particles++);

        pos       = (glm::vec3 *) realloc(pos, num_particles * sizeof(glm::vec3));
        deltapos  = (glm::vec3 *) realloc(deltapos, num_particles * sizeof(glm::vec3));
        prevpos   = (glm::vec3 *) realloc(prevpos, num_particles * sizeof(glm::vec3));
        vel       = (glm::vec3 *) realloc(vel, num_particles * sizeof(glm::vec3));
        nextvel   = (glm::vec3 *) realloc(nextvel, num_particles * sizeof(glm::vec3));
        vorticity = (glm::vec3 *) realloc(vorticity, num_particles * sizeof(glm::vec3));
        lambda        = (float *) realloc(lambda, num_particles * sizeof(float));
        neighbors     = (size_t *)   realloc(neighbors, num_particles * MAX_NEIGHBORS * sizeof(size_t));
        num_neighbors = (size_t *)   realloc(num_neighbors, num_particles * sizeof(size_t));

        glm::vec3 p = glm::vec3(x - PARCEL_R, y - PARCEL_R, z - PARCEL_R);
        PSYSTEM_PARCEL_SPAWN(init_particle(old_num_particles++, glm::vec3(i, j, k) + p));
    }
};