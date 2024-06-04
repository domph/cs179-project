#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "constants.h"


/* Represents a box with a xybound x xybound square base and height zbound */
struct Box {
    int xybound;
    int zbound;
    std::vector<std::vector<std::vector<std::vector<int>>>> partitions;

    int x_partitions;
    int y_partitions;
    int z_partitions;
    int total_partitions;

    Box(int xybound, int zbound) : xybound(xybound), zbound(zbound) {
        x_partitions = (float)xybound / P_H + 1;
        y_partitions = (float)xybound / P_H + 1;
        z_partitions = (float)zbound / P_H + 1;
        total_partitions = x_partitions * y_partitions * z_partitions;

        partitions.resize(x_partitions);
        for (int x = 0; x < x_partitions; x++) {
            partitions[x].resize(y_partitions);
            for (int y = 0; y < y_partitions; y++) {
                partitions[x][y].resize(z_partitions);
            }
        }
    }

    void add_particle(int id, glm::vec3 pos) {
        int x = (pos.x + EPS) / P_H;
        int y = (pos.y + EPS) / P_H;
        int z = (pos.z + EPS) / P_H;

        if (x >= 0 && (size_t)x < partitions.size() &&
            y >= 0 && (size_t)y < partitions[x].size() &&
            z >= 0 && (size_t)z < partitions[x][y].size()) {
            partitions[x][y][z].push_back(id);
        } else {
            printf("Box error: particle out of bounds!\n");
            printf("loc:   x: %d, y: %d, z: %d\n", x, y, z);
            printf("bound: x: %zu, y: %zu, z: %zu\n", partitions.size(), partitions[x].size(), partitions[x][y].size());
        }
    }

    void clear_partitions() {
        for (auto &x_layer : partitions) {
            for (auto &y_layer : x_layer) {
                for (auto &z_layer : y_layer) {
                    z_layer.clear();
                }
            }
        }
    }

    // For debug only
    void print_members() {
        for (int x = 0; x < x_partitions; x++) {
            for (int y = 0; y < y_partitions; y++) {
                for (int z = 0; z < z_partitions; z++) {
                    printf("partition [%d][%d][%d] contains ", x, y, z);
                    for (int i : partitions[x][y][z]) {
                        printf("%d ", i);
                    }
                    printf("\n");
                }
            }
        }
    }
};

struct ParticleSystem {
    int num_particles;

    glm::vec3 *pos;
    glm::vec3 *deltapos;
    glm::vec3 *prevpos;
    
    glm::vec3 *vel;
    glm::vec3 *nextvel;

    glm::vec3 *vorticity;
    float     *lambda;

    // Arranged as:
    // [p1n1, p2n1, ..., pNn1, p1n2, p2n2, ..., pNnM]
    int       *neighbors;
    
    int       *num_neighbors;

    Box *box;

    void init_particle(int id, glm::vec3 p) {
        pos[id] = p;
        prevpos[id] = p;
        vel[id] = glm::vec3(0.0f);
    }

    void spawn_init() {
        int id = 0;
        PSYSTEM_INIT_SPAWN(init_particle(id++, glm::vec3(i, j, k*INIT_STEP)));
    }

    ParticleSystem(int xybound, int zbound) {
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
        neighbors     = (int *)   malloc(num_particles * MAX_NEIGHBORS * sizeof(int));
        num_neighbors = (int *)   malloc(num_particles * sizeof(int));

        box = new Box(xybound, zbound);

        spawn_init();
    }

    void respawn() {
        // Recalculate num_particles
        num_particles = 0;
        PSYSTEM_INIT_SPAWN(num_particles++);
        
        spawn_init();
    }

    void spawn_parcel(float x, float y, float z) {
        int old_num_particles = num_particles;

        // Initial calculation of num_particles
        PSYSTEM_PARCEL_SPAWN(num_particles++);

        pos       = (glm::vec3 *) realloc(pos, num_particles * sizeof(glm::vec3));
        deltapos  = (glm::vec3 *) realloc(deltapos, num_particles * sizeof(glm::vec3));
        prevpos   = (glm::vec3 *) realloc(prevpos, num_particles * sizeof(glm::vec3));
        vel       = (glm::vec3 *) realloc(vel, num_particles * sizeof(glm::vec3));
        nextvel   = (glm::vec3 *) realloc(nextvel, num_particles * sizeof(glm::vec3));
        vorticity = (glm::vec3 *) realloc(vorticity, num_particles * sizeof(glm::vec3));
        lambda        = (float *) realloc(lambda, num_particles * sizeof(float));
        neighbors     = (int *)   realloc(neighbors, num_particles * MAX_NEIGHBORS * sizeof(int));
        num_neighbors = (int *)   realloc(num_neighbors, num_particles * sizeof(int));

        glm::vec3 p = glm::vec3(x - PARCEL_R, y - PARCEL_R, z - PARCEL_R);
        PSYSTEM_PARCEL_SPAWN(init_particle(old_num_particles++, glm::vec3(i, j, k) + p));
    }
};