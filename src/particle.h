#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "constants.h"


/* Represents a box with a w x w square base and height h */
struct Box {
    int h;
    int w;
    std::vector<std::vector<std::vector<std::vector<int>>>> partitions;

    int x_partitions;
    int y_partitions;
    int z_partitions;
    int total_partitions;

    Box(int h, int w) : h(h), w(w) {
        x_partitions = (float)w / P_H + 1;
        y_partitions = (float)w / P_H + 1;
        z_partitions = (float)h / P_H + 1;
        total_partitions = x_partitions * y_partitions * z_partitions;

        printf("x partitions: %d\n", x_partitions);
        printf("y partitions: %d\n", y_partitions);
        printf("z partitions: %d\n", z_partitions);
        printf("total partitions: %d\n", total_partitions);

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

        if (x >= 0 && x < partitions.size() &&
            y >= 0 && y < partitions[x].size() &&
            z >= 0 && z < partitions[x][y].size()) {
            partitions[x][y][z].push_back(id);
        } else {
            printf("Box error: particle out of bounds!\n");
            printf("loc:   x: %d, y: %d, z: %d\n", x, y, z);
            printf("bound: x: %d, y: %d, z: %d\n", partitions.size(), partitions[x].size(), partitions[x][y].size());
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

    ParticleSystem(int h, int w, int num_particles) : num_particles(num_particles) {
        pos       = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        deltapos  = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        prevpos   = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        vel       = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        nextvel   = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        vorticity = (glm::vec3 *) malloc(num_particles * sizeof(glm::vec3));
        lambda    = (float *)     calloc(num_particles, sizeof(float));

        neighbors = (int *) malloc(num_particles * MAX_NEIGHBORS * sizeof(int));
        num_neighbors = (int *) calloc(num_particles, sizeof(int));
        
        for (int i = 0; i < num_particles; i++) {
            pos[i] = glm::vec3(0.0f);
            deltapos[i] = glm::vec3(0.0f);
            prevpos[i] = glm::vec3(0.0f);
            vel[i] = glm::vec3(0.0f);
            nextvel[i] = glm::vec3(0.0f);
            vorticity[i] = glm::vec3(0.0f);
        }
        for (int i = 0; i < MAX_NEIGHBORS; i++) {
            for (int j = 0; j < num_particles; j++) {
                neighbors[i*num_particles + j] = 0;
            }
        }

        box = new Box(h, w);
    }

    void init_particle(int id, glm::vec3 p) {
        pos[id] = p;
        prevpos[id] = p;
    }
};