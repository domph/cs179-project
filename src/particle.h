#pragma once

#include <glm/glm.hpp>
#include <vector>
#include "constants.h"

struct particle_t {
    int id;
    bool fixed;
    
    glm::vec3 pos;
    glm::vec3 deltapos;
    glm::vec3 prevpos;
    
    glm::vec3 vel;
    glm::vec3 nextvel;

    glm::vec3 vorticity;

    float lambda;

    particle_t (bool fixed = false): fixed(fixed) {
        id = 0;
        pos = glm::vec3(0.0f);
        deltapos = glm::vec3(0.0f);
        prevpos = glm::vec3(0.0f);
        vel = glm::vec3(0.0f);
        nextvel = glm::vec3(0.0f);
        vorticity = glm::vec3(0.0f);
        lambda = 0.0f;
    }
};


/* Represents a box with a w x w square base and height h */
struct box_t {
    int h;
    int w;
    std::vector<std::vector<std::vector<std::vector<particle_t *>>>> partitions;

    int x_partitions;
    int y_partitions;
    int z_partitions;
    int total_partitions;

    box_t(int w, int h) : w(w), h(h) {
        x_partitions = (float)w / H + 1;
        y_partitions = (float)w / H + 1;
        z_partitions = (float)h / H + 1;
        total_partitions = x_partitions * y_partitions * z_partitions;

        partitions.resize(x_partitions);
        for (int x = 0; x < x_partitions; x++) {
            partitions[x].resize(y_partitions);
            for (int y = 0; y < y_partitions; y++) {
                partitions[x][y].resize(z_partitions);
            }
        }
    }

    void add_particle(particle_t *particle) {
        int x = (particle->pos.x + EPS) / H;
        int y = (particle->pos.y + EPS) / H;
        int z = (particle->pos.z + EPS) / H;

        if (x >= 0 && x < partitions.size() &&
            y >= 0 && y < partitions[x].size() &&
            z >= 0 && z < partitions[x][y].size()) {
            partitions[x][y][z].push_back(particle);
        }
    }

    void clear_partitions() {
        for (auto& x_layer : partitions) {
            for (auto& y_layer : x_layer) {
                for (auto& z_layer : y_layer) {
                    z_layer.clear();
                }
            }
        }
    }

};

struct particle_system_t {
    particle_t *particles;
    int         num_particles;

    particle_t *neighbors;
    int        *num_neighbors;

    box_t      *box;
};