#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "constants.h"

struct Particle {
    int id;
    bool fixed;
    
    glm::vec3 pos;
    glm::vec3 deltapos;
    glm::vec3 prevpos;
    
    glm::vec3 vel;
    glm::vec3 nextvel;

    glm::vec3 vorticity;
    float lambda;

    std::vector<Particle *> neighbors;

    Particle (int id, bool fixed, glm::vec3 pos): id(id), fixed(fixed), pos(pos) {
        deltapos = glm::vec3(0.0f);
        prevpos = pos;
        vel = glm::vec3(0.0f);
        nextvel = glm::vec3(0.0f);
        vorticity = glm::vec3(0.0f);
        lambda = 0.0f;

        neighbors.reserve(MAX_NEIGHBORS);
    }
};


/* Represents a box with a w x w square base and height h */
struct Box {
    int h;
    int w;
    std::vector<std::vector<std::vector<std::vector<Particle *>>>> partitions;

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

    void add_particle(Particle *particle) {
        int x = (particle->pos.x + EPS) / P_H;
        int y = (particle->pos.y + EPS) / P_H;
        int z = (particle->pos.z + EPS) / P_H;

        if (x >= 0 && x < partitions.size() &&
            y >= 0 && y < partitions[x].size() &&
            z >= 0 && z < partitions[x][y].size()) {
            partitions[x][y][z].push_back(particle);
        } else {
            printf("Box error: particle out of bounds!\n");
            printf("loc:   x: %d, y: %d, z: %d\n", x, y, z);
            printf("bound: x: %d, y: %d, z: %d\n", partitions.size(), partitions[x].size(), partitions[x][y].size());
        }
    }

    void clear_partitions() {
        for (auto x_layer : partitions) {
            for (auto y_layer : x_layer) {
                for (auto z_layer : y_layer) {
                    z_layer.clear();
                }
            }
        }
    }

};

struct ParticleSystem {
    std::vector<Particle *> particles;
    Box *box;

    ParticleSystem(int h, int w) {
        box = new Box(h, w);
    }

    void add_particle(bool fixed, glm::vec3 pos) {
        Particle *particle = new Particle(particles.size(), fixed, pos);
        particles.push_back(particle);
    }
};