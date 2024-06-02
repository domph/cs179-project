#pragma once

#include <iostream>
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>

#include "particle.h"


struct Proxy {
    glm::vec3 pos;
    glm::vec3 vel;
};

class ParticleSimulator {
private:
    size_t num_particles;
    float *data;
    GLuint vao;
    GLuint vbo;
    GLsync sync_obj;

    void create_vbo(size_t size);
    void delete_vbo();

public:
    ParticleSimulator();
    // ~ParticleSimulator();
    void render(ParticleSystem *system);
};