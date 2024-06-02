#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include "particle.h"

/* Handles rendering of particles */
class ParticleSimulator {
private:
    size_t num_particles = 0;
    float *buffer        = nullptr;
    GLuint vao           = 0;
    GLuint vbo           = 0;
    GLsync sync_obj      = 0;

    void check_size(size_t new_size);
    void create_vbo(size_t size);
    // Set check_gl_errors to false to avoid checking for GL errors; used in destructor
    // since it's possible that the context has already been destroyed
    void delete_vbo(bool check_gl_errors = true);
public:
    ~ParticleSimulator();
    void render(ParticleSystem *system);
};