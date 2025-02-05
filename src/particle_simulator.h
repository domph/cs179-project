#pragma once

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGl/gl3ext.h>
#else
#include <GL/glew.h>
#endif
#include <glm/glm.hpp>
#include "particle.h"
#include "buffer_manager.h"

struct ParticleData {
    float pos[3];
    float vel[3];
};

/* Handles rendering of particles */
class ParticleSimulator {
private:
    size_t num_particles = 0;
    ParticleData *buffer = nullptr;
    GLuint vao           = 0;
    GLuint vbo           = 0;
    GLsync sync_obj      = 0;
    int viewport_width   = 0;
    int viewport_height  = 0;

    BufferManager buffer_manager;

    void check_size(size_t new_size);
    void create_vbo(size_t size);
    void delete_vbo();
public:
    ParticleSimulator(int viewport_width, int viewport_height) : buffer_manager(viewport_width, viewport_height) {}
    ~ParticleSimulator();
    GLuint render(ParticleSystem *system);
    void update_viewport(int width, int height) {
        viewport_width = width;
        viewport_height = height;
    }
};