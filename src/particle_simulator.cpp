#define GL_SILENCE_DEPRECATION
#include "particle_simulator.h"
#include "physics.cuh"

void ParticleSimulator::create_vbo(size_t size) {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Create persistent mapped buffer
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(float) * 6, 0, GL_DYNAMIC_DRAW);
    buffer = new float[size * 6];
    // glBufferStorage(GL_ARRAY_BUFFER, buffer_size, 0,
    //     GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT);
    // buffer = (float*)glMapBufferRange(GL_ARRAY_BUFFER, 0, buffer_size,
    //     GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT);

    // if (buffer == nullptr) {
    //     std::cerr << "Failed to map buffer" << std::endl;
    //     throw std::runtime_error("Failed to map buffer");
    // }

    // Set attributes
    // Position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
    glEnableVertexAttribArray(0);

    // Velocity
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)(sizeof(float) * 3 * size));
    glEnableVertexAttribArray(1);

    // Unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0); 
    glBindVertexArray(0);
}

// Caller is responsible for checking that vao/vbo exist
void ParticleSimulator::delete_vbo() {
    // glBindBuffer(GL_ARRAY_BUFFER, vbo); 
    // if (glUnmapBuffer(GL_ARRAY_BUFFER) == GL_FALSE) {
    //     std::cerr << "Failed to unmap buffer" << std::endl;
    //     throw std::runtime_error("Failed to unmap buffer");
    // }
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    delete buffer;
}

ParticleSimulator::~ParticleSimulator() {
    // The only time we don't delete the buffer/delete the sync object is when no render()
    // call has ever been made, which we can check for by simply seeing if vbo is 0
    if (vbo) {
        delete_vbo();
        glDeleteSync(sync_obj);
    }
}

void ParticleSimulator::check_size(size_t new_size) {
    if (new_size != num_particles) {
        if (vbo) {
            delete_vbo();
        }
        create_vbo(new_size);
        num_particles = new_size;
    }
}

GLuint ParticleSimulator::render(ParticleSystem *system, ParticleSystem *gpu_system, bool using_gpu) {
    // Wait for GPU to finish using the buffer
    if (sync_obj) {
        GLenum waitReturn = GL_UNSIGNALED;
        while (waitReturn != GL_ALREADY_SIGNALED && waitReturn != GL_CONDITION_SATISFIED) {
            waitReturn = glClientWaitSync(sync_obj, GL_SYNC_FLUSH_COMMANDS_BIT, 1);
        }
    }

    // Resize texture/RBO if necessary
    buffer_manager.rescale(viewport_width, viewport_height);

    // Resize VBO if necessary
    check_size(system->num_particles);

    // Update data
    if (using_gpu) {
        cudaCopyToRenderBuffer(system, gpu_system, (float*)buffer);
    } else {
        for (size_t i = 0; i < num_particles; ++i) {
            buffer[i * 3 + 0] = system->pos[i].x;
            buffer[i * 3 + 1] = system->pos[i].y;
            buffer[i * 3 + 2] = system->pos[i].z;
        }

        float *vel_buffer = buffer + num_particles * 3;
        for (size_t i = 0; i < num_particles; ++i) {
            vel_buffer[i * 3 + 0] = system->vel[i].x;
            vel_buffer[i * 3 + 1] = system->vel[i].y;
            vel_buffer[i * 3 + 2] = system->vel[i].z;
        }
    }
    
    // Update data
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, num_particles * sizeof(float) * 6, buffer);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Draw
    buffer_manager.bind();
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, (GLsizei)num_particles);
    glBindVertexArray(0);
    buffer_manager.unbind();

    // Create sync object
    if (sync_obj) {
        glDeleteSync(sync_obj);
    }
    sync_obj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);

    return buffer_manager.get_texture_id();
}