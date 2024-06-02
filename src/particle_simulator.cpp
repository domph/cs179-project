#include "particle_simulator.h"

// Struct used internally for data size calculations
struct ParticleData {
    float pos[3];
    float vel[3];
};

void ParticleSimulator::create_vbo(size_t size) {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Create persistent mapped buffer
    const size_t buffer_size = size * sizeof(ParticleData);
    glBufferStorage(GL_ARRAY_BUFFER, buffer_size, 0,
        GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT);
    data = (float*)glMapBufferRange(GL_ARRAY_BUFFER, 0, buffer_size,
        GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT);

    if (data == nullptr) {
        std::cerr << "Failed to map buffer" << std::endl;
        throw std::runtime_error("Failed to map buffer");
    }

    // Set attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleData), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleData), (void*)(offsetof(ParticleData, vel)));
    glEnableVertexAttribArray(1);

    // Unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0); 
    glBindVertexArray(0);
}

// Caller is responsible for checking that vao/vbo exist
void ParticleSimulator::delete_vbo() {
    glBindBuffer(GL_ARRAY_BUFFER, vbo); 
    if (glUnmapBuffer(GL_ARRAY_BUFFER) == GL_FALSE) {
        std::cerr << "Failed to unmap buffer" << std::endl;
        throw std::runtime_error("Failed to unmap buffer");
    }
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    data = nullptr;
}


ParticleSimulator::ParticleSimulator() {
    num_particles = 0;
    data = nullptr;
    vao = 0;
    vbo = 0;
    sync_obj = 0;
}

// ParticleSimulator::~ParticleSimulator() {
//     if (data != nullptr) {
//         delete_vbo();
//     }
//     if (sync_obj) {
//         glDeleteSync(sync_obj);
//     }
// }

void ParticleSimulator::render(ParticleSystem *system) {
    // Wait for GPU to finish using the buffer
    if (sync_obj != 0) {
        GLenum waitReturn = GL_UNSIGNALED;
        while (waitReturn != GL_ALREADY_SIGNALED && waitReturn != GL_CONDITION_SATISFIED) {
            waitReturn = glClientWaitSync(sync_obj, GL_SYNC_FLUSH_COMMANDS_BIT, 1);
        }
    }

    auto &particles = system->particles;

    // Resize if necessary
    if (particles.size() != num_particles) {
        if (data != nullptr) {
            delete_vbo();
        }
        num_particles = particles.size();
        create_vbo(num_particles);
    }

    // Update data
    for (size_t i = 0; i < num_particles; ++i) {
        data[i * 6 + 0] = particles[i]->pos.x;
        data[i * 6 + 1] = particles[i]->pos.y;
        data[i * 6 + 2] = particles[i]->pos.z;
        data[i * 6 + 3] = particles[i]->vel.x;
        data[i * 6 + 4] = particles[i]->vel.y;
        data[i * 6 + 5] = particles[i]->vel.z;
    }

    // Draw
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, num_particles);
    glBindVertexArray(0);

    // Create sync object
    if (sync_obj) {
        glDeleteSync(sync_obj);
    }
    sync_obj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);


}
