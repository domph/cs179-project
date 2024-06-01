#include <box.h>

#include "glm/glm.hpp"

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

struct particle_system_t {
    particle_t *particles;
    int         num_particles;

    particle_t *neighbors;
    int        *num_neighbors;

    box_t      *box;
};