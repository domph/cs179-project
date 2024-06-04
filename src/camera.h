#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GLFW/glfw3.h>

class Camera {
public:
    Camera(glm::vec3 position, glm::vec3 up, float look_x = 0.0f, float look_y = 0.0f, float mouse_sensitivity = 0.05f);
    void process_inputs(GLFWwindow *window, float dt);
    void stop();
    glm::mat4 get_view() { return view; }
    const glm::vec3 &get_position() const { return position; }
    const glm::vec3 &get_view_dir() const { return view_dir; }
    float get_speed() const { return glm::length(velocity); }
private:
    glm::vec3 position;
    glm::vec3 view_dir;
    glm::vec3 up;
    glm::vec3 velocity;
    float mouse_sensitivity;
    glm::mat4 view;
    glm::vec2 look;
    double last_mouse_x;
    double last_mouse_y;
    bool first_mouse;
};

#endif // CAMERA_H
