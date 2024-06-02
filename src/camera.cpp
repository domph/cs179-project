#define GLM_ENABLE_EXPERIMENTAL
#include "Camera.h"
#include <cmath>
#include <iostream>
#include <glm/gtx/string_cast.hpp>

Camera::Camera(glm::vec3 position, glm::vec3 up, float look_x, float look_y, float mouse_sensitivity, float speed)
    : position(position), up(up), mouse_sensitivity(mouse_sensitivity), movement_speed(speed), velocity(glm::vec3(0.0f)) {
    look.x = look_x;
    look.y = look_y;

    view_dir.x = cos(glm::radians(look.x)) * cos(glm::radians(look.y));
    view_dir.z = sin(glm::radians(look.y));
    view_dir.y = -sin(glm::radians(look.x)) * cos(glm::radians(look.y));

    view = glm::lookAt(position, position + view_dir, up);
    first_mouse = true;
}

void Camera::stop() {
    velocity = glm::vec3(0.0f);
    first_mouse = true;
}

void Camera::process_inputs(GLFWwindow *window, float dt) {
    double mouse_x, mouse_y;
    glfwGetCursorPos(window, &mouse_x, &mouse_y);

    if (first_mouse) {
        last_mouse_x = mouse_x;
        last_mouse_y = mouse_y;
        first_mouse = false;
        return;
    }

    double delta_x = mouse_x - last_mouse_x;
    double delta_y = last_mouse_y - mouse_y;

    last_mouse_x = mouse_x;
    last_mouse_y = mouse_y;

    look.x += delta_x * mouse_sensitivity;
    look.y += delta_y * mouse_sensitivity;

    if (look.y >= 90.0f) look.y = 89.999f;
    if (look.y <= -90.0f) look.y = -89.999f;

    view_dir.x = cos(glm::radians(look.x)) * cos(glm::radians(look.y));
    view_dir.z = sin(glm::radians(look.y));
    view_dir.y = -sin(glm::radians(look.x)) * cos(glm::radians(look.y));

    glm::vec3 acceleration(0.0f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        acceleration += view_dir;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        acceleration -= view_dir;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        acceleration += glm::normalize(glm::cross(view_dir, up));
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        acceleration -= glm::normalize(glm::cross(view_dir, up));

    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        velocity += acceleration * 5.0f;
    else if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
        velocity += acceleration * 0.35f;
    else
        velocity += acceleration;

    if (glm::length(velocity) < 0.01f)
        velocity = glm::vec3(0.0f);

    position += velocity * dt;
    velocity *= 0.95f;

    // std::cout << "position: " << glm::to_string(position) << std::endl;
    // std::cout << "look: " << glm::to_string(look) << std::endl;
    view = glm::lookAt(position, position + view_dir, up);
}