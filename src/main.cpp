#include "shaders.h"

#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <staplegl/staplegl.hpp>

//----------------------------------------
// CONSTANTS
//----------------------------------------

constexpr int START_WIDTH = 1600;
constexpr int START_HEIGHT = 1200;

//----------------------------------------
// GLOBALS
//----------------------------------------

glm::mat4 g_proj_matrix;

//----------------------------------------
// FUNCTIONS
//----------------------------------------

// Callback for when the window is resized
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    std::cout << "Window resized to " << width << " x " << height << std::endl;
    glViewport(0, 0, width, height);
    g_proj_matrix = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100.0f);
}

//----------------------------------------
// MAIN
//----------------------------------------

int main() {
    if (glfwInit() == GLFW_FALSE) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Set the OpenGL version to 4.6
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    // Create a window
    GLFWwindow *window = glfwCreateWindow(START_WIDTH, START_HEIGHT, "CS179 Project", nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Print some specs
    std::cout << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "OpenGL Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    // Window resize callback
    framebuffer_size_callback(window, START_WIDTH, START_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Render settings
    glPointSize(1.1f);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSwapInterval(1);    // Enable vsync

    // Create shader program
    staplegl::shader_program shader_program("main", {
        { staplegl::shader_type::vertex, VERTEX_SHADER },
        { staplegl::shader_type::fragment, FRAGMENT_SHADER }
    });

    // Main loop
    while(!glfwWindowShouldClose(window)) {
        // Physics


        // Rendering


        glfwSwapBuffers(window);
        glfwPollEvents();    
    }

    // Clean up
    glfwTerminate();
    return 0;
}