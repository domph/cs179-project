#define GLM_ENABLE_EXPERIMENTAL
#include "shaders.h"
#include "physics.h"
#include "particle_simulator.h"


#include <iostream>
#include <algorithm>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

//----------------------------------------
// CONSTANTS
//----------------------------------------

constexpr int START_WIDTH = 1600;
constexpr int START_HEIGHT = 1200;

//----------------------------------------
// GLOBALS
//----------------------------------------

glm::mat4 g_proj_matrix;
ParticleSimulator g_particle_simulator;

//----------------------------------------
// FUNCTIONS
//----------------------------------------

/* Reads in a shader's source from a file and compiles it, returning a new GLenum corresponding to
   the shader. Throws an exception and outputs the error log upon compilation failure. */
GLuint compile_shader(const std::string &source, GLenum shader_type) {
    const char *source_cstr { source.c_str() };

    GLuint new_shader { glCreateShader(shader_type) };
    glShaderSource(new_shader, 1, &source_cstr, NULL);
    glCompileShader(new_shader);

    GLint compiled {};
    glGetShaderiv(new_shader, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE) {
        // Print logs and throw exception
        GLint max_length {};
        glGetShaderiv(new_shader, GL_INFO_LOG_LENGTH, &max_length);

        std::vector<GLchar> error_log(max_length);
        glGetShaderInfoLog(new_shader, max_length, NULL, &error_log[0]);

        for (size_t i {}; i < error_log.size(); ++i) {
            std::cerr << error_log[i];
        }
        glDeleteShader(new_shader);

        throw std::runtime_error{ "Failed to compile shader: " + std::string(source) };
    }

    return new_shader;
}

GLuint compile_program() {
    GLenum shader_program { glCreateProgram() };
    GLenum vertex_shader { compile_shader(VERTEX_SHADER, GL_VERTEX_SHADER) };
    GLenum fragment_shader { compile_shader(FRAGMENT_SHADER, GL_FRAGMENT_SHADER) };

    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    if (glGetError() != GL_NO_ERROR) {
        throw std::runtime_error{ std::string{ "Error linking: " }
            + reinterpret_cast<const char *>(gluErrorString(glGetError())) };
    }

    glLinkProgram(shader_program);
    GLint link_status {};
    glGetProgramiv(shader_program, GL_LINK_STATUS, &link_status);
    if (link_status != GL_TRUE) {
        // Print the log
        GLint max_length {};
        glGetProgramiv(shader_program, GL_INFO_LOG_LENGTH, &max_length);

        std::vector<GLchar> program_log(max_length);
        glGetProgramInfoLog(shader_program, max_length, NULL, &program_log[0]);

        for (size_t i {}; i < program_log.size(); ++i) {
            std::cerr << program_log[i];
        }
        glDeleteProgram(shader_program);
        throw std::runtime_error{ "Failed to link program" };
    }

    return shader_program;
    //glUseProgram(shader_program);
}

// Callback for when the window is resized
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    std::cout << "Window resized to " << width << " x " << height << std::endl;
    glViewport(0, 0, width, height);
    g_proj_matrix = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100.0f);
}
void debug_message_callback(GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar* message,
    const void* userParam)
{
    if (type == GL_DEBUG_TYPE_OTHER || type == GL_DEBUG_TYPE_PERFORMANCE)
        return;

    fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x,\nmessage = %s\n", // NOLINT
        (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
        type, severity, message);
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
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

    // Create a window
    GLFWwindow *window = glfwCreateWindow(START_WIDTH, START_HEIGHT, "CS179 Project", nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Print some specs
    std::cout << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "OpenGL Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    // Window resize callback
    framebuffer_size_callback(window, START_WIDTH, START_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Debug callback
    glDebugMessageCallback(debug_message_callback, nullptr);

    // Render settings
    glPointSize(10.0f);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSwapInterval(1);    // Enable vsync

    // Create shader program
    GLuint shader_program = compile_program();
    glUseProgram(shader_program);

    // Initialize physics
    std::vector<Proxy> particles;
    particles.push_back({ glm::vec3(5.0f, 5.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f) });

    glm::mat4 view_matrix = glm::lookAt(glm::vec3(0.0f, -20.0f, 10.0f), glm::vec3(10.0f, 10.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    glm::mat4 view_proj_matrix = g_proj_matrix * view_matrix;
    glUniformMatrix4fv(2, 1, GL_FALSE, &view_proj_matrix[0][0]);

    // Initialize physics
    float xybound = 10.0f;
    float zbound  = 40.0f;
    float xystep  = 2.0f;
    int klevels = 4;

    int total = 0;
    for (float i = 0; i < xybound; i += xystep) {
        for (float j = 0; j < xybound; j += xystep) {
            for (int k = 0; k < i/klevels; k++) {
                total++;
            }
        }
    }
    printf("total: %d\n", total);

    ParticleSystem *psystem = new ParticleSystem(zbound, xybound, total);
    int id = 0;
    for (float i = 0; i < xybound; i += xystep) {
        for (float j = 0; j < xybound; j += xystep) {
            for (int k = 0; k < i/klevels; k++) {
                psystem->init_particle(id++, glm::vec3(i, j, k*xystep));
            }
        }
    }

    // Main loop
    while(!glfwWindowShouldClose(window)) {
        // Clear; use a black background
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

        // Physics
        update(psystem);
        // std::cout << "particle pos: " << psystem->pos[0].x << ", " << psystem->pos[0].y << ", " << psystem->pos[0].z << std::endl;
        // std::cout << "particle vel: " << psystem->pos[0].x << ", " << psystem->pos[0].y << ", " << psystem->pos[0].z << std::endl;


        
        // shader_program.upload_uniform_mat4f("proj_view_matrix", (view_matrix * g_proj_matrix).data());

        // Rendering
        g_particle_simulator.render(psystem);

        glfwSwapBuffers(window);
        glfwPollEvents();    
    }
    

    // Clean up
    glfwTerminate();
    return 0;
}