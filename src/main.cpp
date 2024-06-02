#define GLM_ENABLE_EXPERIMENTAL
#include "shaders.h"
#include "physics.h"
#include "particle_simulator.h"
#include "timer.h"
#include "camera.h"

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

// This is global since the framebuffer_size_callback needs to update it
glm::mat4 g_proj_matrix;

bool g_first_person = false;
bool g_disable_physics = false;
Camera camera(glm::vec3(13.432151, -18.409569, 15.710428), glm::vec3(0, 0, 1), 245.800140, -21.050030);

//----------------------------------------
// FUNCTIONS
//----------------------------------------

/* Reads in a shader's source from a file and compiles it, returning a new GLuint corresponding to
   the shader. Throws an exception and outputs the error log upon compilation failure. */
GLuint compile_shader(const std::string &source, GLenum shader_type) {
    const char *source_cstr = source.c_str();

    GLuint new_shader = glCreateShader(shader_type);
    glShaderSource(new_shader, 1, &source_cstr, NULL);
    glCompileShader(new_shader);

    GLint compiled;
    glGetShaderiv(new_shader, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE) {
        // Print logs and throw exception
        GLint max_length;
        glGetShaderiv(new_shader, GL_INFO_LOG_LENGTH, &max_length);

        std::vector<GLchar> error_log(max_length);
        glGetShaderInfoLog(new_shader, max_length, NULL, &error_log[0]);

        for (size_t i = 0; i < error_log.size(); ++i) {
            std::cerr << error_log[i];
        }
        glDeleteShader(new_shader);

        throw std::runtime_error("Failed to compile shader: " + std::string(source));
    }

    return new_shader;
}

GLuint compile_program() {
    GLuint shader_program = glCreateProgram();
    GLuint vertex_shader = compile_shader(VERTEX_SHADER, GL_VERTEX_SHADER);
    GLuint fragment_shader = compile_shader(FRAGMENT_SHADER, GL_FRAGMENT_SHADER);

    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    if (glGetError() != GL_NO_ERROR) {
        throw std::runtime_error(std::string{ "Error linking: " }
            + reinterpret_cast<const char *>(gluErrorString(glGetError())));
    }

    glLinkProgram(shader_program);
    GLint link_status;
    glGetProgramiv(shader_program, GL_LINK_STATUS, &link_status);
    if (link_status != GL_TRUE) {
        // Print the log
        GLint max_length;
        glGetProgramiv(shader_program, GL_INFO_LOG_LENGTH, &max_length);

        std::vector<GLchar> program_log(max_length);
        glGetProgramInfoLog(shader_program, max_length, NULL, &program_log[0]);

        for (size_t i = 0; i < program_log.size(); ++i) {
            std::cerr << program_log[i];
        }
        glDeleteProgram(shader_program);
        throw std::runtime_error("Failed to link program");
    }

    return shader_program;
}

// Callback for when the window is resized
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    std::cout << "Window resized to " << width << " x " << height << std::endl;
    glViewport(0, 0, width, height);
    g_proj_matrix = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100.0f);
}

void debug_message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length,
    const GLchar* message, const void* userParam) {

    if (type == GL_DEBUG_TYPE_OTHER)
        return;

    std::string source_str;
    switch (source) {
        case GL_DEBUG_SOURCE_API: source_str = "API"; break;
        case GL_DEBUG_SOURCE_WINDOW_SYSTEM: source_str = "Window System"; break;
        case GL_DEBUG_SOURCE_SHADER_COMPILER: source_str = "Shader Compiler"; break;
        case GL_DEBUG_SOURCE_THIRD_PARTY: source_str = "Third Party"; break;
        case GL_DEBUG_SOURCE_APPLICATION: source_str = "Application"; break;
        case GL_DEBUG_SOURCE_OTHER: source_str = "Other"; break;
        default: break;
    }

    std::string type_str;
    switch (type) {
        case GL_DEBUG_TYPE_ERROR: type_str = "Error"; break;
        case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: type_str = "Deprecated Behavior"; break;
        case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: type_str = "Undefined Behavior"; break;
        case GL_DEBUG_TYPE_PORTABILITY: type_str = "Portability"; break;
        case GL_DEBUG_TYPE_PERFORMANCE: type_str = "Performance"; break;
        case GL_DEBUG_TYPE_MARKER: type_str = "Marker"; break;
        case GL_DEBUG_TYPE_PUSH_GROUP: type_str = "Push Group"; break;
        case GL_DEBUG_TYPE_POP_GROUP: type_str = "Pop Group"; break;
        case GL_DEBUG_TYPE_OTHER: type_str = "Other"; break;
        default: break;
    }

    std::string severity_str;
    switch (severity) {
        case GL_DEBUG_SEVERITY_HIGH: severity_str = "High"; break;
        case GL_DEBUG_SEVERITY_MEDIUM: severity_str = "Medium"; break;
        case GL_DEBUG_SEVERITY_LOW: severity_str = "Low"; break;
        case GL_DEBUG_SEVERITY_NOTIFICATION: severity_str = "Notification"; break;
        default: break;
    }

    std::cerr << "[OPENGL EVENT] Source: " << source_str
        << " | Type: " << type_str
        << " | Severity: " << severity_str
        << " | Message: " << message << std::endl;
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_E) {
            g_first_person = !g_first_person;
            if (g_first_person) {
                std::cout << "First person mode enabled" << std::endl;
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            } else {
                std::cout << "First person mode disabled" << std::endl;
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                camera.stop();
            }
        } else if (key == GLFW_KEY_T) {
            g_disable_physics = !g_disable_physics;
            if (g_disable_physics) {
                std::cout << "Physics disabled" << std::endl;
            } else {
                std::cout << "Physics enabled" << std::endl;
            }
        }
    } 
}

void process_input(GLFWwindow *window, double dt) {
    if (g_first_person) {
        camera.process_inputs(window, dt);
    }
    glm::mat4 view_proj_matrix = g_proj_matrix * camera.get_view();
    glUniformMatrix4fv(2, 1, GL_FALSE, &view_proj_matrix[0][0]);
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

    // Set callbacks callback
    framebuffer_size_callback(window, START_WIDTH, START_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glDebugMessageCallback(debug_message_callback, nullptr);

    // Render settings
    glPointSize(5.0f);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSwapInterval(0);    // 0 = disable vsync | 1 = enable vsync

    // Create shader program
    GLuint shader_program = compile_program();
    glUseProgram(shader_program);

    // Initial camera
    glm::mat4 view_proj_matrix = g_proj_matrix * camera.get_view();
    glUniformMatrix4fv(2, 1, GL_FALSE, &view_proj_matrix[0][0]);

    // Initialize physics
    float xybound = 15.0f;
    float zbound  = 40.0f;
    float xystep  = 0.5f;
    int klevels = 4;

    int total = 0;
    for (float i = 0; i < xybound; i += xystep) {
        for (float j = 0; j < xybound; j += xystep) {
            for (int k = 0; k < i/klevels; k++) {
                total++;
            }
        }
    }
    printf("total particles: %d\n", total);

    ParticleSystem *psystem = new ParticleSystem(zbound, xybound, total);
    int id = 0;
    for (float i = 0; i < xybound; i += xystep) {
        for (float j = 0; j < xybound; j += xystep) {
            for (int k = 0; k < i/klevels; k++) {
                psystem->init_particle(id++, glm::vec3(i, j, k*xystep));
            }
        }
    }

    ParticleSimulator particle_simulator;
    FPSTimer fps_timer;
    Timer frame_timer;

    // Main loop
    while(!glfwWindowShouldClose(window)) {
        // Update FPS
        glfwSetWindowTitle(window, ("CS179 Project | FPS: " + std::to_string(fps_timer.get_fps())).c_str());

        // Clear screen; use a black background
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

        // FPS
        fps_timer.update();

        // User input
        double dt = frame_timer.elapsed_s();
        frame_timer.reset();
        process_input(window, dt);

        // Physics
        // Timer physics_timer;
        if (!g_disable_physics) {
            update(psystem);
        }
        
        // std::cout << "physics_time: " << physics_timer.elapsed_s() << " s" << std::endl;

        // Rendering
        Timer render_timer;
        particle_simulator.render(psystem);
        
        // Swap buffers and poll events
        glfwSwapBuffers(window);
        // std::cout << "render_time: " << render_timer.elapsed_s() << " s" << std::endl;
        glfwPollEvents();    
    }

    // Clean up
    glfwTerminate();
    return 0;
}