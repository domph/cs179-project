#define GLM_ENABLE_EXPERIMENTAL
#define GL_SILENCE_DEPRECATION
#include "shaders.h"
#include "physics.h"
#include "particle_simulator.h"
#include "timer.h"
#include "camera.h"
#include "fonts.h"

#include <iostream>
#include <algorithm>
#include <format>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGl/gl3ext.h>
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include <imgui.h>
#include <imgui_internal.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

//----------------------------------------
// CONSTANTS
//----------------------------------------

constexpr int START_WIDTH = 1800;
constexpr int START_HEIGHT = 1300;

//----------------------------------------
// GLOBALS
//----------------------------------------

ParticleSystem *g_psystem;
ParticleSimulator *g_particle_simulator;
glm::mat4 g_proj_matrix;
GLint g_proj_matrix_loc;
Timer g_timer;

bool g_enable_camera = false;
bool g_enable_physics = true;
bool g_shake = false;
bool g_vsync = true;
float g_point_size = 4.0f;
Camera g_camera(glm::vec3(-15, -15, 15), glm::vec3(0, 0, 1), 315, -12);

constexpr float BASE_FONT_SIZE = 15.0f;

//----------------------------------------
// FUNCTIONS
//----------------------------------------

// https://github.com/glfw/glfw/issues/1699#issuecomment-723692566
bool glfw_get_window_monitor(GLFWmonitor **monitor, GLFWwindow *window) {
    bool success = false;

    int window_rectangle[4] = {0};
    glfwGetWindowPos(window, &window_rectangle[0], &window_rectangle[1]);
    glfwGetWindowSize(window, &window_rectangle[2], &window_rectangle[3]);

    int monitors_size = 0;
    GLFWmonitor **monitors = glfwGetMonitors(&monitors_size);

    GLFWmonitor *closest_monitor = NULL;
    int max_overlap_area = 0;

    for (int i = 0; i < monitors_size; ++i) {
        int monitor_position[2] = {0};
        glfwGetMonitorPos(monitors[i], &monitor_position[0], &monitor_position[1]);

        const GLFWvidmode *monitor_video_mode = glfwGetVideoMode(monitors[i]);

        int monitor_rectangle[4] = {
            monitor_position[0],
            monitor_position[1],
            monitor_video_mode->width,
            monitor_video_mode->height,
        };

        if (
            !(
                ((window_rectangle[0] + window_rectangle[2]) < monitor_rectangle[0]) ||
                (window_rectangle[0] > (monitor_rectangle[0] + monitor_rectangle[2])) ||
                ((window_rectangle[1] + window_rectangle[3]) < monitor_rectangle[1]) ||
                (window_rectangle[1] > (monitor_rectangle[1] + monitor_rectangle[3]))
            )
        ) {
            int intersection_rectangle[4] = {0};

            // x, width
            if (window_rectangle[0] < monitor_rectangle[0]) {
                intersection_rectangle[0] = monitor_rectangle[0];

                if ((window_rectangle[0] + window_rectangle[2]) < (monitor_rectangle[0] + monitor_rectangle[2])) {
                    intersection_rectangle[2] = (window_rectangle[0] + window_rectangle[2]) - intersection_rectangle[0];
                }
                else {
                    intersection_rectangle[2] = monitor_rectangle[2];
                }
            }
            else {
                intersection_rectangle[0] = window_rectangle[0];

                if ((monitor_rectangle[0] + monitor_rectangle[2]) < (window_rectangle[0] + window_rectangle[2])) {
                    intersection_rectangle[2] = (monitor_rectangle[0] + monitor_rectangle[2]) - intersection_rectangle[0];
                }
                else {
                    intersection_rectangle[2] = window_rectangle[2];
                }
            }

            // y, height
            if (window_rectangle[1] < monitor_rectangle[1]) {
                intersection_rectangle[1] = monitor_rectangle[1];

                if ((window_rectangle[1] + window_rectangle[3]) < (monitor_rectangle[1] + monitor_rectangle[3])) {
                    intersection_rectangle[3] = (window_rectangle[1] + window_rectangle[3]) - intersection_rectangle[1];
                }
                else {
                    intersection_rectangle[3] = monitor_rectangle[3];
                }
            }
            else {
                intersection_rectangle[1] = window_rectangle[1];

                if ((monitor_rectangle[1] + monitor_rectangle[3]) < (window_rectangle[1] + window_rectangle[3])) {
                    intersection_rectangle[3] = (monitor_rectangle[1] + monitor_rectangle[3]) - intersection_rectangle[1];
                }
                else {
                    intersection_rectangle[3] = window_rectangle[3];
                }
            }

            int overlap_area = intersection_rectangle[2] * intersection_rectangle[3];
            if (overlap_area > max_overlap_area) {
                closest_monitor = monitors[i];
                max_overlap_area = overlap_area;
            }
        }
    }

    if (closest_monitor) {
        *monitor = closest_monitor;
        success = true;
    }

    // true: monitor contains the monitor the window is most on
    // false: monitor is unmodified
    return success;
}

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
#ifndef __APPLE__
            + reinterpret_cast<const char *>(gluErrorString(glGetError())));
#else
            + std::to_string(glGetError()));
#endif
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

#ifndef __APPLE__
void debug_message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length,
    const GLchar* message, const void* userParam) {
    (void)id;
    (void)length;
    (void)userParam;

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
#endif

void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    (void)window;

    glViewport(0, 0, width, height);
}  

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    (void)mods;
    (void)scancode;

    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_E) {
            g_enable_camera = !g_enable_camera;
            if (g_enable_camera) {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            } else {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                g_camera.stop();
            }
        }
    } 
}

void process_input(GLFWwindow *window, double dt) {
    if (g_enable_camera) {
        g_camera.process_inputs(window, dt);
    }
    glm::mat4 view_proj_matrix = g_proj_matrix * g_camera.get_view();
    glUniformMatrix4fv(g_proj_matrix_loc, 1, GL_FALSE, &view_proj_matrix[0][0]);
}

ImGuiStyle new_imgui_style() {
    ImGuiStyle style = ImGuiStyle();
    ImGui::StyleColorsDark(&style);

    auto &colors = style.Colors;
    colors[ImGuiCol_WindowBg] = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);

    // colors[ImGuiCol_Header] = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);
    // colors[ImGuiCol_HeaderHovered] = ImVec4(0.3f, 0.3f, 0.3f, 1.0f);
    // colors[ImGuiCol_HeaderActive] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);

    // colors[ImGuiCol_Button] = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);
    // colors[ImGuiCol_ButtonHovered] = ImVec4(0.3f, 0.3f, 0.3f, 1.0f);
    // colors[ImGuiCol_ButtonActive] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);

    // colors[ImGuiCol_FrameBg] = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);
    // colors[ImGuiCol_FrameBgHovered] = ImVec4(0.3f, 0.3f, 0.3f, 1.0f);
    // colors[ImGuiCol_FrameBgActive] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);

    // colors[ImGuiCol_Tab] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);
    // colors[ImGuiCol_TabHovered] = ImVec4(0.38f, 0.38f, 0.38f, 1.0f);
    // colors[ImGuiCol_TabActive] = ImVec4(0.28f, 0.28f, 0.28f, 1.0f);
    // colors[ImGuiCol_TabUnfocused] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);
    // colors[ImGuiCol_TabUnfocusedActive] = ImVec4(0.2f, 0.2f, 0.2f, 1.0f);

    colors[ImGuiCol_TitleBg] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);
    colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.15f, 0.15f, 0.15f, 1.0f);

    style.WindowRounding = 0.0f;
    style.TabRounding = 0.0f;
    return style;
}

void build_imgui_dock() {
    static bool init = true;

    ImGuiID dockspace_id = ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport(), ImGuiDockNodeFlags_PassthruCentralNode);
    ImGuiID dock_id_left, dock_id_right;
    if (init) {
        init = false;
        ImGui::DockBuilderRemoveNode(dockspace_id);
        ImGui::DockBuilderAddNode(dockspace_id);
        ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->Size);

        ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.25f, &dock_id_left, &dock_id_right);
        ImGui::DockBuilderDockWindow("Control Panel", dock_id_left);
        ImGui::DockBuilderDockWindow("Scene", dock_id_right);

        ImGui::DockBuilderFinish(dockspace_id);
    }
}

// ImGui helper functions
namespace ImGui {
    void TextWithBackground(const char *text, bool offset_half_char = false) {
        ImDrawList *draw_list = ImGui::GetWindowDrawList();
        ImGuiStyle &style = ImGui::GetStyle();

        ImVec2 cursor_pos = ImGui::GetCursorScreenPos();
        ImVec2 text_pos(cursor_pos.x + style.FramePadding.x, cursor_pos.y + style.FramePadding.y);
        if (offset_half_char) {
            text_pos.x += ImGui::CalcTextSize(" ").x / 2;
        }

        ImVec2 text_size = ImGui::CalcTextSize(text);
        ImVec2 rect_size(text_size.x + style.FramePadding.x * 2, text_size.y + style.FramePadding.y * 2);

        ImVec2 rect_min = cursor_pos;
        ImVec2 rect_max(rect_min.x + rect_size.x, rect_min.y + rect_size.y);

        ImU32 rect_color = ImGui::GetColorU32(ImGui::GetColorU32(style.Colors[ImGuiCol_Button]), style.DisabledAlpha);
        ImU32 text_color = ImGui::GetColorU32(style.Colors[ImGuiCol_Text]);

        draw_list->AddRectFilled(rect_min, rect_max, rect_color);
        draw_list->AddText(text_pos, text_color, text);

        ImGui::Dummy(rect_size);
    }
    void OnOff(const char *label, bool val) {
        ImGui::Text("%s", label);
        ImGui::SameLine();
        if (val) {
            ImGui::TextColored(ImVec4(0, 1, 0, 1), "ON");
        } else {
            ImGui::TextColored(ImVec4(1, 0, 0, 1), "OFF");
        }
    }
}

void build_control_panel() {
    ImGui::Begin("Control Panel");
    ImGui::TextWrapped("This is a simulation of fluid particles in an invisible container. The color of the particles is determined by their speed; slower particles are bluer, and faster particles are redder. The position and direction of the camera can be adjusted via the controls below.");
    if (ImGui::CollapsingHeader("Statistics", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("Particles: %zu", g_psystem->num_particles);
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::Text("Global time: %.2f s", g_timer.elapsed_s());
        ImGui::Text("Physics time: %.2f s", g_psystem->t);
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Settings", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Physics", &g_enable_physics);
        ImGui::Checkbox("Shaking", &g_shake);
        
#ifndef __APPLE__
        if (ImGui::Checkbox("VSync", &g_vsync)) {
            glfwSwapInterval(g_vsync ? 1 : 0);
        }
#endif
        if (ImGui::SliderFloat("Particle size", &g_point_size, 0.1f, 10.0f)) {
            glPointSize(g_point_size);
        }

        if (ImGui::Button("Reset Container", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
            std::cout << "Resetting container" << std::endl;
            g_psystem->respawn();
        }
        ImGui::Spacing();
    }
       
    if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::OnOff("Camera control:", g_enable_camera);
        ImGui::Spacing();

        ImGui::SeparatorText("Controls");
        ImGui::TextWithBackground("     E     ");
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Toggle camera control");

        ImGui::TextWithBackground("   Mouse   ");
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Change camera direction");

        ImGui::TextWithBackground("   WASD    ", true);
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Change camera position");

        ImGui::TextWithBackground("Hold LShift");
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Faster camera movement");

        ImGui::TextWithBackground("Hold LCtrl ", true);
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Slower camera movement");

        ImGui::Spacing();
        ImGui::SeparatorText("Info");

        ImGui::Text("Camera position: (%.2f, %.2f, %.2f)", g_camera.get_position().x, g_camera.get_position().y, g_camera.get_position().z);
        ImGui::Text("Camera direction: (%.2f, %.2f, %.2f)", g_camera.get_view_dir().x, g_camera.get_view_dir().y, g_camera.get_view_dir().z);
        ImGui::Text("Camera speed: %.2f", g_camera.get_speed());
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Add Particles", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::TextWrapped("This allows you to add parcels of particles of a desired size to the simulation at a desired location and velocity in the box.");
    
        static float pos_x = PARCEL_DEFAULT_XY;
        static float pos_y = PARCEL_DEFAULT_XY;
        static float pos_z = PARCEL_DEFAULT_Z;
        static float vel_z = PARCEL_DEFAULT_Z_VEL;
        static float r     = PARCEL_R_DEFAULT;

        pos_x = glm::max(pos_x, PARCEL_MIN_XY);
        pos_x = glm::min(pos_x, PARCEL_MAX_XY);
        pos_y = glm::max(pos_y, PARCEL_MIN_XY);
        pos_y = glm::min(pos_y, PARCEL_MAX_XY);
        pos_z = glm::max(pos_z, PARCEL_MIN_Z);
        pos_z = glm::min(pos_z, PARCEL_MAX_Z);

        ImGui::SeparatorText("Radius");
        ImGui::SliderFloat("", &r, PARCEL_R_MIN, PARCEL_R_MAX);
        ImGui::SeparatorText("Position");
        ImGui::SliderFloat("X", &pos_x, PARCEL_MIN_XY, PARCEL_MAX_XY);
        ImGui::SliderFloat("Y", &pos_y, PARCEL_MIN_XY, PARCEL_MAX_XY);
        ImGui::SliderFloat("Z", &pos_z, PARCEL_MIN_Z, PARCEL_MAX_Z);
        ImGui::SeparatorText("Velocity");
        ImGui::SliderFloat("Z ", &vel_z, -PARCEL_MIN_Z_VEL, -PARCEL_MAX_Z_VEL);\
        ImGui::Spacing();
        // ImGui::SeparatorText("");
        ImGui::Spacing();

        if (ImGui::Button("Add Parcel", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
            // std::cout << "Adding parcel to pos: " << pos_x << ", " << pos_y << std::endl;
            g_psystem->spawn_parcel(pos_x, pos_y, pos_z, vel_z, r);
        }
    }
    ImGui::End();
}

void build_scene() {
    ImGui::Begin("Scene");

    ImVec2 size = ImGui::GetContentRegionAvail();
    if (size.x <= 0 || size.y <= 0) {
        ImGui::End();
        return;
    }

    g_proj_matrix = glm::perspective(glm::radians(45.0f), (float)size.x / (float)size.y, 0.1f, 100.0f);
    if (g_particle_simulator == nullptr) {
        g_particle_simulator = new ParticleSimulator(size.x, size.y);
    }
    g_particle_simulator->update_viewport(size.x, size.y);

    // Render
    ImGui::Image(
        (void*)(int64_t)g_particle_simulator->render(g_psystem),
        size,
        ImVec2(0, 1), 
        ImVec2(1, 0)
    );
    ImGui::End();
}

void check_dpi(GLFWwindow *window) {
    static float dpi_scale = 1.0;
    GLFWmonitor *monitor;
    if (glfw_get_window_monitor(&monitor, window)) {
        float xscale, yscale;
        glfwGetMonitorContentScale(monitor, &xscale, &yscale);

        // xscale should be equal to yscale
        if (xscale != dpi_scale) {  
            dpi_scale = xscale;
            std::cout << "Changing DPI scale to: " << dpi_scale << std::endl;

            // Clear existing fonts
            ImGuiIO &io = ImGui::GetIO();

            io.Fonts->Clear();

            float font_size = BASE_FONT_SIZE * xscale;
            ImFontConfig font_cfg;
            font_cfg.FontDataOwnedByAtlas = false;
            io.Fonts->AddFontFromMemoryTTF(proggyvector_regular_ttf, proggyvector_regular_ttf_len, font_size, &font_cfg);
            
            ImGui_ImplOpenGL3_CreateFontsTexture();

            ImGuiStyle style = new_imgui_style();
            style.ScaleAllSizes(xscale);
            ImGui::GetStyle() = style;
        }
    }
}

//----------------------------------------
// MAIN
//----------------------------------------

int main() {
    if (glfwInit() == GLFW_FALSE) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Set the OpenGL version to 4.1
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
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

#ifndef __APPLE__
    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }
#endif

    // Print some specs
    std::cout << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "OpenGL Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    // Set callbacks callback
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);  
    glfwSetKeyCallback(window, key_callback);
#ifndef __APPLE__
    glDebugMessageCallback(debug_message_callback, nullptr);
#endif
    glViewport(0, 0, START_WIDTH, START_HEIGHT);

    // Initialize IMGUI
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    ImGui_ImplGlfw_InitForOpenGL(window, true);          // install_callback=true will install GLFW callbacks and chain to existing ones
    ImGui_ImplOpenGL3_Init();

#ifdef __APPLE__
// https://github.com/ocornut/imgui/issues/5301#issuecomment-1122067363
    io.Fonts->Clear();

    ImFontConfig font_cfg;
    font_cfg.FontDataOwnedByAtlas = false;
    io.Fonts->AddFontFromMemoryTTF(proggyvector_regular_ttf, proggyvector_regular_ttf_len, BASE_FONT_SIZE * 2, &font_cfg);
    io.FontGlobalScale = 0.5f;
    ImGui_ImplOpenGL3_CreateFontsTexture();
#endif
    ImGui::GetStyle() = new_imgui_style();

    // Render settings
    glPointSize(g_point_size);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSwapInterval(g_vsync ? 1 : 0);
    if (g_enable_camera) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    } else {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        g_camera.stop();
    }

    // Create shader program
    GLuint shader_program = compile_program();
    glUseProgram(shader_program);

    // Initial camera
    glm::mat4 view_proj_matrix = g_proj_matrix * g_camera.get_view();
    g_proj_matrix_loc = glGetUniformLocation(shader_program, "proj_view_matrix");
    glUniformMatrix4fv(g_proj_matrix_loc, 1, GL_FALSE, &view_proj_matrix[0][0]);

    // Initialize physics
    g_psystem = new ParticleSystem();

    Timer frame_timer;

    // Main loop
    while(!glfwWindowShouldClose(window)) {
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT);
        
        glfwSetWindowTitle(window, "CS179 Project");
#ifndef __APPLE__
        check_dpi(window);  // must be called before imgui::newframe()
#endif
        // User input
        double dt = frame_timer.elapsed_s();
        frame_timer.reset();
        process_input(window, dt);

        // Physics
        if (g_enable_physics) {
            update(g_psystem, g_shake);
        }
        
        // Imgui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        build_imgui_dock();
        build_control_panel();
        build_scene();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();    
    }

    // Clean up
    delete g_psystem;
    delete g_particle_simulator;

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}