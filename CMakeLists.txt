cmake_minimum_required(VERSION 3.10)
project(cs179-project)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
#set(CMAKE_BUILD_TYPE RELEASE)

set(SRC_LIST
    "src/main.cpp"
    "src/particle_simulator.cpp"
    "src/physics.cpp"
    "src/camera.cpp"
    "src/buffer_manager.cpp"
)

set(IMGUI_SRC_LIST
    "lib/imgui-docking-1.90.7/imgui.cpp"
    "lib/imgui-docking-1.90.7/imgui_draw.cpp"
    "lib/imgui-docking-1.90.7/imgui_impl_glfw.cpp"
    "lib/imgui-docking-1.90.7/imgui_impl_opengl3.cpp"
    "lib/imgui-docking-1.90.7/imgui_stdlib.cpp"
    "lib/imgui-docking-1.90.7/imgui_tables.cpp"
    "lib/imgui-docking-1.90.7/imgui_widgets.cpp"
)

# Source files
add_executable(${PROJECT_NAME} ${SRC_LIST} ${IMGUI_SRC_LIST})

# Add ImGUI include files
include_directories(lib/imgui-docking-1.90.7)

# Add GLM include files
include_directories(lib/glm-1.0.1)

# Add and build GLEW
set(BUILD_UTILS off)
add_subdirectory(lib/glew-2.1.0/build/cmake EXCLUDE_FROM_ALL)
include_directories(lib/glew-2.1.0/include)
target_link_libraries(${PROJECT_NAME} glew_s)

# Add and build GLFW
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_DOCS OFF)
set(GLFW_LIBRARY_TYPE STATIC)
add_subdirectory(lib/glfw-3.4 EXCLUDE_FROM_ALL)
include_directories(lib/glfw-3.4/include)
target_link_libraries(${PROJECT_NAME} glfw)

# Add OpenGL library
find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES})

# Platform-specific settings
if (WIN32)
    target_link_libraries(${PROJECT_NAME} opengl32)
elseif (APPLE)
    # TODO
    #find_library(COCOA_LIBRARY Cocoa)
    #find_library(IOKIT_LIBRARY IOKit)
    #find_library(COREVIDEO_LIBRARY CoreVideo)
    #target_link_libraries(${PROJECT_NAME} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY})
endif()
