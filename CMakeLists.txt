cmake_minimum_required(VERSION 3.20)
if (NOT APPLE)
    project(cs179-project LANGUAGES C CXX CUDA)
else()
    project(cs179-project LANGUAGES C CXX)
endif()

##### CONFIG #####
set(ENABLE_DEBUG_MODE OFF)       # ON = Includes debugging symbols, but less optimizations | OFF = No debug symbols, but more optimizations 
set(ENABLE_STRICT_MODE OFF)      # ON = Shows more warnings/errors | OFF = Shows less warnings/errors
##################


# Set CUDA architectures
set(CMAKE_CUDA_ARCHITECTURES 75 86 89)

# Statically link runtime libraries for MSVC
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Debug/Release
if (ENABLE_DEBUG_MODE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
else()
    set(CMAKE_BUILD_TYPE Release)
endif() 

set(SRC_LIST
    "src/main.cpp"
    "src/particle_simulator.cpp"
    "src/physics.cpp"
    "src/camera.cpp"
    "src/buffer_manager.cpp"
)

set(CUDA_SRC_LIST
    "src/physics.cu"
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
if (NOT APPLE)
    add_executable(${PROJECT_NAME} ${SRC_LIST} ${CUDA_SRC_LIST})
else()
    add_executable(${PROJECT_NAME} ${SRC_LIST})
endif()

# Compiler flags for strict mode
if (ENABLE_STRICT_MODE)
    if (MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/W4 /WX>
            #$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4 -Xcompiler=/WX -dc -diag-suppress 20012>
        )
    else()
        target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if (NOT APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress 20012>
    )
endif()

# CUDA properties
set_target_properties(${PROJECT_NAME}
   PROPERTIES
       CUDA_SEPARABLE_COMPILATION ON
       CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# Add CUDA include files
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Add GLM include files
include_directories(lib/glm-1.0.1)

# Add and build ImGui
add_library(imgui ${IMGUI_SRC_LIST})
include_directories(lib/imgui-docking-1.90.7)
target_link_libraries(${PROJECT_NAME} imgui)

# Add and build GLEW (not needed on macOS)
if (NOT APPLE)
    set(BUILD_UTILS off)
    add_subdirectory(lib/glew-2.1.0/build/cmake EXCLUDE_FROM_ALL)
    include_directories(lib/glew-2.1.0/include)
    target_link_libraries(${PROJECT_NAME} glew_s)
endif()

# Add and build GLFW
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_DOCS OFF)
set(GLFW_LIBRARY_TYPE STATIC)
add_subdirectory(lib/glfw-3.4 EXCLUDE_FROM_ALL)
include_directories(lib/glfw-3.4/include)
target_link_libraries(${PROJECT_NAME} glfw)

# Add OpenGL
find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES})