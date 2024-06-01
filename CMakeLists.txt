cmake_minimum_required(VERSION 3.10)
project(cs179-project)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE RELEASE)

# Source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})
include_directories(src)

# Add GLM
include_directories(lib/glm-1.0.1)

# Add and build freeglut
# set(FREEGLUT_BUILD_DEMOS OFF)
# set(FREEGLUT_BUILD_SHARED_LIBS OFF)
# add_subdirectory(lib/freeglut-3.4.0 EXCLUDE_FROM_ALL)
# include_directories(lib/freeglut-3.4.0/include)
# target_link_libraries(${PROJECT_NAME} freeglut_static)

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
