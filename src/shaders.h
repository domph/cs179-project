#pragma once
#include <string>

const std::string VERTEX_SHADER = R"(
#version 410 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 velocity;
layout (location = 2) uniform mat4 proj_view_matrix;

out vec3 vertex_color;

void main() {
    // Color particle based on its velocity
    const float red = 0.0045 * dot(velocity, velocity);
    const float green = clamp(0.08 * max(velocity.x, max(velocity.y, velocity.z)), 0.2, 0.5);
    const float blue = 0.7 - red;
    vertex_color = vec3(red, green, blue);

    // Apply projection/view matrices to particle's position
    gl_Position = proj_view_matrix * vec4(position, 1.0);
}
)";

const std::string FRAGMENT_SHADER = R"(
#version 410 core

in vec3 vertex_color;
out vec4 color;  

void main() {
    color = vec4(vertex_color, 1.0);
}
)";