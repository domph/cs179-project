#pragma once

#include <GL/glew.h>

class BufferManager {
private:
    GLuint fbo;
    GLuint rbo;
    GLuint texture_id;
    int width;
    int height;
    void init(int width, int height);
    void destroy();
public:
    BufferManager(int width, int height);
    ~BufferManager();
    void bind();
    void unbind();
    void rescale(int width, int height);
    GLuint get_texture_id() const { return texture_id; }
};