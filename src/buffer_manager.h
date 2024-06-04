#pragma once

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGl/gl3ext.h>
#else
#include <GL/glew.h>
#endif

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