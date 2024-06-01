/**
 * @file staplegl.hpp

 * @brief StapleGL, a C++20 wrapper for OpenGL.
 *
 * @authors Christian Panov, Dario Loi.
 * @date 2023-04-28
 *
 * @details StapleGL is a C++20 wrapper for OpenGL. It is designed to be simple and easy to use.
 * simply import this header file and provide your own OpenGL loader in modules/gl_functions.h.
 *
 *
 * @see modules/gl_functions.h
 * @copyright MIT License
 *
 */

#pragma once

#include "modules/cubemap.hpp"
#include "modules/framebuffer.hpp"
#include "modules/index_buffer.hpp"
#include "modules/shader.hpp"
#include "modules/texture.hpp"
#include "modules/uniform_buffer.hpp"
#include "modules/vertex_array.hpp"
#include "modules/vertex_buffer.hpp"
#include "modules/vertex_buffer_inst.hpp"
#include "modules/vertex_buffer_layout.hpp"
#include "modules/renderbuffer.hpp"

/*

The rest of the file is reserved for documentation purposes. (mainly for doxygen's main page)

*/

/**
 * @mainpage StapleGL
 *
 * @section intro_sec Introduction
 *
 * StapleGL is a C++20 wrapper for OpenGL. It is designed with the principal goal of having
 * zero dependencies and being as lightweight as possible. You can integrate it
 * with any OpenGL project by providing your own OpenGL loader in modules/gl_functions.h,
 * and then simply importing this header file.
 *
 * StapleGL offers modern RAII interfaces for OpenGL objects, such as Vertex Buffer Objects,
 * abstracting away the need to manually manage OpenGL objects, and lessening the amount
 * of boilerplate code you need to write in order to get started with OpenGL.
 *
 * @section install_sec Installation
 *
 * To install StapleGL, simply clone the repository and copy the contents of the include folder
 * into your project's include folder. Then, provide your own OpenGL loader in modules/gl_functions.h.
 *
 * @section usage_sec Usage
 *
 * Each of the modules in StapleGL has its own documentation page which detail its usage, you can find the documentation
 * pages in the /docs/html folder, by opening index.html.
 *
 * In general, due to StapleGL being very minimalistic, its very easy to get started with it,
 * and you can find a number of examples in the examples folder.
 * 
 * The examples can be compiled through the CMakeLists.txt file in the root directory.
 *
 */