
/*
BSD 3-Clause License

Copyright (c) 2013-2019, Benjamin Allen and Joshua Scott
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <vector>
#include <algorithm>

// include glew.h before (instead of) gl.h, or anything that includes gl.h
// glew.h replaces gl.h and sets up OpenGL functions in a cross-platform manner
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "util.hpp"

// prepend glsl source string literals with this to get line numbers
// e.g. const char *source = CGU_GLSL_LINE_START R"( ...
// there must be no line break before start of string literal
#define CGU_GLSL_LINE_START "#line " CGU_STRINGIFY(__LINE__) "\n"

namespace cgu {

	// name of the uniform used for view-space far plane depth, or null
	// defaults to "u_zfar"
	extern const char * glsl_uniform_zfar_name;

	// set this before using any cgu draw functions if required
	// should provide definition of "float frag_depth()"
	// by default, returns value of gl_FragCoord.z
	// viewspace depth can be determined by (1.0 / gl_FragCoord.w)
	// the uniform identified by glsl_uniform_zfar_name will be set
	// CGU_NEED_FRAG_DEPTH is defined
	extern const char * glsl_frag_depth_source;

	// provides the setup for a 'fullscreen' shader
	// use draw_dummy with 1 vertex
	// uses vertex, geometry and fragment shaders
	// should do everything apart from define fragment main()
	// should turn 1 vertex into a fullscreen triangle
	// should define CGU_FULLSCREEN
	// should provide these uniforms:
	// - uniform int u_layer_offset = 0;
	// - uniform float u_fragz = 0;
	// should provide these fragment interface variables:
	// - in vec2 v_tex_coord;
	// - layout(location=0) out vec4 f_color;
	// gl_InstanceID + u_layer_offset should be written to gl_Layer
	// generally no need to replace this
	extern const char * glsl_fullscreen_source;

	// draws an empty VAO
	// can be used for shaders that do all the work
	// not multi-context-thread safe
	void draw_dummy(int instances = 1, int vertices = 1, GLenum mode = GL_POINTS);

	// 1x1 black 2d texture
	// not multi-context-thread safe
	GLuint dummy_texture2d();

	// binds and returns a dummy shader program that draws a mesh in the specified color
	// expects vertex positions as attribute 0
	// not multi-context-thread safe
	GLuint use_dummy_shaderprog(glm::mat4 view, glm::mat4 proj, float zfar, const glm::vec4 &color);

	// blit 2d color and depth textures to the current viewport
	// will bind textures to active texture unit + {0, 1}
	// either texture can be 0 to omit
	// uses specified fragment z [0,1] if no depth texture
	// note that if the depth test is disabled, no depth can be written
	// not multi-context-thread safe
	void draw_texture2d(GLuint tex_color, GLuint tex_depth, float fragz);

	struct draw_axes_params {
		float axislength = 1000;
		glm::vec4 color_pos[3]{{1.0f, 0, 0, 0.8f}, {0, 1.0f, 0, 0.8f}, {0, 0, 1.0f, 0.8f}};
		glm::vec4 color_neg[3]{{0.5f, 0, 0, 0.8f}, {0, 0.5f, 0, 0.8f}, {0, 0, 0.5f, 0.8f}};
	};

	struct draw_grid_params {
		glm::vec2 extent{10};
		glm::vec2 interval{1};
		glm::ivec2 divisions{10};
		glm::vec4 color_major{0.5f, 0.5f, 0.5f, 0.6f};
		glm::vec4 color_minor{0.5f, 0.5f, 0.5f, 0.3f};
		float minor_depth_max = 20;
		float major_depth_max = 100;
	};

	// draw lines for x/y/z axes
	// not multi-context-thread safe
	void draw_axes(glm::mat4 view, glm::mat4 proj, float zfar, const draw_axes_params & = draw_axes_params{});

	// draw grid on x/z plane
	// not multi-context-thread safe
	void draw_grid(glm::mat4 view, glm::mat4 proj, float zfar, const draw_grid_params & = draw_grid_params{});

	// move-only RAII wrapper for GL objects.
	// stores the GL name (id) of the object.
	// deletes the object (name) on destruction with the stored deleter.
	class [[nodiscard]] gl_object {
	public:
		using deleter_t = void (APIENTRY *)(GLsizei, const GLuint *);

	private:
		GLuint m_id = 0;
		deleter_t m_dtor;

		void destroy() noexcept {
			if (m_id) {
				m_dtor(1, &m_id);
				m_id = 0;
			}
		}

	public:
		// empty object
		gl_object() { }

		// takes an existing GL name (id) and a pointer to a deleter function
		gl_object(GLuint id_, deleter_t dtor_) : m_id(id_), m_dtor(dtor_) { }

		// remove copy ctors
		gl_object(const gl_object &) = delete;
		gl_object & operator=(const gl_object &) = delete;

		gl_object(gl_object &&other) noexcept {
			m_id = other.m_id;
			m_dtor = other.m_dtor;
			other.m_id = 0;
		}

		gl_object & operator=(gl_object &&other) noexcept {
			destroy();
			m_id = other.m_id;
			m_dtor = other.m_dtor;
			other.m_id = 0;
			return *this;
		}

		// returns the GL name of the object
		GLuint get() const noexcept {
			return m_id;
		}

		// implicit GLuint converter
		// returns the GL name of the object
		operator GLuint() const noexcept {
			return m_id;
		}

		// explicit boolean converter
		// true IFF GL name is not zero
		explicit operator bool() const noexcept {
			return m_id;
		}

		// true IFF GL name is zero
		bool operator!() const noexcept {
			return !m_id;
		}

		// relinquishes ownership of the GL object
		// returns the GL name of the object and zeros the stored name
		GLuint release() noexcept {
			GLuint id = m_id;
			m_id = 0;
			return id;
		}

		~gl_object() {
			destroy();
		}

		static gl_object gen_buffer() {
			GLuint o;
			glGenBuffers(1, &o);
			return {o, glDeleteBuffers};
		}

		static gl_object gen_vertex_array() {
			GLuint o;
			glGenVertexArrays(1, &o);
			return {o, glDeleteVertexArrays};
		}

		static gl_object gen_texture() {
			GLuint o;
			glGenTextures(1, &o);
			return {o, glDeleteTextures};
		}

		static gl_object gen_framebuffer() {
			GLuint o;
			glGenFramebuffers(1, &o);
			return {o, glDeleteFramebuffers};
		}

		static gl_object gen_shader(GLenum type) {
			GLuint o = glCreateShader(type);
			return {o, [](GLsizei, const GLuint *o) { glDeleteShader(*o); }};
		}

		static gl_object gen_program() {
			GLuint o = glCreateProgram();
			return {o, [](GLsizei, const GLuint *o) { glDeleteProgram(*o); }};
		}
	};
	
	struct gl_rendertarget_params {
		GLenum attachment = 0;
		GLenum target = GL_TEXTURE_2D;
		GLenum internalformat = GL_RGBA8;
		GLenum dataformat = GL_RGBA;
		GLenum datatype = GL_FLOAT;
	};

	struct gl_rendertarget {
		// 3d single layer not supported
		gl_object tex;
		glm::ivec3 size{0};
		gl_rendertarget_params params;

		gl_rendertarget() {}

		explicit gl_rendertarget(const gl_rendertarget_params &params_)
			: params(params_)
		{}

		void bind(glm::ivec2 size_) {
			bind(glm::ivec3{size_, 0});
		}

		void bind(glm::ivec3 size_);
	};

	struct gl_framebuffer {
		gl_object fbo;
		
		// draw buffers will be respecified whenever bind() is used for a draw target
		std::vector<GLenum> drawbuffers;

		// to change attachments, destroy the fbo so that it gets recreated
		std::vector<gl_rendertarget> rendertargets;

		gl_framebuffer() {}

		// draw buffers set to all color attachments in init order
		gl_framebuffer(std::initializer_list<gl_rendertarget_params> rendertargets_) {
			for (auto &rtp : rendertargets_) {
				rendertargets.emplace_back(rtp);
				if (rtp.attachment >= GL_COLOR_ATTACHMENT0 && rtp.attachment <= GL_COLOR_ATTACHMENT15) {
					drawbuffers.push_back(rtp.attachment);
				}
			}
		}

		gl_framebuffer(std::initializer_list<gl_rendertarget_params> rendertargets_, std::initializer_list<GLenum> drawbuffers_)
			: drawbuffers(drawbuffers_)
		{
			for (auto &rtp : rendertargets_) {
				rendertargets.emplace_back(rtp);
			}
		}

		const gl_rendertarget & operator[](GLenum attachment) const {
			auto it = std::find_if(
				rendertargets.begin(), rendertargets.end(), 
				[&](auto &rt) { return rt.params.attachment == attachment; }
			);
			assert(it != rendertargets.end());
			return *it;
		}

		gl_rendertarget & operator[](GLenum attachment) {
			auto it = std::find_if(
				rendertargets.begin(), rendertargets.end(), 
				[&](auto &rt) { return rt.params.attachment == attachment; }
			);
			assert(it != rendertargets.end());
			return *it;
		}

		void bind(GLenum target, glm::ivec2 size_) {
			bind(target, glm::ivec3{size_, 0});
		}

		void bind(GLenum target, glm::ivec3 size_);

	};

	const char * gl_debug_source_string(GLenum source);

	const char * gl_debug_severity_string(GLenum severity);

	const char * gl_debug_type_string(GLenum type);

}
