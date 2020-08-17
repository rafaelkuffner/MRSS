
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

#include <sstream>
#include <array>

#include "opengl.hpp"
#include "shader.hpp"
#include "util.hpp"

namespace {

	void set_uzfar(GLuint prog, float zfar) {
		if (cgu::glsl_uniform_zfar_name) {
			glUniform1f(glGetUniformLocation(prog, cgu::glsl_uniform_zfar_name), zfar);
		}
	}

}

namespace cgu {

	const char * glsl_uniform_zfar_name = "u_zfar";

	const char * glsl_frag_depth_source = CGU_GLSL_LINE_START R"(
	float frag_depth() {
		return gl_FragCoord.z;
	}
	)";

	const char * glsl_fullscreen_source = CGU_GLSL_LINE_START R"(
	#define CGU_FULLSCREEN
	uniform int u_layer_offset = 0;
	uniform float u_fragz = 0;
	#ifdef _VERTEX_
	flat out int v_instance_id;
	void main() {
		v_instance_id = gl_InstanceID;
	}
	#endif
	#ifdef _GEOMETRY_
	layout(points) in;
	layout(triangle_strip, max_vertices = 3) out;
	flat in int v_instance_id[];
	out vec2 v_tex_coord;
	void main() {
		int l = v_instance_id[0] + u_layer_offset;
		float z = u_fragz * 2 - 1;
		// output a single triangle that covers the whole screen
		gl_Position = vec4(-1, -1, z, 1);
		gl_Layer = l;
		v_tex_coord = vec2(0, 0);
		EmitVertex();
		gl_Position = vec4(3, -1, z, 1);
		gl_Layer = l;
		v_tex_coord = vec2(2, 0);
		EmitVertex();
		gl_Position = vec4(-1, 3, z, 1);
		gl_Layer = l;
		v_tex_coord = vec2(0, 2);
		EmitVertex();
		EndPrimitive();
	}
	#endif
	#ifdef _FRAGMENT_
	in vec2 v_tex_coord;
	layout(location=0) out vec4 f_color;
	#endif
	)";

	void draw_dummy(int instances, int vertices, GLenum mode) {
		assert(instances == decltype(instances)(GLsizei(instances)));
		assert(vertices == decltype(vertices)(GLsizei(vertices)));
		static GLuint vao = 0;
		if (!vao) glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
		glDrawArraysInstanced(mode, 0, vertices, instances);
		glBindVertexArray(0);
	}

	GLuint dummy_texture2d() {
		static GLuint tex = 0;
		if (!tex) {
			glGenTextures(1, &tex);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, tex);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			float pixel[]{0, 0, 0, 0};
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_FLOAT, pixel);
		}
		return tex;
	}

	GLuint use_dummy_shaderprog(glm::mat4 view, glm::mat4 proj, float zfar, const glm::vec4 &color) {
		const char *shader_source = CGU_GLSL_LINE_START R"(
		uniform mat4 u_projection;
		uniform mat4 u_modelview;
		uniform vec4 u_color;
		#ifdef _VERTEX_
		layout(location=0) in vec3 a_pos_m;
		void main() {
			vec4 pos_v = u_modelview * vec4(a_pos_m, 1.0);
			gl_Position = u_projection * pos_v;
		}
		#endif
		#ifdef _FRAGMENT_
		layout(location=0) out vec4 f_color;
		float frag_depth();
		void main() {
			f_color = u_color;
			gl_FragDepth = frag_depth();
		}
		#endif
		)";
		static GLuint prog = 0;
		if (!prog) {
			std::ostringstream hdr;
			hdr << "#ifdef _FRAGMENT_\n#define CGU_NEED_FRAG_DEPTH\n" << glsl_frag_depth_source << "\n#endif\n";
			prog = cgu::make_shader_program(
				"cgu::dummy", "330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{shader_source, hdr.str()}
			).release();
		}
		glUseProgram(prog);
		set_uzfar(prog, zfar);
		glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, value_ptr(color));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(proj));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(view));
		return prog;
	}

	void draw_texture2d(GLuint tex_color, GLuint tex_depth, float fragz) {
		if (!tex_color && !tex_depth) return;
		const char *shader_source = CGU_GLSL_LINE_START R"(
		#ifdef COLOR
		uniform sampler2D u_sampler_color;
		#endif
		#ifdef DEPTH
		uniform sampler2D u_sampler_depth;
		#endif
		#ifdef _FRAGMENT_
		void main() {
		#ifdef COLOR
			f_color = texture(u_sampler_color, v_tex_coord);
		#endif
		#ifdef DEPTH
			gl_FragDepth = texture(u_sampler_depth, v_tex_coord).r;
		#endif
		}
		#endif
		)";
		static GLuint prog_c = 0, prog_d = 0, prog_cd = 0;
		if (!prog_c) {
			prog_c = cgu::make_shader_program(
				"cgu::texture2d@c", "330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{"#define COLOR", glsl_fullscreen_source, shader_source}
			).release();
			prog_d = cgu::make_shader_program(
				"cgu::texture2d@d", "330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{"#define DEPTH", glsl_fullscreen_source, shader_source}
			).release();
			prog_cd = cgu::make_shader_program(
				"cgu::texture2d@cd", "330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{"#define COLOR", "#define DEPTH", glsl_fullscreen_source, shader_source}
			).release();
		}
		GLint at = 0;
		glGetIntegerv(GL_ACTIVE_TEXTURE, &at);
		// TODO only this target?
		glBindTexture(GL_TEXTURE_2D, tex_color);
		glActiveTexture(at + 1);
		glBindTexture(GL_TEXTURE_2D, tex_depth);
		GLuint prog = std::array<GLuint, 4>{0, prog_c, prog_d, prog_cd}[!!tex_color + 2 * !!tex_depth];
		glUseProgram(prog);
		glUniform1i(glGetUniformLocation(prog, "u_sampler_color"), at - GL_TEXTURE0);
		glUniform1i(glGetUniformLocation(prog, "u_sampler_depth"), at - GL_TEXTURE0 + 1);
		glUniform1f(glGetUniformLocation(prog, "u_fragz"), fragz);
		draw_dummy();
		glBindTexture(GL_TEXTURE_2D, 0);
		glActiveTexture(at);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void draw_axes(glm::mat4 view, glm::mat4 proj, float zfar, const draw_axes_params &params) {
		const char *shader_source = CGU_GLSL_LINE_START R"(
		uniform mat4 u_projection;
		uniform mat4 u_modelview;
		uniform float u_axislength;
		uniform vec4[6] u_colors;
		#ifdef _VERTEX_
		out int v_id;
		void main() {
			v_id = gl_VertexID;
		}
		#endif
		#ifdef _GEOMETRY_
		layout(points) in;
		layout(line_strip, max_vertices = 2) out;
		in int v_id[];
		flat out vec4 v_color;
		const vec3 dir[] = vec3[](
			vec3(1, 0, 0),
			vec3(0, 1, 0),
			vec3(0, 0, 1),
			vec3(-1, 0, 0),
			vec3(0, -1, 0),
			vec3(0, 0, -1)
		);
		void main() {
			vec4 pos0_v = u_modelview * vec4(0.0, 0.0, 0.0, 1.0);
			v_color = u_colors[v_id[0]];
			gl_Position = u_projection * pos0_v;
			EmitVertex();
			vec4 pos1_v = u_modelview * vec4(dir[v_id[0]] * u_axislength, 1.0);
			v_color = u_colors[v_id[0]];
			gl_Position = u_projection * pos1_v;
			EmitVertex();
			EndPrimitive();
		}
		#endif
		#ifdef _FRAGMENT_
		flat in vec4 v_color;
		layout(location=0) out vec4 f_color;
		float frag_depth();
		void main() {
			f_color = v_color;
			gl_FragDepth = frag_depth() - 0.000001;
		}
		#endif
		)";
		static GLuint prog = 0;
		if (!prog) {
			std::ostringstream hdr;
			hdr << "#ifdef _FRAGMENT_\n#define CGU_NEED_FRAG_DEPTH\n" << glsl_frag_depth_source << "\n#endif\n";
			prog = cgu::make_shader_program(
				"cgu::axes", "330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{shader_source, hdr.str()}
			).release();
		}
		glUseProgram(prog);
		set_uzfar(prog, zfar);
		glUniform1f(glGetUniformLocation(prog, "u_axislength"), params.axislength);
		glUniform4fv(glGetUniformLocation(prog, "u_colors"), 6, value_ptr(params.color_pos[0]));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(proj));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(view));
		cgu::draw_dummy(1, 6);
	}

	void draw_grid(glm::mat4 view, glm::mat4 proj, float zfar, const draw_grid_params &params) {
		const char *shader_source = CGU_GLSL_LINE_START R"(
		uniform mat4 u_projection;
		uniform mat4 u_modelview;
		uniform int u_divisions;
		uniform int u_id_offset;
		uniform float u_minorinterval;
		uniform float u_linelength;
		uniform float u_minor_depth_max;
		uniform float u_major_depth_max;
		uniform vec4 u_color_major;
		uniform vec4 u_color_minor;
		#ifdef _VERTEX_
		out int v_id;
		void main() {
			v_id = gl_VertexID - u_id_offset;
		}
		#endif
		#ifdef _GEOMETRY_
		layout(points) in;
		layout(line_strip, max_vertices = 2) out;
		in int v_id[];
		out VertexData {
			flat vec4 color;
			flat int id;
		} v_out;
		void main() {
			// TODO why is the abs() necessary?
			vec4 color =  abs(v_id[0]) % u_divisions == 0 ? u_color_major : u_color_minor;
			vec4 pos0_v = u_modelview * vec4(v_id[0] * u_minorinterval, 0, -u_linelength, 1);
			v_out.color = color;
			v_out.id = v_id[0];
			gl_Position = u_projection * pos0_v;
			EmitVertex();
			vec4 pos1_v = u_modelview * vec4(v_id[0] * u_minorinterval, 0, u_linelength, 1);
			v_out.color = color;
			v_out.id = v_id[0];
			gl_Position = u_projection * pos1_v;
			EmitVertex();
			EndPrimitive();
		}
		#endif
		#ifdef _FRAGMENT_
		in VertexData {
			flat vec4 color;
			flat int id;
		} v_in;
		layout(location=0) out vec4 f_color;
		float frag_depth();
		void main() {
			float depth_v = 1.0 / gl_FragCoord.w;
			bool ismajor = abs(v_in.id) % u_divisions == 0;
			float depthmax = ismajor ? u_major_depth_max : u_minor_depth_max;
			f_color = v_in.color;
			f_color.a *= 1.f - smoothstep(depthmax * 0.5, depthmax, depth_v);
			if (f_color.a < 0.01f) discard;
			gl_FragDepth = frag_depth() - (ismajor && v_in.id != 0 ? 0.000001 : 0.0);
		}
		#endif
		)";
		static GLuint prog = 0;
		if (!prog) {
			std::ostringstream hdr;
			hdr << "#ifdef _FRAGMENT_\n#define CGU_NEED_FRAG_DEPTH\n" << glsl_frag_depth_source << "\n#endif\n";
			prog = cgu::make_shader_program(
				"cgu::grid", "330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{shader_source, hdr.str()}
			).release();
		}
		const glm::mat4 rot = glm::rotate(glm::mat4(1), -glm::pi<float>() / 2.f, glm::vec3(0, 1, 0));
		glUseProgram(prog);
		set_uzfar(prog, zfar);
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(proj));
		glUniform1f(glGetUniformLocation(prog, "u_minor_depth_max"), params.minor_depth_max);
		glUniform1f(glGetUniformLocation(prog, "u_major_depth_max"), params.major_depth_max);
		glUniform4fv(glGetUniformLocation(prog, "u_color_major"), 1, value_ptr(params.color_major));
		glUniform4fv(glGetUniformLocation(prog, "u_color_minor"), 1, value_ptr(params.color_minor));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(view));
		const int xcount = int(params.divisions.x * params.extent.x / params.interval.x) * 2 + 1;
		glUniform1i(glGetUniformLocation(prog, "u_divisions"), params.divisions.x);
		glUniform1i(glGetUniformLocation(prog, "u_id_offset"), xcount / 2);
		glUniform1f(glGetUniformLocation(prog, "u_minorinterval"), params.interval.x / params.divisions.x);
		glUniform1f(glGetUniformLocation(prog, "u_linelength"), params.extent.y);
		cgu::draw_dummy(1, xcount);
		const int ycount = int(params.divisions.y * params.extent.y / params.interval.y) * 2 + 1;
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(view * rot));
		glUniform1i(glGetUniformLocation(prog, "u_divisions"), params.divisions.y);
		glUniform1i(glGetUniformLocation(prog, "u_id_offset"), ycount / 2);
		glUniform1f(glGetUniformLocation(prog, "u_minorinterval"), params.interval.y / params.divisions.y);
		glUniform1f(glGetUniformLocation(prog, "u_linelength"), params.extent.x);
		cgu::draw_dummy(1, ycount);
	}

	bool gl_rendertarget::bind(glm::ivec3 size_) {
		if (!tex) {
			tex = gl_object::gen_texture();
			glBindTexture(params.target, tex);
			glTexParameteri(params.target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(params.target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(params.target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(params.target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(params.target, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		} else {
			glBindTexture(params.target, tex);
		}
		if (size == size_) return false;
		size = size_;
		if (size.z == 0) {
			glTexImage2D(params.target, 0, params.internalformat, size.x, size.y, 0, params.dataformat, params.datatype, nullptr);
		} else {
			glTexImage3D(params.target, 0, params.internalformat, size.x, size.y, size.z, 0, params.dataformat, params.datatype, nullptr);
		}
		return true;
	}

	bool gl_framebuffer::bind(GLenum target, glm::ivec3 size_) {
		auto bind = [&]() {
			glBindFramebuffer(target, fbo);
			if (target != GL_READ_FRAMEBUFFER) glDrawBuffers(drawbuffers.size(), drawbuffers.data());
		};
		if (!fbo) {
			fbo = cgu::gl_object::gen_framebuffer();
			bind();
			for (auto &rt : rendertargets) {
				// need to gen the texture before we can attach
				rt.bind(size_);
				glFramebufferTexture(target, rt.params.attachment, rt.tex, 0);
			}
			return true;
		} else {
			bool r = false;
			bind();
			for (auto &rt : rendertargets) {
				// resize
				r |= rt.bind(size_);
			}
			return r;
		}
	}

	const char * gl_debug_source_string(GLenum source) {
		switch (source) {
		case GL_DEBUG_SOURCE_API:
			return "API";
		case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
			return "Window System";
		case GL_DEBUG_SOURCE_SHADER_COMPILER:
			return "Shader Compiler";
		case GL_DEBUG_SOURCE_THIRD_PARTY:
			return "Third Party";
		case GL_DEBUG_SOURCE_APPLICATION:
			return "Application";
		case GL_DEBUG_SOURCE_OTHER:
			return "Other";
		default:
			return "n/a";
		}
	}

	const char * gl_debug_severity_string(GLenum severity) {
		switch (severity) {
		case GL_DEBUG_SEVERITY_HIGH:
			return "High";
		case GL_DEBUG_SEVERITY_MEDIUM:
			return "Medium";
		case GL_DEBUG_SEVERITY_LOW:
			return "Low";
		case GL_DEBUG_SEVERITY_NOTIFICATION:
			return "None";
		default:
			return "n/a";
		}
	}

	const char * gl_debug_type_string(GLenum type) {
		switch (type) {
		case GL_DEBUG_TYPE_ERROR:
			return "Error";
		case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
			return "Deprecated Behaviour";
		case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
			return "Undefined Behaviour";
		case GL_DEBUG_TYPE_PORTABILITY:
			return "Portability";
		case GL_DEBUG_TYPE_PERFORMANCE:
			return "Performance";
		case GL_DEBUG_TYPE_OTHER:
			return "Other";
		default:
			return "n/a";
		}
	}

}
