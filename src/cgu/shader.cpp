
/*
BSD 3-Clause License

Copyright (c) 2013-2019, Benjamin Allen
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "shader.hpp"

namespace {

	void print_shader_info_log(GLuint obj) {
		int infologLength = 0;
		int charsWritten = 0;
		glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infologLength);
		if (infologLength > 1) {
			std::vector<char> infoLog(infologLength);
			glGetShaderInfoLog(obj, infologLength, &charsWritten, &infoLog[0]);
			std::cerr << &infoLog[0] << std::endl;
		}
	}

	void print_program_info_log(GLuint obj) {
		int infologLength = 0;
		int charsWritten = 0;
		glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);
		if (infologLength > 1) {
			std::vector<char> infoLog(infologLength);
			glGetProgramInfoLog(obj, infologLength, &charsWritten, &infoLog[0]);
			std::cerr << &infoLog[0] << std::endl;
		}
	}

	cgu::gl_object compile_shader(GLenum type, std::string_view text) {
		cgu::gl_object shader = cgu::gl_object::gen_shader(type);
		const char *text_c = text.data();
		GLint textlen = GLint(text.size());
		assert(text.size() == decltype(text.size())(textlen));
		glShaderSource(shader, 1, &text_c, &textlen);
		glCompileShader(shader);
		GLint compile_status;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);
		// always print, so we can see warnings
		print_shader_info_log(shader);
		if (!compile_status) throw cgu::shader_error("shader compile failed");
		return shader;
	}

	void link_shader_program(GLuint prog) {
		glLinkProgram(prog);
		GLint link_status;
		glGetProgramiv(prog, GL_LINK_STATUS, &link_status);
		// always print, so we can see warnings
		print_program_info_log(prog);
		if (!link_status) throw cgu::shader_error("shader link failed");
	}

}

namespace cgu {

	gl_object make_shader_program(
		std::string_view description,
		std::string_view profile,
		const std::vector<GLenum> &stypes,
		const std::vector<std::string_view> &sources
	) {
		gl_object prog = gl_object::gen_program();
		std::cerr << "building shader program (" << prog.get() << ") " << description << std::endl;

		auto get_define = [](GLenum stype) {
			switch (stype) {
			case GL_VERTEX_SHADER:
				return "_VERTEX_";
			case GL_GEOMETRY_SHADER:
				return "_GEOMETRY_";
			case GL_TESS_CONTROL_SHADER:
				return "_TESS_CONTROL_";
			case GL_TESS_EVALUATION_SHADER:
				return "_TESS_EVALUATION_";
			case GL_FRAGMENT_SHADER:
				return "_FRAGMENT_";
			default:
				return "_UNKNOWN_";
			}
		};

		for (auto stype : stypes) {
			std::ostringstream oss;
			oss << "#version " << profile << std::endl;
			oss << "#define " << get_define(stype) << std::endl;
			for (auto &source : sources) {
				oss << "\n#line 1 " << (&source - sources.data()) << std::endl;
				oss << source;
			}
			std::cerr << "[shader " << get_define(stype) << " " << description << "]" << std::endl;
			auto shader = compile_shader(stype, oss.str());
			glAttachShader(prog, shader);
		}

		std::cerr << "[program " << description << "]" << std::endl;
		link_shader_program(prog);
		std::cerr << "shader program compiled and linked successfully" << std::endl;
		return prog;
	}

	gl_object make_shader_program_from_files(
		std::string_view profile,
		const std::vector<GLenum> &stypes,
		const std::vector<std::string_view> &headers,
		const std::vector<std::filesystem::path> &srcpaths
	) {
		if (srcpaths.empty()) throw shader_error("no source files");
		std::string desc;
		std::vector<std::string_view> source_views = headers;
		for (size_t i = 0; i < headers.size(); i++) {
			desc += 'h';
			desc += std::to_string(i);
			desc += ';';
		}
		std::string source;
		for (const auto &srcpath : srcpaths) {
			std::cerr << "sourcing shader from " << srcpath << std::endl;
			desc += srcpath.filename().string();
			desc += ';';
			std::ifstream ifs{srcpath};
			if (!ifs.is_open()) throw shader_error("failed to open file");
			source += "#line 1 ";
			source += std::to_string(&srcpath - srcpaths.data() + headers.size());
			source += '\n';
			while (ifs.good()) {
				std::string line;
				std::getline(ifs, line);
				source += line;
				source += '\n';
			}
		}
		source_views.push_back(source);
		desc.pop_back();
		return make_shader_program(desc, profile, stypes, source_views);
	}

}
