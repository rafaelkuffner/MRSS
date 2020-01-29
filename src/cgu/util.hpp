
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

#pragma once

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define CGU_STRINGIFY1(a) #a
#define CGU_STRINGIFY(a) CGU_STRINGIFY1(a)

namespace cgu {

	template <typename ...ArgTs>
	inline std::string stringify(const ArgTs &...args) {
		std::ostringstream oss;
		(oss << ... << args);
		return oss.str();
	}

	struct orbital_camera {
		glm::vec3 focus{0};
		float cam_pitch = 0;
		float cam_yaw = 0;
		float cam_distance = 5;
		float zoom_factor = 1.1f;
		glm::vec2 last_mouse_pos{0};

		void update(bool mouse_down, glm::vec2 mouse_pos, glm::ivec2 winsize) {
			if (mouse_down) {
				glm::vec2 winhalfsize = glm::vec2(winsize) * 0.5f;
				// clamp the pitch to [-pi/2, pi/2]
				cam_pitch += glm::acos(glm::clamp((last_mouse_pos.y - winhalfsize.y) / winhalfsize.y, -1.f, 1.f))
					- glm::acos(glm::clamp((mouse_pos.y - winhalfsize.y) / winhalfsize.y, -1.f, 1.f));
				cam_pitch = glm::clamp(cam_pitch, -glm::pi<float>() * 0.5f, glm::pi<float>() * 0.5f);
				// wrap the yaw to [-pi, pi]
				cam_yaw += glm::acos(glm::clamp((last_mouse_pos.x - winhalfsize.x) / winhalfsize.x, -1.f, 1.f))
					- glm::acos(glm::clamp((mouse_pos.x - winhalfsize.x) / winhalfsize.x, -1.f, 1.f));
				if (cam_yaw > glm::pi<float>()) {
					cam_yaw -= 2 * glm::pi<float>();
				} else if (cam_yaw < -glm::pi<float>()) {
					cam_yaw += 2 * glm::pi<float>();
				}
			}
			last_mouse_pos = mouse_pos;
		}

		void zoom(float scroll) {
			cam_distance *= pow(zoom_factor, -scroll);
		}

		glm::mat4 view() {
			glm::mat4 m{1};
			m = glm::translate(m, glm::vec3(0, 0, -cam_distance));
			m = glm::rotate(m, cam_pitch, glm::vec3(1, 0, 0));
			m = glm::rotate(m, cam_yaw, glm::vec3(0, 1, 0));
			m = glm::translate(m, -focus);
			return m;
		}

	};

}

namespace glm {

	template <typename T, int N, glm::precision P>
	inline std::ostream & operator<<(std::ostream &out, const glm::vec<N, T, P> &v) {
		out << '(';
		if (N > 0) out << v[0];
		for (size_t i = 1; i < N; ++i) {
			out << ", " << v[i];
		}
		out << ')';
		return out;
	}

	template <typename T, size_t Cols, size_t Rows, glm::precision P>
	inline std::ostream & operator<<(std::ostream &out, const glm::mat<Cols, Rows, T, P> &m) {
		static const char *toplinestart = Rows > 1 ? "/ " : "[ ";
		static const char *toplineend = Rows > 1 ? " \\" : " ]";
		static const char *midlinestart = "| ";
		static const char *midlineend = " |";
		static const char *bottomlinestart = Rows > 1 ? "\\ " : "[ ";
		static const char *bottomlineend = Rows > 1 ? " /" : " ]";
		// calculate column widths
		size_t maxwidths[Cols];
		std::ostringstream tempss;
		tempss.precision(out.precision());
		for (size_t i = 0; i < Cols; ++i) {
			maxwidths[i] = 0;
			for (size_t j = 0; j < Rows; ++j) {
				// msvc: seeking to 0 on empty stream sets failbit
				tempss.seekp(0);
				tempss.clear();
				tempss << m[i][j];
				maxwidths[i] = size_t(tempss.tellp()) > maxwidths[i] ? size_t(tempss.tellp()) : maxwidths[i];
			}
		}
		// print
		for (size_t i = 0; i < Rows; ++i) {
			if (i == 0) {
				out << toplinestart;
			} else if (i + 1 == Rows) {
				out << bottomlinestart;
			} else {
				out << midlinestart;
			}
			for (size_t j = 0; j < Cols; ++j) {
				out.width(maxwidths[j]);
				out << m[j][i];
				if (j + 1 < Cols) out << ", ";
			}
			if (i == 0) {
				out << toplineend;
			} else if (i + 1 == Rows) {
				out << bottomlineend;
			} else {
				out << midlineend;
			}
			if (i + 1 < Rows) out << '\n';
		}
		return out;
	}

}
