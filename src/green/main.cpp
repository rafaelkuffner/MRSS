
#include <iostream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <thread>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <glm/gtc/matrix_transform.hpp>

#include <cgu/opengl.hpp>
#include <cgu/util.hpp>
#include <cgu/mesh.hpp>
#include <cgu/shader.hpp>

#include "model.hpp"

using namespace std;
using namespace green;

namespace {
	void drop_callback(GLFWwindow *win, int count, const char **paths);
	void cursor_pos_callback(GLFWwindow *win, double xpos, double ypos);
	void mouse_button_callback(GLFWwindow *win, int button, int action, int mods);
	void scroll_callback(GLFWwindow *win, double xoffset, double yoffset);
	void key_callback(GLFWwindow *win, int key, int scancode, int action, int mods);
	void char_callback(GLFWwindow *win, unsigned int c);
	void APIENTRY gl_debug_callback(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar *, const GLvoid *);

	const char *glsl_depth_env = R"(
	#ifndef DEPTH_ENV
	#define DEPTH_ENV
	uniform float u_zfar;
	float encode_depth(float depth_v) {
		const float c = 0.01;
		float fc = 1.0 / log(u_zfar * c + 1.0);
		return log(depth_v * c + 1.0) * fc;
	}
	float decode_depth(float depth_p) {
		const float c = 0.01;
		float fc = 1.0 / log(u_zfar * c + 1.0);
		return (exp(depth_p / fc) - 1.0) / c;
	}
	#ifdef _FRAGMENT_
	float frag_depth() {
		return encode_depth(1.0 / gl_FragCoord.w);
	}
	#endif
	#endif
	)";

	ImGuiContext *imguictx;

	float zfar = 1000000;
	cgu::orbital_camera cam;
	
	cgu::gl_framebuffer fb_scene{
		cgu::gl_rendertarget_params{GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT, GL_FLOAT},
		cgu::gl_rendertarget_params{GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, GL_RGBA16F, GL_RGBA, GL_FLOAT},
		cgu::gl_rendertarget_params{GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, GL_RG32I, GL_RG_INTEGER, GL_INT}
	};

	std::vector<std::unique_ptr<Entity>> entities;

	GLsync sync_read_ids = nullptr;
	cgu::gl_object buf_read_ids, buf_read_depths;
	glm::ivec2 pos_read_ids{0};
	const int size_read_ids = 64;
	selection cur_sel;
	float cur_depth = zfar;

	float decode_depth(float depth_p) {
		const float c = 0.01;
		float fc = 1.0 / log(zfar * c + 1.0);
		return (exp(depth_p / fc) - 1.0) / c;
	}

	std::pair<glm::vec3, glm::vec3> unproject(const glm::mat4 &proj_inv, const glm::vec2 &pos_p, float depth_v) {
		// near plane position
		glm::vec3 pos0_p = glm::vec3(pos_p, -1.0);
		glm::vec4 pos0_vh = proj_inv * glm::vec4(pos0_p, 1.0);
		glm::vec3 pos0_v = glm::vec3(pos0_vh) / pos0_vh.w;
		// not-far plane position
		// TODO better not-far-plane?
		glm::vec3 pos1_p = glm::vec3(pos_p, 0.0);
		glm::vec4 pos1_vh = proj_inv * glm::vec4(pos1_p, 1.0);
		glm::vec3 pos1_v = glm::vec3(pos1_vh) / pos1_vh.w;
		glm::vec3 dir = normalize(pos1_v - pos0_v);
		glm::vec3 pos = pos0_v + dir * ((depth_v + pos0_v.z) / -dir.z);
		return {pos, dir};
	}

	void update_hover(const glm::ivec2 &mouse_pos_fb) {
		if (!buf_read_ids) buf_read_ids = cgu::gl_object::gen_buffer();
		if (!buf_read_depths) buf_read_depths = cgu::gl_object::gen_buffer();
		if (sync_read_ids) {
			if (glClientWaitSync(sync_read_ids, 0, 0) != GL_TIMEOUT_EXPIRED) {
				glDeleteSync(sync_read_ids);
				sync_read_ids = nullptr;
				glBindBuffer(GL_ARRAY_BUFFER, buf_read_ids);
				auto *iddata = reinterpret_cast<glm::ivec2 *>(
					glMapBufferRange(GL_ARRAY_BUFFER, 0, sizeof(glm::ivec2) * size_read_ids * size_read_ids, GL_MAP_READ_BIT)
				);
				float best_entity_dist = 9001;
				float best_vertex_dist = 9001;
				cur_sel.hover_entity = -1;
				cur_sel.hover_vertex = -1;
				for (int y = 0; y < size_read_ids; y++) {
					for (int x = 0; x < size_read_ids; x++) {
						glm::ivec2 ids = iddata[size_read_ids * y + x];
						if (ids.x != -1) {
							float dist = glm::length(glm::vec2{pos_read_ids} + glm::vec2{x, y} - glm::vec2{mouse_pos_fb});
							if (dist < best_entity_dist && cur_sel.hover_vertex == -1) {
								best_entity_dist = dist;
								cur_sel.hover_entity = ids.x;
							}
							if (dist < best_vertex_dist && ids.y != -1) {
								best_entity_dist = dist;
								best_vertex_dist = dist;
								cur_sel.hover_entity = ids.x;
								cur_sel.hover_vertex = ids.y;
							}
						}
					}
				}
				glUnmapBuffer(GL_ARRAY_BUFFER);
				glBindBuffer(GL_ARRAY_BUFFER, buf_read_depths);
				auto *depthdata = reinterpret_cast<float *>(
					glMapBufferRange(GL_ARRAY_BUFFER, 0, sizeof(float) * size_read_ids * size_read_ids, GL_MAP_READ_BIT)
				);
				auto cur_xy = clamp(mouse_pos_fb - pos_read_ids, glm::ivec2(0), glm::ivec2(size_read_ids - 1));
				float depth = depthdata[size_read_ids * cur_xy.y + cur_xy.x];
				cur_depth = decode_depth(depth);
				glUnmapBuffer(GL_ARRAY_BUFFER);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				if (ImGui::Begin("Hover")) {
					ImGui::Text("Entity %d; Vertex %d", cur_sel.hover_entity, cur_sel.hover_vertex);
					ImGui::Text("Depth %f", cur_depth);
				}
				ImGui::End();
			}
		}
		if (!sync_read_ids) {
			// need fbo already setup
			glBindFramebuffer(GL_READ_FRAMEBUFFER, fb_scene.fbo);
			pos_read_ids = mouse_pos_fb - size_read_ids / 2;
			glBindBuffer(GL_PIXEL_PACK_BUFFER, buf_read_ids);
			glBufferData(GL_PIXEL_PACK_BUFFER, sizeof(glm::ivec2) * size_read_ids * size_read_ids, nullptr, GL_STREAM_READ);
			glReadBuffer(GL_COLOR_ATTACHMENT1);
			glReadPixels(pos_read_ids.x, pos_read_ids.y, size_read_ids, size_read_ids, GL_RG_INTEGER, GL_INT, nullptr);
			glBindBuffer(GL_PIXEL_PACK_BUFFER, buf_read_depths);
			glBufferData(GL_PIXEL_PACK_BUFFER, sizeof(glm::ivec2) * size_read_ids * size_read_ids, nullptr, GL_STREAM_READ);
			// read buffer only specified color, using depth component will read depth
			glReadBuffer(GL_NONE);
			glReadPixels(pos_read_ids.x, pos_read_ids.y, size_read_ids, size_read_ids, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
			sync_read_ids = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
			glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
			glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
		}
	}

	void render(GLFWwindow *window) {

		glm::ivec2 fbsize;
		glfwGetFramebufferSize(window, &fbsize.x, &fbsize.y);
		glViewport(0, 0, fbsize.x, fbsize.y);

		// FIXME when fb size != win size (hidpi)
		glm::dvec2 mouse_pos{0};
		glfwGetCursorPos(window, &mouse_pos.x, &mouse_pos.y);
		glm::ivec2 mouse_pos_fb = {int(mouse_pos.x), fbsize.y - int(mouse_pos.y) - 1};

		fb_scene.bind(GL_DRAW_FRAMEBUFFER, fbsize);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glClearBufferfv(GL_DEPTH, 0, value_ptr(glm::vec4{1.f}));
		glClearBufferfv(GL_COLOR, 0, value_ptr(glm::vec4{0.1f, 0.1f, 0.2f, 1.0f}));
		glClearBufferiv(GL_COLOR, 1, value_ptr(glm::ivec4{-1}));
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glm::mat4 proj = glm::perspective(1.f, float(fbsize.x) / fbsize.y, 0.1f, zfar);
		glm::mat4 view = cam.view();

		auto proj_inv = inverse(proj);
		auto view_inv = inverse(view);

		auto [view_pos, view_dir] = unproject(proj_inv, glm::vec2(mouse_pos_fb) / glm::vec2(fbsize) * 2.f - 1.f, cur_depth);
		auto world_pos = glm::vec3(view_inv * glm::vec4(view_pos, 1));
		if (ImGui::Begin("Hover")) {
			ImGui::Text("Position x=%.2f y=%.2f z=%.2f", world_pos.x, world_pos.y, world_pos.z);
		}
		ImGui::End();

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glPointSize(3);

		bool can_hover = false;
		for (auto it = entities.begin(); it != entities.end(); ) {
			can_hover = true;
			auto &e = *it;
			e->draw(view, proj, zfar, cur_sel);
			if (e->dead()) {
				it = entities.erase(it);
			} else {
				++it;
			}
		}


		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glDepthFunc(GL_ALWAYS);
		glEnable(GL_FRAMEBUFFER_SRGB);
		cgu::draw_texture2d(fb_scene[GL_COLOR_ATTACHMENT0].tex, fb_scene[GL_DEPTH_ATTACHMENT].tex, 0);
		glDisable(GL_FRAMEBUFFER_SRGB);

		glDepthFunc(GL_LEQUAL);

		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		cgu::draw_grid(view, proj, zfar);
		cgu::draw_axes(view, proj, zfar);
		glDisable(GL_BLEND);

		glDisable(GL_DEPTH_TEST);

		// clearing the id buffer may not actually execute if we dont draw anything?! (driver bug?)
		if (can_hover) update_hover(mouse_pos_fb);

	}

}

int main() {

	if (!glfwInit()) {
		cerr << "Error: Could not initialize GLFW" << endl;
		abort();
	}

	// GL 3.3 core context
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// disallow legacy (for OS X)
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// request a debug context so we get debug callbacks
	// remove this for possible GL performance increases
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);

	// request a window that can perform gamma correction
	glfwWindowHint(GLFW_SRGB_CAPABLE, true);

	GLFWwindow *window = glfwCreateWindow(800, 600, "GREEN Mesh Saliency", nullptr, nullptr);
	if (!window) {
		cerr << "Error: Could not create GLFW window" << endl;
		abort();
	}

	glfwMakeContextCurrent(window);

	// required for full GLEW functionality for GL3+
	glewExperimental = GL_TRUE;

	// initialize GLEW
	// must be done after making a GL context current
	GLenum glew_err = glewInit();
	if (GLEW_OK != glew_err) {
		cerr << "Error: Could not initialize GLEW: " << glewGetErrorString(glew_err) << endl;
		abort();
	}

	// print out versions
	cout << "OpenGL : " << glGetString(GL_VERSION) << endl;
	cout << "GLEW   : " << glewGetString(GLEW_VERSION) << endl;
	cout << "GLFW   : " << glfwGetVersionString() << endl;
	cout << "ImGui  : " << ImGui::GetVersion() << endl;

	// enable GL_ARB_debug_output if available
	if (glfwExtensionSupported("GL_ARB_debug_output")) {
		// this allows the error location to be determined from a stacktrace
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
		// setup up the callback
		glDebugMessageCallbackARB(gl_debug_callback, nullptr);
		glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, true);
		cout << "GL_ARB_debug_output callback installed" << endl;
	} else {
		cout << "GL_ARB_debug_output not available" << endl;
	}

	// initialize ImGui
	imguictx = ImGui::CreateContext();
	ImGui::SetCurrentContext(imguictx);
	ImGui_ImplGlfw_InitForOpenGL(window, false);
	ImGui_ImplOpenGL3_Init();

	// attach input callbacks to window
	glfwSetCursorPosCallback(window, cursor_pos_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetCharCallback(window, char_callback);
	glfwSetDropCallback(window, drop_callback);

	cgu::glsl_frag_depth_source = glsl_depth_env;

	auto frame_duration = 8ms;
	auto time_next_frame = chrono::steady_clock::now();

	glDisable(GL_DITHER);

	// loop until the user closes the window
	while (!glfwWindowShouldClose(window)) {

		// frame rate limit
		this_thread::sleep_until(time_next_frame);
		time_next_frame += frame_duration;

		glfwPollEvents();
		
		//if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, true);

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		if (ImGui::Begin("Models")) {
			ImGui::Text("Drag'n'Drop to load");
		}
		ImGui::End();

		render(window);

		ImGui::ShowMetricsWindow();

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	glfwTerminate();

}

namespace {

	void drop_callback(GLFWwindow *win, int count, const char **paths) {
		for (int i = 0; i < count; i++) {
			auto e = std::make_unique<ModelEntity>();
			e->load(std::filesystem::u8path(paths[i]));
			entities.push_back(std::move(e));
		}
	}

	void cursor_pos_callback(GLFWwindow *win, double xpos, double ypos) {
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) return;
		glm::ivec2 winsize;
		glfwGetWindowSize(win, &winsize.x, &winsize.y);
		cam.update(glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS, {float(xpos), float(ypos)}, winsize);
	}

	void mouse_button_callback(GLFWwindow *win, int button, int action, int mods) {
		ImGui_ImplGlfw_MouseButtonCallback(win, button, action, mods);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) return;
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			cur_sel.select_entity = cur_sel.hover_entity;
			cur_sel.select_vertex = cur_sel.hover_vertex;
		}
	}

	void scroll_callback(GLFWwindow *win, double xoffset, double yoffset) {
		ImGui_ImplGlfw_ScrollCallback(win, xoffset, yoffset);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) return;
		cam.zoom(float(yoffset));
	}

	void key_callback(GLFWwindow *win, int key, int scancode, int action, int mods) {
		ImGui_ImplGlfw_KeyCallback(win, key, scancode, action, mods);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureKeyboard) return;
	}

	void char_callback(GLFWwindow *win, unsigned int c) {
		ImGui_ImplGlfw_CharCallback(win, c);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantTextInput) return;
	}

	void APIENTRY gl_debug_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei, const GLchar *message, const GLvoid *) {
		// don't report notification messages
		if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;

		// nvidia: avoid debug spam about attribute offsets
		if (id == 131076) return;

		cerr << "GL [" << cgu::gl_debug_source_string(source) << "] " << cgu::gl_debug_type_string(type) << ' ' << id << " : ";
		cerr << message << " (Severity: " << cgu::gl_debug_severity_string(severity) << ')' << endl;

		if (type == GL_DEBUG_TYPE_ERROR_ARB) throw runtime_error("GL Error: "s + message);
	}

}
