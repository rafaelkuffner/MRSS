
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

	ImGuiContext *imguictx;

	float zfar = 1000;
	cgu::orbital_camera cam;
	
	cgu::gl_framebuffer fb_scene{
		cgu::gl_rendertarget_params{GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT, GL_FLOAT},
		cgu::gl_rendertarget_params{GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, GL_RGBA16F, GL_RGBA, GL_FLOAT}
	};

	std::vector<std::unique_ptr<Entity>> entities;
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

	GLFWwindow *window = glfwCreateWindow(800, 600, "CGU", nullptr, nullptr);
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

	auto frame_duration = 8ms;
	auto time_next_frame = chrono::steady_clock::now();

	// loop until the user closes the window
	while (!glfwWindowShouldClose(window)) {

		// frame rate limit
		this_thread::sleep_until(time_next_frame);
		time_next_frame += frame_duration;

		glfwPollEvents();
		
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, true);

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		if (ImGui::Begin("Models")) {
			ImGui::Text("Drag'n'Drop to load");
		}
		ImGui::End();

		glm::ivec2 fbsize;
		glfwGetFramebufferSize(window, &fbsize.x, &fbsize.y);
		glViewport(0, 0, fbsize.x, fbsize.y);

		fb_scene.bind(GL_DRAW_FRAMEBUFFER, fbsize);

		glClearColor(0.1f, 0.1f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 proj = glm::perspective(1.f, float(fbsize.x) / fbsize.y, 0.1f, zfar);
		glm::mat4 view = cam.view();

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glPointSize(3);

		for (auto &e : entities) {
			e->draw(view, proj, zfar);
		}

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
		cam.update(glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS, {float(xpos), float(ypos)}, winsize);
	}

	void mouse_button_callback(GLFWwindow *win, int button, int action, int mods) {
		ImGui_ImplGlfw_MouseButtonCallback(win, button, action, mods);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) return;
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
