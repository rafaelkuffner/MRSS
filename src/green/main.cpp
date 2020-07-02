﻿
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <omp.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <glm/gtc/matrix_transform.hpp>

#include <cgu/opengl.hpp>
#include <cgu/util.hpp>
#include <cgu/mesh.hpp>
#include <cgu/shader.hpp>

#include "main.hpp"
#include "imguiex.hpp"
#include "dialog.hpp"
#include "model.hpp"
#include "saliency.hpp"
#include "clipp.h"
#include "uilocale.hpp"

#include "deferred.glsl.hpp"

#include "gitver.hpp"
#include "about_licences.txt.hpp"

using namespace std;
using namespace green;

namespace {
	void drop_callback(GLFWwindow *win, int count, const char **paths);
	void cursor_pos_callback(GLFWwindow *win, double xpos, double ypos);
	void mouse_button_callback(GLFWwindow *win, int button, int action, int mods);
	void scroll_callback(GLFWwindow *win, double xoffset, double yoffset);
	void key_callback(GLFWwindow *win, int key, int scancode, int action, int mods);
	void char_callback(GLFWwindow *win, unsigned int c);
	void focus_callback(GLFWwindow *win, int focused);
	void APIENTRY gl_debug_callback(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar *, const GLvoid *);
	void set_ui_thread_priority();
	void init_console();
	void init_fonts(ImGuiIO &io, ImFontConfig &fc);
	void load_model(const std::filesystem::path &p);

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

	GLFWwindow *window = nullptr;
	ImGuiContext *imguictx = nullptr;

	float zfar = 1000000;
	cgu::orbital_camera cam;
	float cam_fov = 1.f;
	
	cgu::gl_framebuffer fb_scene{
		cgu::gl_rendertarget_params{GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT, GL_FLOAT},
		cgu::gl_rendertarget_params{GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, GL_RGBA16F, GL_RGBA, GL_FLOAT},
		cgu::gl_rendertarget_params{GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, GL_RG32I, GL_RG_INTEGER, GL_INT}
	};

	std::vector<std::unique_ptr<Entity>> entities;
	ModelEntity *select_model = nullptr;
	ModelEntity *hover_model = nullptr;

	std::vector<std::function<void(bool *p_open)>> persistent_ui;

	std::future<std::vector<std::filesystem::path>> open_paths_future;
	std::future<std::filesystem::path> save_path_future;

	enum class drag_mode {
		plane, axis
	};

	GLsync sync_read_ids = nullptr;
	cgu::gl_object buf_read_ids, buf_read_depths;
	glm::ivec2 pos_read_ids{0};
	const int size_read_ids = 64;

	entity_selection cur_sel;
	float cur_depth = zfar;
	glm::vec3 cur_world_pos{0}, cur_camdrag_pos{0};
	glm::vec3 drag_world_origin{0}, drag_world_pos{0};
	chrono::steady_clock::time_point time_drag_start = chrono::steady_clock::now();
	drag_mode cur_drag_mode = drag_mode::plane;

	bool need_select = false;
	bool maybe_dragging = false;
	bool dragging_camera = false;
	bool drag_skip_next = false;
	bool focus_lost = false;
	bool focus_gained = false;
	bool show_grid = true;
	bool show_axes = true;
	bool show_focus = false;
	bool about_window_open = false;
	bool about_licences_window_open = false;
	bool controlhelp_window_open = true;
	int fps = 0;

	bool saliency_window_open = false;
	int sal_entity_id = -1;
	saliency_user_params sal_uparams;
	// TODO saliency progress object could (should?) be moved into the model entity
	// TODO and then progress reporting should be done as with decimation
	saliency_progress sal_progress;
	std::future<saliency_result> sal_future;
	bool sal_need_preview = false;
	float sal_preview_spn = 10;

	bool decimation_window_open = false;
	decimate_user_params dec_uparams;

	float decode_depth(float depth_p) {
		const float c = 0.01;
		float fc = 1.0 / log(zfar * c + 1.0);
		return (exp(depth_p / fc) - 1.0) / c;
	}

	struct unproject_result {
		glm::vec3 origin{0};
		glm::vec3 hitpos{0};
		glm::vec3 dir{0};

		friend unproject_result operator*(const glm::mat4 &m, const unproject_result &r) {
			unproject_result rr;
			rr.origin = glm::vec3(m * glm::vec4(r.origin, 1));
			rr.hitpos = glm::vec3(m * glm::vec4(r.hitpos, 1));
			rr.dir = glm::vec3(m * glm::vec4(r.dir, 0));
			return rr;
		}
	};

	unproject_result unproject(const glm::mat4 &proj_inv, const glm::vec2 &pos_p, float depth_v) {
		// near plane position
		glm::vec3 pos0_p = glm::vec3(pos_p, -1.0);
		glm::vec4 pos0_vh = proj_inv * glm::vec4(pos0_p, 1.0);
		glm::vec3 pos0_v = glm::vec3(pos0_vh) / pos0_vh.w;
		// not-far plane position
		// TODO better not-far-plane?
		glm::vec3 pos1_p = glm::vec3(pos_p, 0.0);
		glm::vec4 pos1_vh = proj_inv * glm::vec4(pos1_p, 1.0);
		glm::vec3 pos1_v = glm::vec3(pos1_vh) / pos1_vh.w;
		glm::vec3 dir_v = normalize(pos1_v - pos0_v);
		glm::vec3 pos_v = pos0_v + dir_v * ((depth_v + pos0_v.z) / -dir_v.z);
		return {pos0_v, pos_v, dir_v};
	}

	bool is_dragging() {
		return maybe_dragging && (cur_sel.select_entity >= 0 || dragging_camera) && (chrono::steady_clock::now() - time_drag_start > chrono::milliseconds(100));
	}

	float ray_plane_intersect(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, const glm::vec3 &plane_norm, float plane_d) {
		const float d0 = dot(ray_origin, plane_norm);
		const float dd = dot(ray_dir, plane_norm);
		return (plane_d - d0) / dd;
	}

	glm::vec3 drag_plane_normal(const glm::vec3 &raydir, const glm::vec3 &axis, drag_mode mode) {
		return mode == drag_mode::plane ? axis : normalize(cross(cross(raydir, axis), axis));
	}

	void update_dragging(const unproject_result &world_ray, const glm::vec3 &axis, drag_mode mode) {
		const glm::vec3 plane_norm = drag_plane_normal(world_ray.dir, axis, mode);
		const float plane_d = dot(plane_norm, drag_world_origin);
		const float k = ray_plane_intersect(world_ray.origin, world_ray.dir, plane_norm, plane_d);
		if (!isfinite(k) || k < 0.1f || abs(dot(plane_norm, world_ray.dir)) < 0.05f) {
			// -ve k would mean we were dragging on the 'other side' of the drag place (which works but is not nice to use)
			// skip next drag to avoid sudden jumps when returning to the draggable zone
			drag_skip_next = true;
			return;
		}
		if (drag_skip_next) {
			drag_skip_next = false;
			return;
		}
		const auto p = world_ray.origin + k * world_ray.dir;
		auto delta = p - drag_world_pos;
		if (mode == drag_mode::axis) delta = axis * dot(delta, axis);
		drag_world_pos = p;
		if (dragging_camera) {
			cam.focus -= delta;
			drag_world_pos -= delta;
		} else {
			for (auto &e : entities) {
				if (e->id() == cur_sel.select_entity) {
					e->move_by(delta);
				}
			}
		}
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
				cur_sel.hover_entity_dist = 9001;
				cur_sel.hover_vertex_dist = 9001;
				cur_sel.hover_entity = -1;
				cur_sel.hover_vertex = -1;
				// first find hovered vertex
				for (int y = 0; y < size_read_ids; y++) {
					for (int x = 0; x < size_read_ids; x++) {
						glm::ivec2 ids = iddata[size_read_ids * y + x];
						if (ids.x != -1 && ids.y != -1) {
							float dist = glm::length(glm::vec2{pos_read_ids} + glm::vec2{x, y} - glm::vec2{mouse_pos_fb});
							if (dist < cur_sel.hover_vertex_dist) {
								cur_sel.hover_entity_dist = dist;
								cur_sel.hover_vertex_dist = dist;
								cur_sel.hover_entity = ids.x;
								cur_sel.hover_vertex = ids.y;
							}
						}
					}
				}
				// then find hovered entity, constrained to match hovered vertex
				for (int y = 0; y < size_read_ids; y++) {
					for (int x = 0; x < size_read_ids; x++) {
						glm::ivec2 ids = iddata[size_read_ids * y + x];
						if (ids.x != -1 && (ids.x == cur_sel.hover_entity || cur_sel.hover_vertex == -1)) {
							float dist = glm::length(glm::vec2{pos_read_ids} + glm::vec2{x, y} - glm::vec2{mouse_pos_fb});
							if (dist < cur_sel.hover_entity_dist) {
								cur_sel.hover_entity_dist = dist;
								cur_sel.hover_entity = ids.x;
							}
							
						}
					}
				}
				// ensure dist 0 when not hovering anything (nice for camera drag)
				if (cur_sel.hover_entity < 0) cur_sel.hover_entity_dist = 0;
				if (cur_sel.hover_vertex < 0) cur_sel.hover_vertex_dist = 0;
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
				//if (ImGui::Begin("Hover")) {
				//	ImGui::Text("Entity %d; Vertex %d", cur_sel.hover_entity, cur_sel.hover_vertex);
				//	ImGui::Text("Depth %f", cur_depth);
				//}
				//ImGui::End();
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
			// read buffer only specifies color, using depth component will read depth
			glReadBuffer(GL_NONE);
			glReadPixels(pos_read_ids.x, pos_read_ids.y, size_read_ids, size_read_ids, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
			sync_read_ids = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
			glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
			glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
		}
	}

	void draw_window_saliency() {
		using namespace ImGui;
		if (!saliency_window_open) return;
		const auto winsize = GetIO().DisplaySize;
		SetNextWindowBgAlpha(0.5f);
		SetNextWindowSize({500, 500}, ImGuiCond_Appearing);
		SetNextWindowPos({winsize.x / 2.f, winsize.y / 2.f}, ImGuiCond_Appearing, {0.5f, 0.5f});
		if (Begin("Saliency", &saliency_window_open)) {
			int eid = sal_entity_id >= 0 ? sal_entity_id : cur_sel.select_entity;
			Entity *ep = nullptr;
			for (auto &e : entities) {
				if (e->id() == eid) ep = e.get();
			}
			if (!ep) sal_entity_id = -1;
			if (ep) {
				if (Button("Reselect")) {
					cur_sel.select_entity = eid;
					cur_sel.select_vertex = -1;
				}
				SameLine();
				Text("%s", ep->name().c_str());
			}
			sal_need_preview |= Checkbox("Preview", &sal_uparams.preview);
			SetHoveredTooltip("Enable interactive preview\nUse on small models only! (< ~500k vertices)\nCoarse subsampling will be activated.");
			SameLine();
			bool go = false;
			if (sal_future.valid()) {
				if (Button("Cancel", {-1, 0})) sal_progress.should_cancel = true;
			} else {
				if (ep) {
					if (sal_uparams.preview) {
						ButtonDisabled("GO", {-1, 0});
					} else {
						if (Button("GO", {-1, 0})) go = true;
					}
				} else {
					TextDisabled("Select a model");
				}
			}
			Separator();
			sal_need_preview |= edit_saliency_params(sal_uparams);
			Separator();
			draw_saliency_progress(sal_progress);
			if (Button("Clear", {-1, 0})) sal_progress = {};
			if (sal_need_preview && sal_uparams.preview) {
				// user params edited and preview enabled, try launch preview
				if (sal_future.valid()) {
					sal_progress.should_cancel = true;
				} else {
					go = true;
					sal_need_preview = false;
				}
			}
			if (go && ep) {
				// launc calculation
				auto real_uparams = sal_uparams;
				if (real_uparams.preview) {
					real_uparams.subsample_auto = true;
					real_uparams.subsample_manual = false;
					real_uparams.samples_per_neighborhood = sal_preview_spn;
				}
				sal_entity_id = eid;
				sal_progress = {};
				sal_progress.levels.resize(real_uparams.levels);
				sal_future = ep->compute_saliency_async(real_uparams, sal_progress);
			}
		}
		End();
	}

	void draw_window_decimation() {
		using namespace ImGui;
		if (!decimation_window_open) return;
		const auto winsize = GetIO().DisplaySize;
		SetNextWindowBgAlpha(0.5f);
		SetNextWindowSize({500, 500}, ImGuiCond_Appearing);
		SetNextWindowPos({winsize.x / 2.f, winsize.y / 2.f}, ImGuiCond_Appearing, {0.5f, 0.5f});
		if (Begin("Decimation", &decimation_window_open)) {
			PushStyleColor(ImGuiCol_Text, {0.9f, 0.4f, 0.4f, 1});
			Text("-- Work in Progress --");
			PopStyleColor();
			if (select_model) {
				Text("%s", select_model->name().c_str());
				if (select_model->decimated()) {
					PushStyleColor(ImGuiCol_Text, {0.9f, 0.4f, 0.4f, 1});
					Text("Warning: model has already been decimated");
					PopStyleColor();
				}
				auto *sd = select_model->selected_saliency();
				if (sd) Text("Saliency: %s", sd->str().c_str());
				bool go = false;
				if (sd && dec_uparams.use_saliency) {
					if (Button("GO with saliency", {-1, 0})) go = true;
				} else if (!dec_uparams.use_saliency) {
					if (Button("GO without saliency", {-1, 0})) go = true;
				} else {
					TextDisabled("Select a saliency result or uncheck 'use saliency'");
				}
				if (go) {
					auto e = select_model->decimate_async(dec_uparams);
					if (e) entities.push_back(std::move(e));
				}
			} else {
				TextDisabled("Select a model");
			}
			Separator();
			edit_decimate_params(dec_uparams);
			Separator();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
		}
	}

	void render_early_ui() {
		using namespace ImGui;
		SetNextWindowPos({10, 30}, ImGuiCond_FirstUseEver);
		SetNextWindowSize({350, 250}, ImGuiCond_FirstUseEver);
		if (Begin("Models", nullptr, ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
			// ensure the window exists			
		}
		End();

		SetNextWindowPos({10, 290}, ImGuiCond_FirstUseEver);
		SetNextWindowSize({350, 500}, ImGuiCond_FirstUseEver);
		if (Begin("Selection", nullptr, ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
			// ensure the window exists
		}
		End();

		// draw main decimation window before models add progress to it
		draw_window_decimation();
	}

	void render_main_ui() {

		using namespace ImGui;
		const auto winsize = GetIO().DisplaySize;

		// check for saliency completion
		if (sal_future.valid()) {
			if (sal_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				sal_entity_id = -1;
				// getting the result will trigger cleanup (from its dtor)
				if (!sal_future.get()) {
					// cancelled
				}
			}
		}

		// check for files to open
		if (open_paths_future.valid()) {
			if (open_paths_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				auto paths = open_paths_future.get();
				for (auto &p : paths) {
					load_model(p);
				}
			}
		}

		// main menu bar
		const ImVec4 menubg{0.15f, 0.15f, 0.15f, 1};
		const auto normal_frame_padding = GetStyle().FramePadding;
		const auto normal_item_spacing = GetStyle().ItemSpacing;
		PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 10));
		PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6, 6));
		PushStyleColor(ImGuiCol_MenuBarBg, menubg);
		PushStyleColor(ImGuiCol_PopupBg, menubg);
		if (BeginMainMenuBar()) {
			SetCursorPosX(GetCursorPosX() + 30);
			Separator();
			SetNextWindowSizeConstraints({200, 0}, {9001, 9001});
			if (BeginMenu("File")) {
				PushStyleVar(ImGuiStyleVar_FramePadding, normal_frame_padding);
				if (Selectable("Open...")) {
					open_paths_future = open_file_dialog(window, "", true);
				}
				Separator();
				if (select_model) {
					if (Selectable(("Close " + select_model->name() + "...").c_str())) {
						select_model->try_kill();
					}
					if (Selectable(("Export " + select_model->name() + "...").c_str())) {
						select_model->try_export();
					}
				} else {
					TextDisabled("Close");
					TextDisabled("Export");
				}
				Separator();
				PopStyleVar();
				EndMenu();
			}
			if (BeginMenu("View")) {
				PushStyleVar(ImGuiStyleVar_FramePadding, normal_frame_padding);
				Checkbox("Show Grid", &show_grid);
				Checkbox("Show Axes", &show_axes);
				Checkbox("Show Focus", &show_focus);
				Separator();
				PushStyleVar(ImGuiStyleVar_ItemSpacing, normal_item_spacing);
				Text("Camera");
				InputFloat3("Focus", value_ptr(cam.focus));
				SliderAngle("Yaw", &cam.cam_yaw, -180, 180);
				SliderAngle("Pitch", &cam.cam_pitch, -90, 90);
				SliderFloat("Distance", &cam.cam_distance, 0, 10);
				SliderAngle("Vertical FoV", &cam_fov, 0, 170);
				PopStyleVar(2);
				Separator();
				EndMenu();
			}
			if (BeginMenu("Window")) {
				PushStyleVar(ImGuiStyleVar_FramePadding, normal_frame_padding);
				Checkbox("Saliency", &saliency_window_open);
				Checkbox("Decimation", &decimation_window_open);
				Checkbox("Control Help", &controlhelp_window_open);
				Separator();
				PopStyleVar();
				EndMenu();
			}
			if (BeginMenu("Help")) {
				PushStyleVar(ImGuiStyleVar_FramePadding, normal_frame_padding);
				if (Selectable("About")) about_window_open = true;
				TextDisabled("Run with '-h' for command line help");
				Separator();
				PopStyleVar();
				EndMenu();
			}
			// TODO maybe an actual status bar?
			Separator();
			TextDisabled("%3d FPS", fps);
			Separator();
			if (focus_lost) {
				// doing this with a modal popup was kinda nice,
				// but has the side effect of killing other popups
				TextDisabled("Click to regain focus");
				Separator();
			}
			EndMainMenuBar();
		}
		PopStyleColor(2);
		PopStyleVar(2);

		// NOTE BeginMenu respects constraints but doesnt consume them (bug?)
		SetNextWindowSizeConstraints({0, 0}, {9001, 9001});

		draw_window_saliency();

		if (about_window_open) {
			SetNextWindowPos({winsize.x / 2, winsize.y / 2}, ImGuiCond_Appearing, {0.5f, 0.5f});
			if (Begin("About", &about_window_open, ImGuiWindowFlags_AlwaysAutoResize)) {
				Text("Multi-Resolution Subsampled Saliency");
				Separator();
				Text("Copyright 2020\nVictoria University of Wellington\nComputational Media Innovation Centre\nAll rights reserved.");
				Separator();
				Text("Version:\n%s (%s)\n%s", git_describe(), git_revision(), git_timestamp());
				Separator();
				if (Button("Library Licenses")) about_licences_window_open = true;
			}
			End();
		}

		if (about_licences_window_open) {
			SetNextWindowPos({winsize.x / 2, winsize.y / 2}, ImGuiCond_Appearing, {0.5f, 0.5f});
			SetNextWindowSize({winsize.x / 2, 2 * winsize.y / 3}, ImGuiCond_Appearing);
			if (Begin("Library Licences", &about_licences_window_open)) {
				Text("This program uses libraries released under the following licenses:");
				Separator();
				// imgui wants a non-const pointer (because its a text edit field, just in readonly mode)
				// +3 to skip BOM (because VS likes to use it for utf8 files)
				static std::string text = cgu::strings::about_licences + 3;
				InputTextMultiline("", text.data(), text.size(), {-1, -1}, ImGuiInputTextFlags_ReadOnly);
			}
			End();
		}

		if (controlhelp_window_open) {
			SetNextWindowPos({winsize.x - 20, 30}, ImGuiCond_Always, {1.f, 0.f});
			if (Begin("ControlHelp", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoFocusOnAppearing)) {
				PushStyleVar(ImGuiStyleVar_Alpha, 0.5f);
				Text("Left Mouse: Move\n [Object, Horizontal]\n + Shift: Vertical\n + Alt: Camera");
				Text("Right Mouse: Rotate camera");
				Text("Scroll: Zoom");
				PopStyleVar();
			}
			End();
		}

		if (Begin("Models")) {
			TextDisabled("Open with [File > Open] or drag-and-drop");
		}
		End();

	}

	void render_deferred() {

		static GLuint prog = 0;
		if (!prog) {
			prog = cgu::make_shader_program(
				"green::deferred",
				"330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{cgu::glsl_fullscreen_source, cgu::strings::glsl_green_deferred}
			).release();
		}

		glUseProgram(prog);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, fb_scene[GL_COLOR_ATTACHMENT0].tex);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, fb_scene[GL_DEPTH_ATTACHMENT].tex);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, fb_scene[GL_COLOR_ATTACHMENT1].tex);

		glUniform1i(glGetUniformLocation(prog, "u_sampler_color"), 0);
		glUniform1i(glGetUniformLocation(prog, "u_sampler_depth"), 1);
		glUniform1i(glGetUniformLocation(prog, "u_sampler_id"), 2);
		glUniform4iv(glGetUniformLocation(prog, "u_selection"), 1, (GLint *) &cur_sel);
		glm::ivec4 show_sel{hover_model && hover_model->show_entity_outline(), 1, select_model && select_model->show_entity_outline(), 1};
		glUniform4iv(glGetUniformLocation(prog, "u_show_selection"), 1, value_ptr(show_sel));

		cgu::draw_dummy();

	}

	void render() {

		// TODO dont redraw models when nothing has changed

		glm::ivec2 fbsize;
		glfwGetFramebufferSize(window, &fbsize.x, &fbsize.y);
		glViewport(0, 0, fbsize.x, fbsize.y);

		// FIXME when fb size != win size (hidpi)
		glm::dvec2 mouse_pos{0};
		glfwGetCursorPos(window, &mouse_pos.x, &mouse_pos.y);
		glm::ivec2 mouse_pos_fb = {int(mouse_pos.x), fbsize.y - int(mouse_pos.y) - 1};

		fb_scene.bind(GL_DRAW_FRAMEBUFFER, fbsize);

		if (glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
			// nvidia: fbo is incomplete when window is minimized
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
			return;
		}

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glClearBufferfv(GL_DEPTH, 0, value_ptr(glm::vec4{1.f}));
		glClearBufferfv(GL_COLOR, 0, value_ptr(glm::vec4{0.1f, 0.1f, 0.2f, 1.0f}));
		glClearBufferiv(GL_COLOR, 1, value_ptr(glm::ivec4{-1}));
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glm::mat4 proj = glm::perspective(cam_fov, float(fbsize.x) / fbsize.y, 0.1f, zfar);
		glm::mat4 view = cam.view();
		auto proj_inv = inverse(proj);
		auto view_inv = inverse(view);
		auto world_ray = view_inv * unproject(proj_inv, glm::vec2(mouse_pos_fb) / glm::vec2(fbsize) * 2.f - 1.f, cur_depth);
		cur_world_pos = world_ray.hitpos;
		// TODO dragging on cam.focus (instead of y=0) is more generic but unintuitive to use
		// forcing y=0 assumes the plane is y=0 or parallel to the y axis
		// NOTE when dragging vertically, the point being dragged is on the focus plane not the ground plane,
		// so the apparent movement of the ground plane can seem incorrect (in magnitude)
		auto camdrag_norm = drag_plane_normal(world_ray.dir, {0, 1, 0}, cur_drag_mode);
		cur_camdrag_pos = world_ray.origin + world_ray.dir * ray_plane_intersect(world_ray.origin, world_ray.dir, camdrag_norm, dot(camdrag_norm, {cam.focus.x, 0, cam.focus.z}));

		if (need_select) {
			// select entity/vertex
			// NOTE the event handler runs before we've updated the selected entity
			// so we delay actually doing the selection until here
			need_select = false;
			cur_sel.select_entity = cur_sel.hover_entity;
			cur_sel.select_vertex = cur_sel.hover_vertex;
			if (cur_sel.hover_entity_dist == 0) {
				// begin maybe dragging; dist 0 to ensure click was actually on the entity,
				// because we rely on the click position on the entity for dragging
				maybe_dragging = true;				
				time_drag_start = chrono::steady_clock::now();
				// need to delay setting origin until cur_camdrag_pos is set with the correct drag mode!
				if (dragging_camera) {
					// drag camera
					drag_world_origin = cur_camdrag_pos;
					drag_world_pos = cur_camdrag_pos;
				} else {
					// drag entity
					drag_world_origin = cur_world_pos;
					drag_world_pos = cur_world_pos;
				}
			}
		}

		if (is_dragging()) {
			update_dragging(world_ray, {0, 1, 0}, cur_drag_mode);
		}

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glPointSize(1);

		select_model = nullptr;
		hover_model = nullptr;
		bool can_hover = false;
		for (auto it = entities.begin(); it != entities.end(); ) {
			can_hover = true;
			auto &e = *it;
			e->draw(view, proj, zfar);
			if (e->dead()) {
				it = entities.erase(it);
				if (static_cast<Entity *>(select_model) == e.get()) select_model = nullptr;
				if (static_cast<Entity *>(hover_model) == e.get()) hover_model = nullptr;
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
		//cgu::draw_texture2d(fb_scene[GL_COLOR_ATTACHMENT0].tex, fb_scene[GL_DEPTH_ATTACHMENT].tex, 0);
		render_deferred();

		glDisable(GL_FRAMEBUFFER_SRGB);

		glDepthFunc(GL_LEQUAL);

		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (show_grid) cgu::draw_grid(view, proj, zfar);
		if (show_axes) cgu::draw_axes(view, proj, zfar);
		if (show_focus) {
			cgu::draw_axes_params axp;
			for (int i = 0; i < 3; i++) {
				axp.color_pos[i].a = 0.5f;
				axp.color_neg[i].a = 0.5f;
			}
			axp.axislength = 1;
			auto m = view;
			m = glm::translate(m, cam.focus);
			cgu::draw_axes(m, proj, zfar, axp);
		}
		glDisable(GL_BLEND);

		glDisable(GL_DEPTH_TEST);

		// clearing the id buffer may not actually execute if we dont draw anything?! (driver bug?)
		if (can_hover) update_hover(mouse_pos_fb);

	}

	void center_window(GLFWwindow *window, GLFWmonitor *monitor) {
		if (!monitor) return;

		const GLFWvidmode *mode = glfwGetVideoMode(monitor);
		if (!mode) return;

		int monitorX, monitorY;
		glfwGetMonitorPos(monitor, &monitorX, &monitorY);

		int windowWidth, windowHeight;
		glfwGetWindowSize(window, &windowWidth, &windowHeight);

		glfwSetWindowPos(
			window,
			monitorX + (mode->width - windowWidth) / 2,
			monitorY + (mode->height - windowHeight) / 2
		);
	}

	void main_gui() {

		if (!glfwInit()) {
			cerr << "Error: Could not initialize GLFW" << endl;
			abort();
		}

		set_ui_thread_priority();

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

		window = glfwCreateWindow(1280, 800, (std::string("Multi-Resolution Subsampled Saliency (") + git_describe() + ")").c_str(), nullptr, nullptr);
		if (!window) {
			cerr << "Error: Could not create GLFW window" << endl;
			abort();
		}

		center_window(window, glfwGetPrimaryMonitor());

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

		{
			ImGuiStyle &sty = ImGui::GetStyle();
			//sty.FramePadding.x += 1;
			//sty.FramePadding.y += 1;
			sty.SelectableTextAlign.y = 0.5f;
			sty.Colors[ImGuiCol_WindowBg].w = 0.6f;
			sty.Colors[ImGuiCol_ScrollbarBg].w = 0.6f;
			sty.Colors[ImGuiCol_ScrollbarGrab].w = 0.6f;
		}

		{
			ImGuiIO& io = ImGui::GetIO();
			// disable imgui.ini
			io.IniFilename = nullptr;
			// fonts
			io.Fonts->Flags |= ImFontAtlasFlags_NoPowerOfTwoHeight;
			ImFontConfig fc;
			// dont oversample to save texture space
			fc.PixelSnapH = true;
			fc.OversampleH = 1;
			// use imgui default font for ascii
			io.Fonts->AddFontDefault(&fc);
			// load other fonts
			fc.MergeMode = true;
			init_fonts(io, fc);
			// NOTE we are not able to achieve full unicode rendering because
			// a) iirc imgui only supports the BMP, and is certainly not capable of fancy layout
			// b) the font texture atlas would be intractably large
			io.Fonts->Build();
			// hacks to increase line height to make japanese etc text easier
			io.Fonts->Fonts[0]->FontSize += 3;
			io.Fonts->Fonts[0]->DisplayOffset.y += 1;
			std::cout << "imgui font atlas " << io.Fonts->TexWidth << "x" << io.Fonts->TexHeight << std::endl;
		}

		// attach input callbacks to window
		glfwSetCursorPosCallback(window, cursor_pos_callback);
		glfwSetMouseButtonCallback(window, mouse_button_callback);
		glfwSetScrollCallback(window, scroll_callback);
		glfwSetKeyCallback(window, key_callback);
		glfwSetCharCallback(window, char_callback);
		glfwSetDropCallback(window, drop_callback);
		glfwSetWindowFocusCallback(window, focus_callback);

		cgu::glsl_frag_depth_source = glsl_depth_env;

		auto time_next_frame = chrono::steady_clock::now();
		auto time_last_fps = time_next_frame;
		int fps_counter = 0;

		cam.focus.y = 1;
		cam.cam_yaw = -glm::pi<float>() / 4;
		cam.cam_pitch = glm::pi<float>() / 8;

		glDisable(GL_DITHER);

		const auto poll_events = []() {
			glfwPollEvents();
			if (focus_gained) {
				focus_lost = false;
				focus_gained = false;
			}
		};

		// loop until the user closes the window
		while (!glfwWindowShouldClose(window)) {

			poll_events();

			// frame rate limiter and counter
			auto now = chrono::steady_clock::now();
			while (time_next_frame - now >= 1ms) {
				this_thread::sleep_for(1ms);
				// keep polling events when limiting frame rate to ensure responsiveness with e.g. win32 dialogs
				poll_events();
				now = chrono::steady_clock::now();
			}
			const auto frame_duration = focus_lost ? 100ms : 6ms;
			time_next_frame = max(time_next_frame, now - frame_duration) + frame_duration;
			fps_counter++;
			if (now - time_last_fps >= 1s) {
				fps = fps_counter;
				fps_counter = 0;
				time_last_fps = now;
			}

			//if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, true);

			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();

			// render scene before main ui so that current model can be determined
			render_early_ui();
			render();
			render_main_ui();

			for (auto it = persistent_ui.begin(); it != persistent_ui.end(); ) {
				bool r = true;
				(*it)(&r);
				if (r) {
					++it;
				} else {
					it = persistent_ui.erase(it);
				}
			}

			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

			glfwSwapBuffers(window);
		}

		glfwTerminate();
		sal_progress.should_cancel = true;
		if (sal_future.valid()) sal_future.wait();

	}

	// expects utf-8 encoded args
	void main_cli(int argc, const char *argv[]) {
		using namespace clipp;
		using namespace uistrings;

		const uilocale &loc = uilocale_en();

		string infile, outfile;
		string decprop = "quality";
		bool show_gui = false;
		bool do_version = false, do_help = false, do_sal = false, do_dec = false;
		bool save_ascii = false;
		int threads = 0;

		auto alt_opts = group{
			option("--version").set(do_version)
				.doc(loc[help_cli_version].clone()),
			option("-h", "--help").set(do_help)
				.doc(loc[help_cli_help].clone())
		};
		
		auto common_opts = group{
			(required("-i", "--input") & value("infile", infile))
				.doc(loc[help_input].clone()),
			(option("-o", "--output") & value("outfile", outfile))
				.doc(loc[help_output].clone()),
			(option("-j", "--threads") & integer("threads", threads))
				.doc(loc[help_threads].clone()),
			option("--ascii").set(save_ascii)
				.doc(loc[help_cli_ascii].clone()),
			option("--gui").set(show_gui)
				.doc(loc[help_cli_gui].clone()),
			option("--noprogress").call([]{
					sal_uparams.show_progress = false;
					// TODO dec progress
				})
				.doc(loc[help_cli_noprogress].clone())
		}.doc(loc[help_cli_opts_common].clone());

		auto sal_opts = group{
			option("--saliency").set(do_sal)
				.doc(loc[help_sal_go].clone()),
			(option("-a", "--area") & number("area", sal_uparams.area))
				.doc(loc[help_sal_area].clone()),
			(option("-l", "--levels") & integer("levels", sal_uparams.levels))
				.doc(loc[help_sal_levels].clone()),
			(option("-r", "--contrast") & number("contrast", sal_uparams.normal_power))
				.doc(loc[help_sal_normpower].clone()),
			(option("-c", "--contour") & number("contour", sal_uparams.curv_weight))
				.doc(loc[help_sal_curvweight].clone()),
			(option("-n", "--noisefilter").set(sal_uparams.normalmap_filter))
				.doc(loc[help_sal_noisefilter].clone()),
			(option("-e", "--noiseheight") & number("height", sal_uparams.noise_height))
				.doc(loc[help_sal_noiseheight].clone()),
			(option("-s", "--subsampling") & number("samples-per-neighborhood", sal_uparams.samples_per_neighborhood))
				.doc(loc[help_sal_samplespern].clone()),
			option("--fullsampling").call([]{ sal_uparams.subsample_auto = false; })
				.doc(loc[help_sal_full].clone())
		}.doc(loc[help_cli_opts_saliency].clone());

		auto dec_opts = group{
			option("--decimate").set(do_dec)
				.doc(loc[help_dec_go].clone()),
			(option("-t", "--targetverts") & integer("verts", dec_uparams.targetverts))
				.doc(loc[help_dec_targetverts].clone()),
			(option("-w", "--salweight") & number("weight", dec_uparams.weight))
				.doc(loc[help_dec_weight].clone()),
			(option("-p", "--binpower") & number("power", dec_uparams.power))
				.doc(loc[help_dec_power].clone()),
			(option("--decbins") & integer("bins", dec_uparams.nbins))
				.doc(loc[help_dec_bins].clone()),
			(option("--decprop") & value("propname", decprop))
				.doc(loc[help_dec_propname].clone()),
			option("--decimate-nosaliency").call([]{ dec_uparams.use_saliency = false; })
				.doc(loc[help_dec_nosaliency].clone())
		}.doc(loc[help_cli_opts_decimate].clone());

		auto opts = (common_opts, sal_opts, dec_opts) | alt_opts;

		if (auto res = parse(argv + 1, argv + argc, opts); !res) {
			for (auto &m : res.missing()) {
				if (m.param()->flags().size()) {
					cerr << "missing required argument " << m.param()->flags()[0] << endl;
				} else {
					cerr << "missing required argument <" << m.param()->label() << "> after '" << argv[1 + m.after_index()] << "'" << endl;
				}
			}
			for (auto &m : res) {
				if (m.any_error()) {
					cerr << "argument '" << m.arg() << "' : ";
					if (m.blocked()) cerr << "blocked";
					if (m.conflict()) cerr << "conflict";
					cerr << endl;
				}
			}
			cerr << "usage:\n" << usage_lines(opts) << endl;
			exit(1);
		} 
		
		if (do_version) {
			cout << "Multi-Resolution Subsampled Saliency\n" << git_describe() << " (" << git_revision() << ")" << endl;
		}

		if (do_help) {
			doc_formatting fmt;
			fmt.first_column(4);
			const auto exename = std::filesystem::u8path(argv[0]).filename().replace_extension("").u8string();
			auto man = make_man_page(opts, exename, fmt);
			man.prepend_section(
				"DESCRIPTION",
				"    Multi-Resolution Subsampled Saliency\n"
				"    Copyright 2020\n    Victoria University of Wellington\n"
				"    Computational Media Innovation Centre\n    All rights reserved."
			);
			man.append_section(
				"VERSION", 
				string("    ") + git_describe() + " (" + git_revision() + ")\n    " + git_timestamp()
			);
			cout << man << endl;
		}
		
		if (!(do_sal || do_dec)) {
			if (!do_help && !do_version) cerr << "nothing to do" << endl;
			exit(0);
		}

		sal_uparams.thread_count = threads;
		sal_uparams.sanitize();
		dec_uparams.sanitize();
		sal_progress.levels.resize(sal_uparams.levels);

		const auto inpath = std::filesystem::u8path(infile);

		Model m;
		Model md;
		Model *pmr = &m;

		try {
			m = Model{inpath};
		} catch (exception &e) {
			cerr << "failed to load model: " << e.what() << endl;
			exit(1);
		}

		decimate_progress dec_progress;
		model_saliency_data sd;

		if (!do_sal && do_dec && dec_uparams.use_saliency) {
			if (m.trimesh().get_property_handle(sd.prop_saliency, decprop)) {
				sd.filename = inpath.filename().u8string();
				sd.propname = decprop;
			} else {
				cerr << "couldn't find saliency property '" << decprop << "' for decimation" << endl;
				exit(1);
			}
		}

		if (do_sal) {

			saliency_mesh_params mparams;
			mparams.mesh = &m.trimesh();
			mparams.prop_vertex_area = m.prop_vertex_area();
			mparams.prop_edge_length = m.prop_edge_length();

			// create properties
			mparams.prop_saliency_levels.resize(sal_uparams.levels);
			mparams.mesh->add_property(mparams.prop_curvature);
			mparams.mesh->add_property(mparams.prop_saliency);
			mparams.mesh->add_property(mparams.prop_sampled);
			for (int i = 0; i < sal_uparams.levels; i++) {
				mparams.mesh->add_property(mparams.prop_saliency_levels[i]);
			}

			// calculate saliency
			if (!compute_saliency(mparams, sal_uparams, sal_progress)) {
				cerr << "saliency cancelled" << endl;
				exit(1);
			}

			// destroy temp properties
			mparams.mesh->remove_property(mparams.prop_curvature);
			for (int i = 0; i < sal_uparams.levels; i++) {
				mparams.mesh->remove_property(mparams.prop_saliency_levels[i]);
			}

			sd.prop_saliency = mparams.prop_saliency;
			sd.prop_sampled = mparams.prop_sampled;
			sd.uparams = sal_uparams;
			sd.progress = sal_progress;

			m.original_saliency().push_back(sd);
		}

		if (do_dec) {
			md = m.prepare_decimate(sd);
			if (!md.decimate(dec_uparams, dec_progress)) {
				cerr << "decimation cancelled" << endl;
				exit(1);
			}
			pmr = &md;
		}

		if (outfile.size()) {
			try {
				pmr->save(std::filesystem::u8path(outfile), sd.prop_saliency, !save_ascii);
			} catch (exception &e) {
				cerr << "failed to save model: " << e.what() << endl;
				if (!show_gui) exit(1);
			}
		}

		if (show_gui) {
			auto e = std::make_unique<ModelEntity>();
			ModelEntity *ep = e.get();
			e->load(inpath, std::move(m));
			entities.push_back(std::move(e));
			if (do_dec) {
				auto ed = std::make_unique<ModelEntity>();
				ed->load(inpath, std::move(md), dec_uparams, dec_progress);
				entities.push_back(std::move(ed));
				// move the pre-decimated result out of the way
				ep->move_by({0, 0, -4});
			}
			main_gui();
		}

	}
}

namespace green {

	entity_selection & ui_selection() {
		return cur_sel;
	}

	saliency_user_params & ui_saliency_user_params() {
		return sal_uparams;
	}

	const saliency_progress & ui_saliency_progress() {
		return sal_progress;
	}

	decimate_user_params & ui_decimate_user_params() {
		return dec_uparams;
	}

	void ui_select_model(ModelEntity *e) {
		select_model = e;
	}

	ModelEntity * ui_select_model() {
		return select_model;
	}

	void ui_hover_model(ModelEntity *e) {
		hover_model = e;
	}

	ModelEntity * ui_hover_model() {
		return hover_model;
	}

	std::future<std::filesystem::path> & ui_save_path(const std::filesystem::path &hint, bool prompt) {
		if (save_path_future.valid()) {
			if (prompt && save_path_future.wait_for(0s) == future_status::ready) {
				// discard existing save path
				save_path_future.get();
			} else {
				return save_path_future;
			}
		}
		if (prompt) save_path_future = save_file_dialog(window, hint);
		return save_path_future;
	}

	void ui_spawn(std::function<void(bool *p_open)> f) {
		persistent_ui.push_back(std::move(f));
	}

	bool ui_saliency_window_open() {
		return saliency_window_open;
	}

	bool ui_decimation_window_open() {
		return decimation_window_open;
	}

	const char * saliency_state_str(saliency_computation_state s) {
		switch (s) {
		case saliency_computation_state::idle:
			return "Idle";
		case saliency_computation_state::curv:
			return "Computing curvature";
		case saliency_computation_state::area:
			return "Computing surface area";
		case saliency_computation_state::nhprep:
			return "Preparing mesh data for neighborhood search";
		case saliency_computation_state::cand:
			return "Preparing subsampling candidates";
		case saliency_computation_state::run_full:
			return "Running (full)";
		case saliency_computation_state::run_sub:
			return "Running (subsampled)";
		case saliency_computation_state::merge:
			return "Merging saliency levels";
		case saliency_computation_state::norm:
			return "Normalizing saliency";
		case saliency_computation_state::done:
			return "Done";
		case saliency_computation_state::cancelled:
			return "Cancelled";
		default:
			return "???";
		}
	}

	const char * decimation_state_str(decimation_state s) {
		switch (s) {
		case decimation_state::idle:
			return "Idle";
		case decimation_state::bins:
			return "Initializing decimation bins";
		case decimation_state::run:
			return "Running";
		case decimation_state::done:
			return "Done";
		case decimation_state::cancelled:
			return "Cancelled";
		default:
			return "???";
		}
	}
}

namespace ImGui {

	void draw_saliency_progress(const green::saliency_progress &progress) {
		for (int i = 0; i < std::min<int>(progress.completed_levels + 1, progress.levels.size()); i++) {
			const auto &level = progress.levels[i];
			Text("%2d/%zu", i + 1, progress.levels.size());
			SameLine();
			char buf[1024];
			const char *normalmap_filter_str = level.normalmap_filter ? "\Noise filter active" : "";
			if (level.subsampled) {
				const float actual_subsampling = progress.total_vertices / float(level.completed_samples);
				std::snprintf(buf, sizeof(buf), "%.0f%% [subsampled %.1fx]", level.completion * 100, actual_subsampling);
			} else {
				std::snprintf(buf, sizeof(buf), "%.0f%% [full]", level.completion * 100);
			}
			ProgressBar(level.completion, {-1, 0}, buf);
			if (IsItemHovered()) SetTooltip("%d samples / %d vertices%s", level.completed_samples, progress.total_vertices, normalmap_filter_str);
		}
		// TODO cancel button here
		Text("%s [%.3fs]", saliency_state_str(progress.state), progress.elapsed_time / std::chrono::duration<double>(1.0));
	}

	bool edit_saliency_params(green::saliency_user_params &uparams) {
		using green::saliency_user_params;
		using namespace green::uistrings;
		const uilocale &loc = uilocale_en();
		const auto uparams0 = uparams;
		// max hardware threads leaving 1 for UI
		static const int defthreads = std::max(1, int(std::thread::hardware_concurrency()) - 1);
		saliency_user_params defparams;
		defparams.thread_count = defthreads;
		if (!uparams.thread_count) uparams.thread_count = defthreads;
		param_widgets widgets{&loc, &defparams, &uparams};
		TextDisabled("Ctrl-click sliders to enter values directly");
		widgets.slider(param_sal_levels, &saliency_user_params::levels, 1, 6);
		SetHoveredTooltip(loc, help_sal_levels);
		widgets.slider(param_sal_area, &saliency_user_params::area, 0.f, 0.05f, "%.5f", 2.f);
		SetHoveredTooltip(loc, help_sal_area);
		widgets.slider(param_sal_curvweight, &saliency_user_params::curv_weight, 0.f, 1.f);
		SetHoveredTooltip(loc, help_sal_curvweight);
		widgets.slider(param_sal_normpower, &saliency_user_params::normal_power, 0.f, 2.f);
		SetHoveredTooltip(loc, help_sal_normpower);
		widgets.checkbox(param_sal_noisefilter, &saliency_user_params::normalmap_filter);
		SetHoveredTooltip(loc, help_sal_noisefilter);
		if (uparams.normalmap_filter) {
			widgets.slider(param_sal_noiseheight, &saliency_user_params::noise_height, 0.f, 0.01f, "%.4f", 2.f);
			SetHoveredTooltip(loc, help_sal_noiseheight);
		}
		if (uparams.preview) {
			Text("Preview mode: subsampling at %.1f S/N", sal_preview_spn);
		} else {
			widgets.checkbox(param_sal_subsample, &saliency_user_params::subsample_auto);
			SetHoveredTooltip(loc, help_sal_subsample);	
			if (uparams.subsample_manual) {
				// no longer exposing this mode of subsampling in the ui
				//widgets.slider("Rate", &saliency_user_params::subsampling_rate, 1.f, 5000.f, "%.1fx", 3.f);
				//SetHoveredTooltip("Subsampling Rate\nMust be tuned for each model.\nHigher: fewer samples, less accurate results.");
			} else if (uparams.subsample_auto) {
				widgets.slider(param_sal_samplespern, &saliency_user_params::samples_per_neighborhood, 1.f, 500.f, "%.1f", 2.f);
				SetHoveredTooltip(loc, help_sal_samplespern);
			}
		}
		widgets.slider(param_threads, &saliency_user_params::thread_count, 1, defthreads);
		SetHoveredTooltip(loc, help_threads);
		// ensure params are valid
		uparams.sanitize();
		return uparams != uparams0;
	}

	void draw_saliency_params(const green::saliency_user_params &uparams) {
		ImGui::Text("Levels: %d", uparams.levels);
		ImGui::Text("Area: %.3f", uparams.area);
		ImGui::Text("Contour: %.3f", uparams.curv_weight);
		ImGui::Text("Contrast: %.3f", uparams.normal_power);
		ImGui::Text("Noise Filter: %s", uparams.normalmap_filter ? "true" : "false");
		ImGui::Text("Noise Height: %f", uparams.noise_height);
		ImGui::Text("Subsampling: %s", uparams.subsample_manual ? "Manual" : uparams.subsample_auto ? "Auto" : "None");
		if (uparams.subsample_manual) {
			ImGui::Text("Subsampling Rate: %.1f", uparams.subsampling_rate);
		} else if (uparams.subsample_auto) {
			ImGui::Text("Samples per Neighbourhood: %.1f", uparams.samples_per_neighborhood);
		}
	}

	void draw_decimate_progress(green::decimate_progress &progress) {
		ProgressBar(float(progress.completed_collapses) / progress.target_collapses);
		SetHoveredTooltip("Decimated %d / %d vertices", progress.completed_collapses, progress.target_collapses);
		if (progress.state < decimation_state::done) {
			if (Button("Cancel")) progress.should_cancel = true;
			SameLine();
		}
		Text("%s [%.3fs]", decimation_state_str(progress.state), progress.elapsed_time / std::chrono::duration<double>(1.0));
	}

	bool edit_decimate_params(green::decimate_user_params &uparams) {
		using green::decimate_user_params;
		using namespace green::uistrings;
		const uilocale &loc = uilocale_en();
		decimate_user_params defparams;
		param_widgets widgets{&loc, &defparams, &uparams};
		TextDisabled("Ctrl-click sliders to enter values directly");
		widgets.checkbox(param_dec_usesaliency, &decimate_user_params::use_saliency);
		SetHoveredTooltip(loc, help_dec_usesaliency);
		widgets.inputint(param_dec_targetverts, &decimate_user_params::targetverts, 100, 1000);
		SetHoveredTooltip(loc, help_dec_targetverts);
		widgets.slider(param_dec_bins, &decimate_user_params::nbins, 1, 10);
		SetHoveredTooltip(loc, help_dec_bins);
		widgets.slider(param_dec_weight, &decimate_user_params::weight, 0.f, 1.f);
		SetHoveredTooltip(loc, help_dec_weight);
		widgets.slider(param_dec_power, &decimate_user_params::power, 0.f, 2.f);
		SetHoveredTooltip(loc, help_dec_power);
		// ensure params are valid
		uparams.sanitize();
		// TODO return true if modified?
		return true;
	}

}

namespace {

	void load_model(const std::filesystem::path &p) {
		auto e = std::make_unique<ModelEntity>();
		e->load(p);
		entities.push_back(std::move(e));
	}

	void drop_callback(GLFWwindow *win, int count, const char **paths) {
		for (int i = 0; i < count; i++) {
			auto p = std::filesystem::u8path(paths[i]);
			if (std::filesystem::is_directory(p)) {
				std::cout << "chdir to " << p << std::endl;
				std::filesystem::current_path(p);
			} else {
				load_model(p);
			}
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
		if (focus_lost) {
			// prevent mouse clicks that refocus the window from doing anything
			focus_lost = false;
			return;
		}
		ImGui_ImplGlfw_MouseButtonCallback(win, button, action, mods);
		// if not captured then foward to application
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) return;
		maybe_dragging = false;
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			need_select = true;
			dragging_camera = mods & GLFW_MOD_ALT;
			cur_drag_mode = (mods & GLFW_MOD_SHIFT) ? drag_mode::axis : drag_mode::plane;
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

	void focus_callback(GLFWwindow *win, int focused) {
		if (focused) {
			focus_gained = true;
		} else {
			focus_lost = true;
		}
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

#ifdef _WIN32
	extern "C" {
		__declspec(dllimport) int __stdcall GetLastError();
		__declspec(dllimport) void * __stdcall GetCurrentThread();
		__declspec(dllimport) int __stdcall SetThreadPriority(void *hThread, int nPriority);
		__declspec(dllimport) bool __stdcall SetConsoleOutputCP(unsigned int wCodePageID);
		__declspec(dllimport) int __stdcall WideCharToMultiByte(unsigned CodePage, int dwFlags, const wchar_t *lpWideCharStr, int cchWideChar, char *lpMultiByteStr, int cbMultiByte, const char *lpDefaultChar, bool *lpUsedDefaultChar);
	}

	void set_ui_thread_priority() {
		constexpr int THREAD_PRIORITY_ABOVE_NORMAL = 1;
		SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
	}

	void init_console() {
		if (!SetConsoleOutputCP(65001)) {
			std::cerr << "failed to set utf8 console codepage" << std::endl;
		}
	}

	void init_fonts(ImGuiIO &io, ImFontConfig &fc) {
		// TODO font names on earlier windows versions
		// japanese
		fc.GlyphOffset = {0, 2};
		fc.RasterizerMultiply = 2.0f;
		io.Fonts->AddFontFromFileTTF("c:\\windows\\fonts\\yugothl.ttc", 18, &fc, io.Fonts->GetGlyphRangesJapanese());
		io.Fonts->AddFontFromFileTTF("c:\\windows\\fonts\\yugothl.ttc", 18, &fc, io.Fonts->GetGlyphRangesChineseSimplifiedCommon());
		// chinese (for characters not covered by japanese)
		fc.GlyphOffset = {0, 1};
		fc.RasterizerMultiply = 1.5f;
		io.Fonts->AddFontFromFileTTF("c:\\windows\\fonts\\msyhl.ttc", 20, &fc, io.Fonts->GetGlyphRangesChineseSimplifiedCommon());
		// korean
		fc.GlyphOffset = {0, 1};
		fc.RasterizerMultiply = 1.5f;
		io.Fonts->AddFontFromFileTTF("c:\\windows\\fonts\\malgunsl.ttf", 18, &fc, io.Fonts->GetGlyphRangesKorean());
	}

#else
	void set_ui_thread_priority() {

	}

	void init_console() {

	}

	void init_fonts(ImGuiIO &io, ImFontConfig &fc) {

	}
#endif

}

#ifdef _WIN32
int wmain(int argc, const wchar_t *wargv[]) {
	init_console();
	if (argc > 1) {
		// need to convert to utf-8
		std::vector<const char *> argv(argc);
		for (int i = 0; i < argc; i++) {
			const int wl = wcsnlen(wargv[i], 32767);
			// pessimistic buffer sizing (note x2 is same size in bytes), +1 for null terminator
			const int cl = wl * 3 + 1;
			char *buf = new char[cl];
			if (!WideCharToMultiByte(65001, 0, wargv[i], wl + 1, buf, cl, nullptr, nullptr)) {
				std::cerr << "failed to convert argv[" << i << "] to utf-8: error " << GetLastError() << std::endl;
				abort();
			}
			argv[i] = buf;
		}
		main_cli(argc, argv.data());
	} else {
		main_gui();
	}
}
#else
int main(int argc, const char *argv[]) {
	init_console();
	if (argc > 1) {
		main_cli(argc, argv);
	} else {
		main_gui();
	}
}
#endif
