﻿
#include "model_entity.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/euler_angles.hpp>

#include <imgui.h>

#include "imguiex.hpp"
#include "main.hpp"

namespace {

	auto & mesh_io_mutex() {
		static std::mutex m;
		return m;
	}

	struct saliency_clipboard_data {
		green::model_saliency_data sd;
		std::vector<float> data;
	};

	auto & saliency_clipboard() {
		static saliency_clipboard_data d;
		return d;
	}

}

namespace green {

	ModelEntity::ModelEntity() {

	}

	void ModelEntity::load(const std::filesystem::path &fpath) {
		// lock in case of premature closure
		// mainly because the pending_load future's dtor will block otherwise
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async([=, lock=std::move(lock)]() {
			// apparently openmesh load is not threadsafe?
			// TODO check assimp too
			std::lock_guard iolock(mesh_io_mutex());
			return std::make_unique<Model>(fpath);
		});
	}

	void ModelEntity::save(const std::filesystem::path &fpath) {
		if (!m_model) return;
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		const auto &saliency = m_model->saliency();
		const auto sdcolor = m_saliency_export_index < saliency.size() ? saliency[m_saliency_export_index] : model_saliency_data{};
		const auto sdbase = m_saliency_baseline_index < saliency.size() ? saliency[m_saliency_baseline_index] : model_saliency_data{};
		m_fpath_save = fpath;
		m_pending_save = std::async([=, lock=std::move(lock)]() {
			// TODO is this also not threadsafe? idk
			// ... can it run parallel with load?
			std::lock_guard iolock(mesh_io_mutex());
			model_save_params sparams;
			sparams.prop_saliency = sdcolor.prop_saliency;
			sparams.prop_saliency_baseline = sdbase.prop_saliency;
			sparams.color_mode = m_exp_color_mode;
			sparams.original_vids = m_save_original_vids;
			sparams.binary = m_save_binary;
			m_model->save(fpath, sparams);
			return true;
		});
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m) {
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async(std::launch::deferred, [m=std::move(m)]() mutable {
			return std::make_unique<Model>(std::move(m));
		});
		// non-blocking wait will never succeed otherwise
		m_pending_load.wait();
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress) {
		load(fpath, std::move(m));
		m_decimated = true;
		m_dec_uparams = dec_uparams;
		m_dec_progress = dec_progress;
	}

	void ModelEntity::move_by(const glm::vec3 &d) {
		m_translation += d;
		if (d != glm::vec3(0)) invalidate_scene();
	}

	glm::mat4 ModelEntity::transform() const {
		if (!m_model) return glm::mat4(1);
		glm::mat4 transform(1);
		transform = glm::translate(transform, m_translation);
		transform *= glm::eulerAngleYXZ(m_rotation_euler_yxz.y, m_rotation_euler_yxz.x, m_rotation_euler_yxz.z);
		transform = glm::scale(transform, glm::vec3(m_scale));
		glm::mat3 basis(m_basis_vectors[m_basis_right].v, m_basis_vectors[m_basis_up].v, m_basis_vectors[m_basis_back].v);
		return transform * glm::mat4(transpose(basis));
	}

	void ModelEntity::update_vbo() {
		if (!m_model) return;
		if (!m_saliency_vbo_dirty) return;
		std::shared_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) return;
		auto &saliency_outputs = m_model->saliency();
		model_color_params cparams;
		if (m_saliency_index < saliency_outputs.size()) cparams.prop_saliency = saliency_outputs[m_saliency_index].prop_saliency;
		if (m_saliency_baseline_index < saliency_outputs.size()) cparams.prop_saliency_baseline = saliency_outputs[m_saliency_baseline_index].prop_saliency;;
		cparams.color_mode = m_disp_color_mode;
		cparams.error_scale = m_saliency_error_scale;
		if (!m_model->update_color(cparams, &m_saliency_errors)) m_disp_color_mode = model_color_mode::none;
		m_saliency_vbo_dirty = false;
		invalidate_scene();
	}

	void ModelEntity::draw_window_models(bool selected) {
		using namespace ImGui;
		if (Begin("Models")) {
			PushID(this);
			bool dirty = false;
			draw_select_header(selected);
			if (m_pending_load.valid()) {
				Text("Loading...");
			} else if (m_model) {
				// ok
			} else {
				Text("Failed to load model");
			}
			dirty |= Checkbox("Faces", &m_show_faces);
			SameLine();
			dirty |= Checkbox("Edges", &m_show_edges);
			SameLine();
			dirty |= Checkbox("Verts", &m_show_verts);
			SameLine();
			SetCursorPosX(GetCursorPosX() + std::max(0.f, GetContentRegionAvail().x - 20));
			PushStyleColor(ImGuiCol_Button, ImVec4{0.6f, 0.3f, 0.3f, 1});
			if (Button("X", {-1, 0}) || m_try_kill) {
				OpenPopup("##close");
				m_try_kill = false;
			}
			if (IsItemHovered()) SetTooltip("Close model");
			PopStyleColor();
			if (BeginPopupModal("##close", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoDecoration)) {
				Text("Close model \"%s\" ?", name().c_str());
				SetHoveredTooltip("%s", make_name_tooltip().c_str());
				if (Button("Close")) {
					std::unique_lock lock(m_modelmtx, std::defer_lock);
					if (lock.try_lock()) {
						m_dead = true;
						invalidate_scene();
					} else {
						spawn_locked_notification();
					}
					CloseCurrentPopup();
				}
				SameLine();
				if (Button("Cancel")) {
					CloseCurrentPopup();
				}
				EndPopup();
			}
			Separator();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
			if (dirty) invalidate_scene();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_selection() {
		using namespace ImGui;
		if (Begin("Selection")) {
			PushID(this);
			PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
			Selectable(name().c_str(), true, 0, {0, GetTextLineHeightWithSpacing()});
			PopStyleColor();
			SetHoveredTooltip("%s", make_name_tooltip().c_str());
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);

			Text("%zd vertices, %zd triangles", m_model->trimesh().n_vertices(), m_model->trimesh().n_faces());
			if (m_decimated) Text(
				"Decimated %d vertices (%.1f%%) in %.3fs",
				m_dec_progress.completed_collapses,
				100.f * m_dec_progress.completed_collapses / float(m_dec_uparams.targetverts + m_dec_progress.completed_collapses),
				m_dec_progress.elapsed_time / std::chrono::duration<double>(1.0)
			);

			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			if (m_pending_save.valid()) {
				if (m_pending_save.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
					try {
						m_save_ok = m_pending_save.get();
					} catch (std::exception &e) {
						std::cerr << "failed to save model: " << e.what() << std::endl;
						m_save_ok = false;
					} catch (...) {
						std::cerr << "failed to save model" << std::endl;
						m_save_ok = false;
					}
				} else {
					TextDisabled("Exporting %s...", m_fpath_save.filename().u8string().c_str());
				}
			} else if (m_save_ok) {
				TextDisabled("Export %s succeeded", m_fpath_save.filename().u8string().c_str());
			} else if (!m_fpath_save.empty()) {
				TextColored(badcol, "Export %s failed", m_fpath_save.filename().u8string().c_str());
			} else {
				TextDisabled("Not exported");
			}
			if (IsItemHovered() && !m_fpath_save.empty()) SetTooltip("%s", m_fpath_save.u8string().c_str());

			if (CollapsingHeader("Transform")) {
				auto pick_basis = [&](const char *label, int *basis) {
					SetNextItemWidth(GetTextLineHeight() * 3.5f);
					Combo(
						label, basis,
						[](void *data, int item, const char **out_text) {
							auto vs = reinterpret_cast<basis_vector *>(data);
							*out_text = vs[item].name;
							return true;
						}, m_basis_vectors, 6
					);
				};
				pick_basis("Right", &m_basis_right);
				SameLine();
				pick_basis("Up", &m_basis_up);
				SameLine();
				pick_basis("Back", &m_basis_back);
				if (Button("Reset##scale")) m_scale = m_model->unit_bound_scale() * 4;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderFloat("Scale", &m_scale, 0, 1000, "%.4f", 8);
				if (Button("Reset##translation")) m_translation = glm::vec3{0};
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				InputFloat3("Translation", value_ptr(m_translation));
				if (Button("Reset##roty")) m_rotation_euler_yxz.y = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Yaw", &m_rotation_euler_yxz.y, -180, 180);
				if (Button("Reset##rotx")) m_rotation_euler_yxz.x = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Pitch", &m_rotation_euler_yxz.x, -180, 180);
				if (Button("Reset##rotz")) m_rotation_euler_yxz.z = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Roll", &m_rotation_euler_yxz.z, -180, 180);
			}

			if (CollapsingHeader("Rendering")) {
				bool r = false;
				if (RadioButton("Smooth Shading", !m_shade_flat)) {
					m_shade_flat = false;
					r = true;
				}
				SameLine();
				if (RadioButton("Flat Shading", m_shade_flat)) {
					m_shade_flat = true;
					r = true;
				}
				r |= Checkbox("Color Faces", &m_color_faces);
				SetHoveredTooltip("Apply color mode to faces");
				SameLine();
				r |= Checkbox("Color Verts", &m_color_verts);
				SetHoveredTooltip("Apply color mode to vertices");
				r |= Checkbox("Cull Faces", &m_cull_faces);
				SameLine();
				r |= Checkbox("Cull Edges", &m_cull_edges);
				r |= SliderInt("Point Size", &m_vert_point_size, 1, 5);
				m_vert_point_size = std::max(1, m_vert_point_size);
				if (r) invalidate_scene();
			}

			Separator();

			//SetNextItemWidth(GetContentRegionAvail().x);
			int cur_color_mode = int(m_disp_color_mode);
			if (Combo("Color Mode", &cur_color_mode, "None\0Vertex Color\0Saliency\0Saliency Comparison\0Curvature (DoN)\0")) {
				m_disp_color_mode = model_color_mode(cur_color_mode);
				invalidate_saliency_vbo();
			}

			Separator();

			std::vector<model_saliency_data> empty_saliency_outputs;
			auto &saliency_outputs = m_model ? m_model->saliency() : empty_saliency_outputs;

			if (Combo(
				"Saliency", &m_saliency_index,
				[](void *data, int item, const char **out_text) {
					auto ds = reinterpret_cast<model_saliency_data *>(data);
					const auto &so = ds[item];
					// TODO better?
					static std::string str;
					str = std::string(so);
					*out_text = str.c_str();
					return true;
				}, 
				saliency_outputs.data(), saliency_outputs.size()
					)) {
				invalidate_saliency_vbo();
			}

			if (m_disp_color_mode == model_color_mode::saliency_comparison && m_saliency_baseline_index < saliency_outputs.size()) {
				Text("Baseline: %s", saliency_outputs[m_saliency_baseline_index].str().c_str());
				if (m_saliency_index < saliency_outputs.size()) {
					auto &err = m_saliency_errors;
					Text("Errors: min=%.3f, max=%.3f, rmse=%.3f", err.min, err.max, err.rms);
					if (SliderFloat("Error Scale", &m_saliency_error_scale, 1, 100, "%.3f", 2)) {
						invalidate_saliency_vbo();
					}
				}
			}

			if (Button("Paste")) {
				auto &clip = saliency_clipboard();
				std::unique_lock lock(m_modelmtx, std::defer_lock);
				if (!lock.try_lock()) {
					spawn_locked_notification();
				} else if (clip.data.size() == m_model->trimesh().n_vertices()) {
					model_saliency_data sd = std::move(clip.sd);
					m_model->trimesh().add_property(sd.prop_saliency);
					m_model->trimesh().property(sd.prop_saliency).data_vector() = std::move(clip.data);
					clip.sd = {};
					clip.data.clear();
					saliency_outputs.push_back(std::move(sd));
					invalidate_saliency_vbo();
					// references/iterators into saliency outputs are invalidated
				} else {
					OpenPopup("Paste Error##pasteerror");
				}
			}

			if (m_saliency_index < saliency_outputs.size()) {
				auto &salout = saliency_outputs[m_saliency_index];
				SameLine();
				if (Button("Copy")) {
					std::shared_lock lock(m_modelmtx, std::defer_lock);
					if (!lock.try_lock()) {
						spawn_locked_notification();
					} else {
						auto &clip = saliency_clipboard();
						clip.sd = salout;
						clip.sd.prop_saliency.reset();
						clip.data = m_model->trimesh().property(salout.prop_saliency).data_vector();
					}
				}
				SameLine();
				if (Button("Remove")) OpenPopup("##remove");
				SameLine();
				if (Button("Rename")) OpenPopup("Rename Saliency##rename");
				SameLine();
				if (Button("Baseline")) m_saliency_baseline_index = m_saliency_index;
				Checkbox("Persistent", &salout.persistent);
				SetHoveredTooltip("Persistent properties will be preserved when decimating\nand will be exported by default");
				Separator();
				if (CollapsingHeader("Saliency Parameters")) {
					if (salout.uparams_known) {
						draw_saliency_params(salout.uparams);
						if (Button("Reload")) ui_saliency_user_params() = salout.uparams;
					} else {
						TextDisabled("Parameters unknown");
					}
				}
				if (salout.filename.empty()) {
					if (CollapsingHeader("Saliency Progress")) {
						draw_saliency_progress(salout.progress);
					}
				}
				if (BeginPopupModal("Paste Error##pasteerror")) {
					auto &clip = saliency_clipboard();
					if (clip.data.size() == m_model->trimesh().n_vertices()) {
						CloseCurrentPopup();
					} else {
						Text("Can't paste saliency for %d vertices into model with %d vertices", clip.data.size(), m_model->trimesh().n_vertices());
						if (Button("OK")) CloseCurrentPopup();
					}
				}
				bool rename_open = true;
				SetNextWindowSize({300, 100}, ImGuiCond_Appearing);
				if (BeginPopupModal("Rename Saliency##rename", &rename_open)) {
					Text("New saliency property name");
					// TODO better?
					static char buf[1024]{};
					if (IsWindowAppearing()) {
						strncpy(buf, salout.dispname.c_str(), sizeof(buf));
						SetKeyboardFocusHere();
					}
					SetNextItemWidth(GetContentRegionAvail().x);
					const char *badchars = " \t@[]";
					if (InputText("", buf, sizeof(buf), ImGuiInputTextFlags_EnterReturnsTrue)) {
						if (std::string_view(buf).find_first_of(badchars) == std::string::npos) {
							salout.dispname = buf;
							CloseCurrentPopup();
						}
					}
					if (std::string_view(buf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of: %s", badchars);
					}
					EndPopup();
				}
				if (BeginPopupModal("##remove")) {
					Text("Remove saliency result?");
					if (Button("Remove")) {
						std::unique_lock lock(m_modelmtx, std::defer_lock);
						if (!lock.try_lock()) {
							spawn_locked_notification();
						} else {
							auto it = saliency_outputs.begin() + m_saliency_index;
							m_model->trimesh().remove_property(it->prop_saliency);
							saliency_outputs.erase(it);
							// references/iterators into saliency outputs are invalidated
							// (so this section should come last)
							if (m_saliency_index >= saliency_outputs.size()) {
								m_saliency_index = std::max(0, m_saliency_index - 1);
								invalidate_saliency_vbo();
							}
						}
						CloseCurrentPopup();
					}
					SameLine();
					if (Button("Cancel")) CloseCurrentPopup();
					EndPopup();
				}
			} else {
				SameLine();
				ButtonDisabled("Copy");
				SameLine();
				ButtonDisabled("Remove");
				SameLine();
				ButtonDisabled("Rename");
				SameLine();
				ButtonDisabled("Baseline");
			}
			Separator();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_export() {
		using namespace ImGui;
		if (!m_model) return;
		auto &saliency_outputs = m_model->saliency();
		if (m_try_export) {
			m_try_export = false;
			m_saliency_export_index = m_saliency_index;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				sd.should_export = sd.persistent;
			}
			OpenPopup("Export##export");
		}
		SetNextWindowSize({500, 400}, ImGuiCond_Appearing);
		bool export_window_open = true;
		if (BeginPopupModal("Export##export", &export_window_open)) {
			PushID(this);
			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
			Selectable(m_fpath_load.filename().u8string().c_str(), true, ImGuiSelectableFlags_DontClosePopups, {0, GetTextLineHeightWithSpacing()});
			PopStyleColor();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
			SetHoveredTooltip("%s", make_name_tooltip().c_str());
			// export path
			const auto &pathhint = m_fpath_save.empty() ? m_fpath_load : m_fpath_save;
			// TODO better?
			static char pathbuf[1024]{};
			if (IsWindowAppearing()) snprintf(pathbuf, sizeof(pathbuf), "%s", pathhint.u8string().c_str());
			auto &fpath_save_fut = ui_save_path(pathhint, Button("Browse"));
			if (fpath_save_fut.valid() && fpath_save_fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				auto s = fpath_save_fut.get().u8string();
				if (!s.empty()) snprintf(pathbuf, sizeof(pathbuf), "%s", s.c_str());
			}
			SameLine();
			InputText("Path", pathbuf, sizeof(pathbuf));
			// export color mode
			int cur_color_mode = int(m_exp_color_mode);
			if (Combo("Color Mode", &cur_color_mode, "None\0Vertex Color\0Saliency\0Saliency Comparison\0Curvature (DoN)\0")) {
				m_exp_color_mode = model_color_mode(cur_color_mode);
			}
			if (m_exp_color_mode == model_color_mode::saliency_comparison) {
				if (m_saliency_baseline_index < m_model->saliency().size()) {
					Text("Baseline: %s", m_model->saliency()[m_saliency_baseline_index].str().c_str());
				} else {
					TextColored(badcol, "Invalid baseline");
				}
			}
			Separator();
			Text("Saliency Properties");
			if (BeginChild("saliency", {0, -115}, false, ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
				for (int i = 0; i < saliency_outputs.size(); i++) {
					PushID(i);
					auto &sd = saliency_outputs[i];
					Checkbox("##shouldexport", &sd.should_export);
					SetHoveredTooltip("Export this saliency property?");
					SameLine();
					if (RadioButton("##colorize", i == m_saliency_export_index)) m_saliency_export_index = i;
					SetHoveredTooltip("Colorize this saliency property?");
					const char *badchars = " \t@[]";
					char propbuf[256]{};
					strncpy(propbuf, sd.expname.c_str(), sizeof(propbuf));
					SameLine();
					SetNextItemWidth(150);
					if (InputText("##expname", propbuf, sizeof(propbuf))) {
						if (std::string_view(propbuf).find_first_of(badchars) == std::string::npos) {
							sd.expname = propbuf;
						}
					}
					SetHoveredTooltip("Export property name.\nDefault when empty is the current display name.");
					SameLine();
					if (sd.expname.empty()) {
						auto p0 = GetCursorPos();
						SetCursorPosX(p0.x - 150);
						TextDisabled("%s", sd.dispname.c_str());
						SameLine();
						SetCursorPosX(p0.x);
					}
					if (Button("<")) sd.expname = sd.dispname;
					SetHoveredTooltip("Copy display name to export name");
					SameLine();
					Text("%s", sd.str().c_str());
					SetHoveredTooltip(sd.uparams_known ? sd.uparams.str().c_str() : "Parameters unknown");
					if (std::string_view(propbuf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of \"%s\"", badchars);
					}
					PopID();
				}
			}
			EndChild();
			Separator();
			// other options
			Checkbox("Binary", &m_save_binary);
			SetHoveredTooltip("Save file as binary if supported");
			SameLine();
			Checkbox("Original Vertex IDs", &m_save_original_vids);
			SetHoveredTooltip("Export extra property with original vertex ids from before decimation");
			TextDisabled("Supported formats are currently PLY and OBJ.\nColorization and property export are only supported for PLY.");
			// validate path and maybe export
			auto fpath = std::filesystem::u8path(pathbuf);
			auto stat = std::filesystem::status(fpath);
			bool cansave = !fpath.empty();
			const char *badchars = "\\/:*?\"<>|";
			if (fpath.is_relative()) {
				TextColored(badcol, "Path must be absolute");
				cansave = false;
			} else if (std::filesystem::is_directory(stat)) {
				TextColored(badcol, "Path is a directory");
				cansave = false;
			} else if (!std::filesystem::is_directory(fpath.parent_path())) {
				TextColored(badcol, "Directory does not exist");
				cansave = false;
			} else if (fpath.filename().u8string().find_first_of(badchars) != std::string::npos) {
				TextColored(badcol, "File name must not contain any of \"%s\"", badchars);
				cansave = false;
			} else if (std::filesystem::exists(stat)) {
				TextColored(badcol, "Path exists! Save will overwrite");
				//TextColored(badcol, u8"保存先が存在します！保存したら上書きします");
			}
			if (Button("Save") && cansave && !m_pending_save.valid()) {
				save(std::filesystem::absolute(fpath));
				CloseCurrentPopup();
			}
			if (IsItemHovered()) {
				if (m_pending_save.valid()) {
					SetTooltip("Export already in progress");
				} else if (!cansave) {
					SetTooltip("Can't save to this path");
				}
			}
			SameLine();
			if (Button("Cancel")) CloseCurrentPopup();
			PopID();
			EndPopup();
		}
	}

	void ModelEntity::draw_window_decimation(bool selected) {
		using namespace ImGui;
		if (m_decimated && ui_decimation_window_open() && (m_dec_progress.state < decimation_state::done)) {
			if (Begin("Decimation")) {
				PushID(this);
				draw_select_header(selected);
				draw_decimate_progress(m_dec_progress);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				PopID();
			}
		}
	}

	bool ModelEntity::draw_select_header(bool selected) {
		using namespace ImGui;
		bool r = false;
		PushStyleColor(ImGuiCol_Header, selected ? ImVec4{0.7f, 0.4f, 0.1f, 1} : GetStyle().Colors[ImGuiCol_Button]);
		// hack - selectable is always 'selected' in order to show highlight, it just changes colour
		if (Selectable(name().c_str(), true, 0, {0, GetTextLineHeightWithSpacing()})) {
			auto &sel = ui_selection();
			sel.select_entity = id();
			sel.select_vertex = -1;
			r = true;
		}
		PopStyleColor();
		SetHoveredTooltip("%s", make_name_tooltip().c_str());
		SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
		return r;
	}

	void ModelEntity::spawn_locked_notification() const {
		using namespace ImGui;
		ui_spawn([nametooltip=make_name_tooltip(), name=name()](bool *p_open) {
			if (*p_open) OpenPopup("Model in use");
			if (BeginPopupModal("Model in use", p_open, ImGuiWindowFlags_AlwaysAutoResize)) {
				Selectable(name.c_str(), true, ImGuiSelectableFlags_DontClosePopups, {0, GetTextLineHeightWithSpacing()});
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				SetHoveredTooltip("%s", nametooltip.c_str());
				Text("The operation cannot be performed because\nthe model is currently in use.");
				if (Button("  OK  ")) {
					*p_open = false;
					CloseCurrentPopup();
				}
				EndPopup();
			}
		});
	}

	std::string ModelEntity::make_name_tooltip() const {
		auto s = m_fpath_load.u8string();
		if (m_decimated) {
			s += "\n";
			s += m_dec_uparams.str();
		}
		return s;
	}

	void ModelEntity::invalidate_saliency_vbo() {
		m_saliency_vbo_dirty = true;
		invalidate_scene();
	}

	void ModelEntity::draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) {
		const auto xform0 = transform();
		if (m_pending_load.valid() && m_pending_load.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
			try {
				// is this the best place to do this?
				m_model = m_pending_load.get();
				// gl stuff has to happen on main thread
				m_model->update_vao();
				m_model->update_vbos();
				m_scale = m_model->unit_bound_scale() * 4;
				invalidate_saliency_vbo();
			} catch (std::exception &e) {
				std::cerr << "failed to load model: " << e.what() << std::endl;
			} catch (...) {
				std::cerr << "failed to load model" << std::endl;
			}
		}
		auto &sel = ui_selection();
		const bool selected = sel.select_entity == id();
		if (!dead() && selected) ui_select_model(this);
		if (!dead() && sel.hover_entity == id()) ui_hover_model(this);
		draw_window_models(selected);
		draw_window_decimation(selected);
		if (m_model) {
			if (selected) draw_window_selection();
			if (selected) draw_window_export();
			if (draw_scene) {
				update_vbo();
				// determine color map to apply in shader
				auto &saliency_outputs = m_model->saliency();
				int color_map = 0;
				const bool sal_valid = m_saliency_index < saliency_outputs.size();
				if (m_disp_color_mode == model_color_mode::saliency && sal_valid) color_map = 3;
				if (m_disp_color_mode == model_color_mode::saliency_comparison && sal_valid && m_saliency_baseline_index < saliency_outputs.size()) color_map = 4;
				if (m_disp_color_mode == model_color_mode::vcolor && m_model->prop_vcolor_original().is_valid()) color_map = 1;
				if (m_disp_color_mode == model_color_mode::doncurv) color_map = 3;
				// prepare to draw
				model_draw_params params;
				params.sel = sel;
				params.entity_id = id();
				auto set_cull_faces = [](bool b) {
					if (b) {
						glEnable(GL_CULL_FACE);
						glCullFace(GL_BACK);
					} else {
						glDisable(GL_CULL_FACE);
					}
				};
				// faces
				params.shade_flat = m_shade_flat;
				params.color = {0.6f, 0.6f, 0.5f, 1};
				params.vert_color_map = m_color_faces ? color_map : 0;
				set_cull_faces(m_cull_faces);
				glColorMaski(1, GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_faces) m_model->draw(view * transform(), proj, zfar, params, GL_FILL);
				// edges
				params.shade_flat = false;
				params.shading = 0;
				params.color = {0.03f, 0.03f, 0.03f, 1};
				params.vert_color_map = 0;
				set_cull_faces(m_cull_edges);
				glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_edges) m_model->draw(view * transform(), proj, zfar, params, GL_LINE);
				// verts
				params.color = {0.5f, 0, 0, 1};
				params.vert_color_map = m_color_verts ? color_map : 0;
				params.sel.hover_entity = -1;
				params.sel.select_entity = -1;
				//params.show_samples = sal_valid && (m_color_mode == color_mode::saliency || m_color_mode == color_mode::saliency_comparison);
				glDisable(GL_CULL_FACE);
				glPointSize(m_vert_point_size);
				glColorMaski(1, GL_TRUE, GL_TRUE, GL_FALSE, GL_FALSE);
				if (m_show_verts) m_model->draw(view * transform(), proj, zfar, params, GL_POINT);
			}
		}
		const auto xform1 = transform();
		if (xform0 != xform1) invalidate_scene();
	}

	std::future<saliency_result> ModelEntity::compute_saliency_async(const saliency_user_params &uparams0, saliency_progress &progress) {
		if (!m_model) return {};
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return {};
		}
		saliency_user_params uparams = uparams0;
		if (uparams.auto_contrast) uparams.normal_power = m_model->auto_contrast();
		saliency_mesh_params mparams;
		mparams.mesh = &m_model->trimesh();
		mparams.prop_vertex_area = m_model->prop_vertex_area();
		mparams.prop_doncurv_raw = m_model->prop_doncurv_raw();
		mparams.prop_edge_length = m_model->prop_edge_length();
		// create properties
		mparams.prop_saliency_levels.resize(uparams.levels);
		mparams.mesh->add_property(mparams.prop_curvature);
		mparams.mesh->add_property(mparams.prop_saliency);
		mparams.mesh->add_property(mparams.prop_sampled);
		for (int i = 0; i < uparams.levels; i++) {
			mparams.mesh->add_property(mparams.prop_saliency_levels[i]);
		}
		// cleanup will be run when the result is received from the future
		// can use shared lock for actual computation because only property contents are being used
		// HACK make_shared is an ugly workaround for std::function needing to be copyable
		lock.unlock();
		mparams.cleanup = [=, pprogress=&progress, slock=std::make_shared<std::shared_lock<std::shared_mutex>>(m_modelmtx)](bool r) mutable noexcept {
			// now need to upgrade to unique lock
			slock->unlock();
			std::unique_lock lock(m_modelmtx);
			if (!m_model) return;
			auto &saliency_outputs = m_model->saliency();
			// destroy temp properties
			mparams.mesh->remove_property(mparams.prop_curvature);
			for (int i = 0; i < uparams.levels; i++) {
				mparams.mesh->remove_property(mparams.prop_saliency_levels[i]);
			}
			if (r) {
				// if preview, remove any previous preview results
				saliency_outputs.erase(std::remove_if(saliency_outputs.begin(), saliency_outputs.end(), [](const auto &sd) { return sd.uparams.preview; }), saliency_outputs.end());
				// save user params, progress output and actual saliency mesh property
				model_saliency_data sd;
				sd.uparams = uparams;
				sd.progress = *pprogress;
				sd.prop_saliency = mparams.prop_saliency;
				sd.prop_sampled = mparams.prop_sampled;
				sd.uparams_known = true;
				saliency_outputs.push_back(sd);
				// give focus to this result
				m_saliency_index = saliency_outputs.size() - 1;
				invalidate_saliency_vbo();
			} else {
				// cancelled, destroy saliency property too
				mparams.mesh->remove_property(mparams.prop_saliency);
				mparams.mesh->remove_property(mparams.prop_sampled);
			}
		};
		return green::compute_saliency_async(mparams, uparams, progress);
	}

	std::unique_ptr<ModelEntity> ModelEntity::decimate_async(const decimate_user_params &uparams) {
		// this function isn't const because it needs to lock the mutex - mutable?
		std::shared_lock lock1(m_modelmtx, std::defer_lock);
		if (!lock1.try_lock()) {
			spawn_locked_notification();
			return {};
		}
		auto e = std::make_unique<ModelEntity>();
		std::unique_lock lock2(e->m_modelmtx);
		e->m_fpath_load = m_fpath_load;
		e->m_decimated = true;
		e->m_dec_uparams = uparams;
		e->m_basis_right = m_basis_right;
		e->m_basis_up = m_basis_up;
		e->m_basis_back = m_basis_back;
		e->m_rotation_euler_yxz = m_rotation_euler_yxz;
		e->m_cull_faces = m_cull_faces;
		e->m_cull_edges = m_cull_edges;
		e->m_pending_load = std::async([=, e=e.get(), lock1=std::move(lock1), lock2=std::move(lock2)]() mutable {
			// prepare mesh copy to decimate
			// this is slow enough on big models to be async
			auto &saliency_outputs = m_model->saliency();
			std::vector<model_saliency_data> sdv;
			saliency_prop_t srcprop;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				if (i == m_saliency_index || sd.persistent) {
					srcprop = sd.prop_saliency;
					sdv.push_back(sd);
				} else if (sd.persistent) {
					sdv.push_back(sd);
				}
			}
			auto m = m_model->prepare_decimate(srcprop, sdv);
			// unlock source model
			lock1.unlock();
			// now decimate
			if (!m.decimate(uparams, e->m_dec_progress)) throw std::runtime_error("decimation was cancelled");
			return std::make_unique<Model>(std::move(m));
		});
		return e;
	}

	ModelEntity::~ModelEntity() {

	}

}