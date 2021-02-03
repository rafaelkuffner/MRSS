#pragma once

#ifndef GREEN_MODEL_ENTITY_HPP
#define GREEN_MODEL_ENTITY_HPP

#include "model.hpp"

namespace green {

	class ModelEntity : public Entity {
	private:
		// unique: adding/removing properties, vertices etc
		// shared: reading/writing contents of existing properties
		std::shared_mutex m_modelmtx;
		std::unique_ptr<Model> m_model;

		std::filesystem::path m_fpath_load;
		std::future<std::unique_ptr<Model>> m_pending_load;
		decimate_user_params m_dec_uparams;
		decimate_progress m_dec_progress;
		bool m_decimated = false;

		std::filesystem::path m_fpath_save;
		bool m_save_binary = true;
		bool m_save_original_vids = false;
		bool m_save_ok = false;
		std::future<bool> m_pending_save;

		int m_saliency_index = 0;
		int m_saliency_export_index = 0;
		int m_saliency_baseline_index = 0;
		bool m_saliency_vbo_dirty = false;

		model_saliency_errors m_saliency_errors;
		float m_saliency_error_scale = 1;

		float m_scale = 1;
		glm::vec3 m_translation{0};
		// TODO store as quat, only use euler for editing?
		glm::vec3 m_rotation_euler_yxz{0};

		struct basis_vector {
			glm::vec3 v;
			const char *name;
		};

		basis_vector m_basis_vectors[6]{
			{{+1, 0, 0}, "+X"},
			{{-1, 0, 0}, "-X"},
			{{0, +1, 0}, "+Y"},
			{{0, -1, 0}, "-Y"},
			{{0, 0, +1}, "+Z"},
			{{0, 0, -1}, "-Z"}
		};

		int m_basis_right = 0, m_basis_up = 2, m_basis_back = 4;

		int m_vert_point_size = 3;
		model_color_mode m_disp_color_mode = model_color_mode::saliency;
		model_color_mode m_exp_color_mode = model_color_mode::saliency;
		bool m_color_faces = true;
		bool m_color_verts = false;
		bool m_show_faces = true;
		bool m_show_edges = false;
		bool m_show_verts = false;
		bool m_cull_faces = true;
		bool m_cull_edges = true;
		bool m_shade_flat = false;

		bool m_try_export = false;
		bool m_try_kill = false;
		bool m_dead = false;

		void update_vbo();
		void draw_window_models(bool selected);
		void draw_window_selection();
		void draw_window_export();
		void draw_window_decimation(bool selected);
		bool draw_select_header(bool selected);
		void spawn_locked_notification() const;
		std::string make_name_tooltip() const;
		void invalidate_saliency_vbo();

	public:
		ModelEntity();

		const Model * model() const {
			return m_model.get();
		}

		Model * model() {
			return m_model.get();
		}

		void load(const std::filesystem::path &fpath);

		void save(const std::filesystem::path &fpath);

		void load(const std::filesystem::path &fpath, Model m);

		void load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress);

		void try_export() {
			m_try_export = true;
		}

		bool show_entity_outline() const {
			return m_show_faces;
		}

		bool decimated() const {
			return m_decimated;
		}

		const model_saliency_data * selected_saliency() const {
			if (!m_model) return nullptr;
			if (m_saliency_index >= m_model->saliency().size()) return nullptr;
			return &m_model->saliency()[m_saliency_index];
		}

		virtual std::string name() const override {
			std::string s = m_fpath_load.filename().u8string();
			if (m_decimated) {
				s += " [dec ";
				s += m_dec_uparams.str();
				s += "]";
			}
			return s;
		}

		virtual void move_by(const glm::vec3 &d) override;

		virtual void try_kill() override {
			m_try_kill = true;
		}

		virtual bool dead() const override {
			return m_dead;
		}

		virtual glm::mat4 transform() const override;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) override;

		virtual std::future<saliency_result> compute_saliency_async(const saliency_user_params &uparams, saliency_progress &progress) override;

		std::unique_ptr<ModelEntity> decimate_async(const decimate_user_params &uparams);

		virtual ~ModelEntity();
	};

}

#endif
