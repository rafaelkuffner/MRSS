
#pragma once

#ifndef GREEN_MAIN_HPP
#define GREEN_MAIN_HPP

#include <future>
#include <filesystem>
#include <functional>

#include "entity.hpp"
#include "model_entity.hpp"
#include "saliency.hpp"
#include "decimate.hpp"

namespace green {

	struct user_options {
		// 0: all available
		int threads = 0;
	};

	void spawn_entity(std::unique_ptr<Entity>);

	entity_selection & ui_selection();

	void ui_select_model(ModelEntity *);

	ModelEntity * ui_select_model();

	void ui_hover_model(ModelEntity *);

	ModelEntity * ui_hover_model();

	std::future<std::filesystem::path> & ui_save_path(const std::filesystem::path &hint, bool prompt);

	// spawn persistent ui component. persists until *p_open is set to false.
	void ui_spawn(std::function<void(bool *p_open)>);

	bool ui_saliency_window_open();

	bool ui_decimation_window_open();

	void invalidate_scene();

}

namespace ImGui {

	void draw_saliency_progress(green::saliency_progress &progress);

	bool edit_saliency_params(green::saliency_user_params &uparams);

	void draw_saliency_params(const green::saliency_user_params &uparams);

	void draw_decimate_progress(green::decimate_progress &progress);

	bool edit_decimate_params(green::decimate_user_params &uparams);

}

#endif
