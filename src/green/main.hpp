
#pragma once

#ifndef GREEN_MAIN_HPP
#define GREEN_MAIN_HPP

#include "entity.hpp"
#include "saliency.hpp"

namespace green {

	entity_selection & ui_selection();

	saliency_user_params & ui_saliency_user_params();

	const saliency_progress & ui_saliency_progress();

}

namespace ImGui {

	void draw_saliency_progress(const green::saliency_progress &progress);

	bool edit_saliency_params(green::saliency_user_params &uparams);

	void draw_saliency_params(const green::saliency_user_params &uparams);

}

#endif
