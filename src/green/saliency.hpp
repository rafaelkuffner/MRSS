
/*
 * Copyright 2020
 * Computational Media Innovation Centre
 * Victoria University of Wellington
 *
 */

#pragma once

#ifndef GREEN_SALIENCY_HPP
#define GREEN_SALIENCY_HPP

#include <chrono>
#include <future>

namespace green {

	class Model;

	struct saliency_params {
		int levels = 5;
		float area = 0.02f;
		float curv_weight = 0.f;
		bool normalmap_filter = false;
		// command line progress output
		bool show_progress = true;
		// TODO etc
	};

	struct saliency_progress {
		std::chrono::milliseconds elapsed_time{0};
		int completed_levels = 0;
		int completed_vertices = 0;
		int total_vertices = 0;
		bool should_cancel = false;
	};

	bool compute_saliency(Model &model, const saliency_params &params, saliency_progress &progress);

	std::future<bool> compute_saliency_async(Model &model, const saliency_params &params, saliency_progress &progress);

}

#endif
