
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
		// TODO etc
	};

	struct saliency_progress {
		std::chrono::milliseconds elapsed_time{0};
		int completed_levels = 0;
		int completed_vertices = 0;
		int total_vertices = 0;
		bool should_cancel = false;
	};

	std::future<bool> compute_saliency_async(Model &model, const saliency_params &params, saliency_progress &progress);

}

#endif
