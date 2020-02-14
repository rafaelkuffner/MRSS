
#include "saliency.hpp"
#include "model.hpp"

#include <thread>
#include <chrono>

namespace green {

	bool compute_saliency(Model &model, const saliency_params &params, saliency_progress &progress);

	std::future<bool> compute_saliency_async(Model &model, const saliency_params &params, saliency_progress &progress) {
		return std::async([&]() { return compute_saliency(model, params, progress); });
	}

	bool compute_saliency(Model &model, const saliency_params &params, saliency_progress &progress) {
		const auto time_start = std::chrono::steady_clock::now();
		auto &mesh = model.trimesh();
		progress.total_vertices = int(mesh.n_vertices());
		for (int level = 0; level < params.levels; level++) {
			for (int i = 0; i < progress.total_vertices; i += 10) {
				if (progress.should_cancel) return false;
				const auto now = std::chrono::steady_clock::now();
				progress.elapsed_time = std::chrono::duration_cast<decltype(progress.elapsed_time)>(now - time_start);
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				progress.completed_vertices = i + 1;
			}
			progress.completed_levels = level + 1;
		}
		return true;
	}

}
