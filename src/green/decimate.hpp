
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#pragma once

#ifndef GREEN_DECIMATE_HPP
#define GREEN_DECIMATE_HPP

#include <cstdio>
#include <cstring>
#include <chrono>
#include <future>
#include <memory>
#include <functional>
#include <string>
#include <algorithm>

#include <imgui.h>

#include "meshutils.hpp"

namespace green {

	struct decimate_user_params {
		int targetverts = 10000;
		int nbins = 5;
		float weight_falloff = 0.1f;

		void sanitize();
	};

	struct decimate_mesh_params {
		TriMesh *mesh = nullptr;
		OpenMesh::VPropHandleT<float> prop_saliency;
	};

	enum class decimation_state {
		idle, bins, run, done, cancelled
	};

	struct decimate_progress {
		std::chrono::milliseconds elapsed_time{0};
		int target_collapses = 0;
		int completed_collapses = 0;
		bool should_cancel = false;
		decimation_state state = decimation_state::idle;
	};

	bool decimate(const decimate_mesh_params &mparams, const decimate_user_params &uparams, decimate_progress &progress);

}

#endif
