
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
		float weight = 1.f;
		float power = 1.f;
		bool use_saliency = true;
		bool show_progress = true;

		std::string str() const {
			char buf[128];
			char *end = buf + sizeof(buf);
			char *p = buf;
			p += snprintf(p, end - p, "v=%d", targetverts);
			if (use_saliency) p += snprintf(p, end - p, ",b=%d,w=%.2f,p=%.2f", nbins, weight, power);
			return {buf};
		}

		explicit operator std::string() const {
			return str();
		}

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
