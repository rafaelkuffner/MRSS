
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
		int targettris = 10000;
		int nbins = 5;
		float bin_weight = 1.f;
		float sal_weight = 20.f;
		float bin_power = 1.f;
		float sal_power = 2.f;
		float max_aspect = 5.f;
		bool use_bins = false;
		bool use_tris = false;
		bool use_saliency = true;
		bool limit_aspect = true;
		bool preserve_seams = true;
		bool prevent_folds = true;
		bool show_progress = true;

		std::string str() const {
			char buf[128];
			char *end = buf + sizeof(buf);
			char *p = buf;
			if (use_tris) {
				p += snprintf(p, end - p, "t=%d", targettris);
			} else {
				p += snprintf(p, end - p, "v=%d", targetverts);
			}
			if (use_saliency) {
				if (use_bins) {
					p += snprintf(p, end - p, ",b=%d,w=%.2f,p=%.2f", nbins, bin_weight, bin_power);
				} else {
					p += snprintf(p, end - p, ",w=%.2f,p=%.2f", sal_weight, sal_power);
				}
			}
			// TODO update with aspect, seams, folds etc
			return {buf};
		}

		explicit operator std::string() const {
			return str();
		}

		void sanitize();
	};

	struct decimate_mesh_params {
		PolyMesh *mesh = nullptr;
		OpenMesh::VPropHandleT<float> prop_saliency;
		// output
		OpenMesh::VPropHandleT<float> prop_dec_error;
	};

	enum class decimation_state {
		idle, init, run, done, cancelled
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
