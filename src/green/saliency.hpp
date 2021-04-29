
/*
 * Copyright 2020
 * Computational Media Innovation Centre
 * Victoria University of Wellington
 *
 */

#pragma once

#ifndef GREEN_SALIENCY_HPP
#define GREEN_SALIENCY_HPP

#include <cstdio>
#include <cstring>
#include <chrono>
#include <future>
#include <vector>
#include <functional>
#include <string>
#include <string_view>

#include <imgui.h>

#include "meshutils.hpp"
#include "Histogram.h"

namespace green {

	using saliency_prop_t = OpenMesh::VPropHandleT<float>;

	class saliency_result {
	private:
		std::function<void(bool)> m_cleanup;
		bool m_success = false;

		void cleanup() noexcept {
			if (m_cleanup) m_cleanup(m_success);
		}

	public:
		saliency_result() = default;

		saliency_result(std::function<void(bool) /*noexcept*/> cleanup_, bool success_)
			: m_cleanup(std::move(cleanup_))
			, m_success(success_)
		{}

		saliency_result(const saliency_result &) = delete;
		saliency_result & operator=(const saliency_result &) = delete;

		saliency_result(saliency_result &&other) noexcept
			: m_cleanup(std::move(other.m_cleanup))
			, m_success(other.m_success)
		{
			other.m_cleanup = {};
		}

		saliency_result & operator=(saliency_result &&other) noexcept {
			cleanup();
			m_cleanup = std::move(other.m_cleanup);
			m_success = other.m_success;
			other.m_cleanup = {};
			return *this;
		}

		explicit operator bool() const {
			return m_success;
		}

		bool operator!() const {
			return !bool(*this);
		}

		~saliency_result() {
			cleanup();
		}

	};

	struct saliency_user_params {
		int levels = 5;
		float area = 0.02f;
		float curv_weight = 0;
		float normal_power = 1;
		// used with manual subsampling
		float subsampling_rate = 1;
		// used with auto subsampling
		float samples_per_neighborhood = 100;
		// used with normalmap filter; relative to sqrt(real_surface_area)
		float noise_height = 0.002f;
		// filter out normalmappable detail
		bool normalmap_filter = false;
		// enable manual subsampling (takes precedence)
		bool subsample_manual = false;
		// enable auto subsampling
		bool subsample_auto = true;
		// use automatic contrast (normal_power ignored if true)
		bool auto_contrast = true;
		// command line progress output
		bool show_progress = true;
		// interactive previewing
		bool preview = false;

		std::string str(bool verbose = false) const;

		bool parse(std::string_view);

		explicit operator std::string() const {
			return str();
		}

		bool operator==(const saliency_user_params &other) const {
			// TODO quick hack assuming pod
			return std::memcmp(this, &other, sizeof(*this)) == 0;
		}

		bool operator!=(const saliency_user_params &other) const {
			return !(*this == other);
		}

		// clamp parameters to valid ranges
		void sanitize();

	};

	struct saliency_mesh_params {
		// TODO const
		PolyMesh *mesh = nullptr;
		// input
		OpenMesh::VPropHandleT<float> prop_vertex_area;
		// input
		OpenMesh::VPropHandleT<float> prop_curv_raw;
		// input
		OpenMesh::EPropHandleT<float> prop_edge_length;
		// input
		float curv_min = 0;
		float curv_max = 1;
		// output: saliency
		saliency_prop_t prop_saliency;
		// output: sample location flags (bool, uchar to avoid vector<bool>)
		OpenMesh::VPropHandleT<unsigned char> prop_sampled;
		// temp
		OpenMesh::VPropHandleT<float> prop_curvature;
		// temp
		std::vector<saliency_prop_t> prop_saliency_levels;
		// cleanup to be run from result
		std::function<void(bool)> cleanup;
	};

	enum class saliency_state {
		idle, curv, area, nhprep, cand, run_full, run_sub, merge, norm, done, cancelled
	};

	struct saliency_progress {
		struct per_level {
			// progress in [0,1]
			float completion = 0;
			// number of samples taken
			int completed_samples = 0;
			// was (is) subsampling active?
			bool subsampled = false;
			// was (is) the normalmap filter active?
			bool normalmap_filter = false;
		};

		// this must be sized correctly before starting and not resized while running
		std::vector<per_level> levels;
		std::chrono::milliseconds elapsed_time{0};
		int completed_levels = 0;		
		int total_vertices = 0;
		bool should_cancel = false;
		saliency_state state = saliency_state::idle;
	};

	saliency_result compute_saliency(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress);

	std::future<saliency_result> compute_saliency_async(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress);

}

#endif
