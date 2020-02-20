
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
#include <vector>
#include <functional>

#include <imgui.h>

#include "meshutils.hpp"
#include "Histogram.h"

namespace green {

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
		bool normalmap_filter = false;
		// enable manual subsampling (takes precedence)
		bool subsample_manual = false;
		// enable auto subsampling
		bool subsample_auto = true;
		// command line progress output
		bool show_progress = true;
	};

	struct saliency_mesh_params {
		// TODO const
		TriMesh *mesh = nullptr;
		// inputs
		OpenMesh::VPropHandleT<float> prop_vertex_area;
		OpenMesh::EPropHandleT<float> prop_edge_length;
		// outputs
		OpenMesh::VPropHandleT<float> prop_curvature;
		OpenMesh::VPropHandleT<float> prop_saliency;
		std::vector<OpenMesh::VPropHandleT<float>> prop_saliency_levels;
		std::function<void(bool)> cleanup;
	};

	struct saliency_progress {
		struct per_level {
			int completed_vertices = 0;
			int desired_subsampling = 1;
			bool subsampled = false;
		};

		// this must be sized correctly before starting and not resized while running
		std::vector<per_level> levels;
		std::chrono::milliseconds elapsed_time{0};
		int completed_levels = 0;		
		int total_vertices = 0;
		bool should_cancel = false;
		
	};

	saliency_result compute_saliency(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress);

	std::future<saliency_result> compute_saliency_async(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress);

}

#endif
