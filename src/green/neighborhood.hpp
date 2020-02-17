
/*
 * Copyright 2020
 * Computational Media Innovation Centre
 * Victoria University of Wellington
 *
 */

#pragma once

#ifndef GREEN_NEIGHBORHOOD_HPP
#define GREEN_NEIGHBORHOOD_HPP

#include <vector>
#include <array>
#include <functional>
#include <chrono>

#include "simd.h"
#include "model.hpp"

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __GNUC__
 // this'll do for now i guess
#define __popcnt __builtin_popcount
#endif

#define USE_VERTEX_AREA_WEIGHTING

namespace green {

	using simd::simd4_int;
	using simd::simd4_float;
	using simd::simd8_int;
	using simd::simd8_float;

	struct simd1_traits {
		static constexpr int simd_size = 1;
		using mask_t = bool;
		using index_t = int;
		using dist_t = float;
	};

	struct simd4_traits {
		static constexpr int simd_size = 4;
		using mask_t = simd4_int;
		using index_t = simd4_int;
		using dist_t = simd4_float;
	};

	struct simd8_traits {
		static constexpr int simd_size = 8;
		using mask_t = simd8_int;
		using index_t = simd8_int;
		using dist_t = simd8_float;
	};

	// simd8 has minimal gain over simd4 when multithreaded
	using simd_traits = simd4_traits;

	struct MeshCache {

		// max right shift for vertex data indices that still gives unique ids
		// vertex data is padded to ensure this
		// (this is checked on construction, just to be sure)
		static constexpr int vdi_shift = 4;

		enum {
			vertex_prop_area, vertex_prop_curv
		};

		struct edge {
			unsigned ndi; // data index, not vertex index!
			float cost;
		};

		struct vertex {
			int vi; // actual vertex index
			unsigned nedges;
			float props[2];
			OpenMesh::Vec3f pos, norm;
			// actual size == nedges (extends past struct end)
			edge edges[0];
		};
		std::vector<unsigned> vi2di;
		std::vector<unsigned> vdis;
		std::vector<unsigned> data;

		MeshCache() = default;

		MeshCache(
			const TriMesh &mesh,
			OpenMesh::EPropHandleT<float> edgeLengthProperty,
			OpenMesh::VPropHandleT<float> vertexAreasProperty,
			OpenMesh::VPropHandleT<float> curvatureMeasure
		);

		const vertex & get_vertex(unsigned vdi) const {
			return reinterpret_cast<const vertex &>(data[vdi]);
		}

	};

	struct CalculationStats {

		int64_t nh_loop_count = 0;
		int64_t nh_loop_occupancy = 0;
		int64_t nh_push_count = 0;
		int64_t nh_push_occupancy = 0;

		std::chrono::nanoseconds nh_duration{};
		std::chrono::nanoseconds sal_duration{};

		std::chrono::steady_clock::time_point time0;

		template <typename MaskT>
		void nh_record_loop(MaskT occupancy) {
			using simd::mask2bits;
			const int obits = mask2bits(occupancy);
			nh_loop_count++;
			nh_loop_occupancy += __popcnt(obits);
		}

		template <typename MaskT>
		void nh_record_push(MaskT occupancy) {
			using simd::mask2bits;
			const int obits = mask2bits(occupancy);
			nh_push_count++;
			nh_push_occupancy += __popcnt(obits);
		}

		void merge(const CalculationStats &other) {
			nh_loop_count += other.nh_loop_count;
			nh_loop_occupancy += other.nh_loop_occupancy;
			nh_push_count += other.nh_push_count;
			nh_push_occupancy += other.nh_push_occupancy;
			nh_duration += other.nh_duration;
			sal_duration += other.sal_duration;
		}

		void timer_begin() {
			time0 = std::chrono::steady_clock::now();
		}

		void nh_timer_end() {
			auto now = std::chrono::steady_clock::now();
			nh_duration += now - time0;
			time0 = now;
		}

		void sal_timer_end() {
			auto now = std::chrono::steady_clock::now();
			sal_duration += now - time0;
			time0 = now;
		}

		void dump_stats(int threadcount) const;

	};

	void getGeodesicNeighborhood(
		const MeshCache &mesh,
		CalculationStats &stats,
		const std::array<unsigned, simd_traits::simd_size> &vdis,
		float radius,
		std::vector<simd_traits::index_t> & neighbors
	);

	std::array<float, simd_traits::simd_size> getGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		const std::array<unsigned, simd_traits::simd_size> &rootvdis,
		float radius,
		float curvmin,
		float curvmax,
		bool minimizeSmallChanges
	);

	void subsampleGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		unsigned rootvdi,
		float radius,
		float curvmin,
		float curvmax,
		bool minimizeSmallChanges,
		float distribution_radius,
		const std::function<void(unsigned vdi, float r, float s)> &distribute_saliency
	);

}

#endif
