
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

		static constexpr unsigned ncurvbins = 256;

		// min vertex size (in unsigned alloc units)
		// vertex data is padded to ensure this
		// (this is checked on construction, just to be sure)
		// reasonable values: {10,12,14,16}
		static constexpr unsigned vert_min_size = 12;

#pragma pack(push)
#pragma pack(1)
		struct edge {
			// data index
			// note: while we could theoretically store this more compactly (in most cases)
			// as an offset from the current vertex, in practice this is extremely difficult
			unsigned ndi;
			// fixed-point sqrt cost in [0,1] (* cost_scale)
			uint16_t costx;
		};
#pragma pack(pop)

		struct vertex {
			// 65535 edge limit
			uint16_t nedges;
			// binned curvature
			// currently only using 1 byte, but more could be useful in future
			uint16_t curvbin;
			// fixed-point sqrt area in [0,1] (* area_scale)
			uint16_t areax;
			// actual size == nedges (extends past struct end)
			edge edges[0];
		};

		struct vertex_aux {
			OpenMesh::Vec3f pos;
			OpenMesh::Vec3f norm;
			// actual vertex index
			int vi;
		};
		
		std::vector<unsigned> vi2di;
		std::vector<unsigned> vdis;

		// TODO would be less naughty if these were vector<std::byte>
		std::vector<unsigned> data;
		std::vector<unsigned> dataaux;

		float area_scale = 1;
		float cost_scale = 1;

		MeshCache() = default;

		MeshCache(
			const PolyMesh &mesh,
			OpenMesh::EPropHandleT<float> edgeLengthProperty,
			OpenMesh::VPropHandleT<float> vertexAreasProperty,
			OpenMesh::VPropHandleT<float> curvatureMeasure
		);

		unsigned vdi_uid(unsigned vdi) const noexcept {
			// note: compiler optimization turns this into a bitshift or reciprocal multiplication
			return vdi / vert_min_size;
		}

		vertex & get_vertex(unsigned vdi) noexcept {
			return reinterpret_cast<vertex &>(data[vdi]);
		}

		const vertex & get_vertex(unsigned vdi) const noexcept {
			return reinterpret_cast<const vertex &>(data[vdi]);
		}

		vertex_aux & get_vertex_aux(unsigned vdi) noexcept {
			return reinterpret_cast<vertex_aux &>(dataaux[vdi_uid(vdi) * (sizeof(vertex_aux) / sizeof(unsigned))]);
		}

		const vertex_aux & get_vertex_aux(unsigned vdi) const noexcept {
			return reinterpret_cast<const vertex_aux &>(dataaux[vdi_uid(vdi) * (sizeof(vertex_aux) / sizeof(unsigned))]);
		}

		void dump_to_file() const;

		uint16_t encode_area(float area) const noexcept {
			// prevent 0 area
			return std::max(unsigned(65535 * sqrt(area / area_scale)), 1u);
		}

		float decode_area(uint16_t areax) const noexcept {
			const float areaxf = areax * (1.f / 65535);
			return area_scale * (areaxf * areaxf);
		}

		uint16_t encode_cost(float cost) const noexcept {
			// prevent 0 cost
			return std::max(unsigned(65535 * sqrt(cost / cost_scale)), 1u);
		}

		float decode_cost(uint16_t costx) const noexcept {
			const float costxf = costx * (1.f / 65535);
			return cost_scale * (costxf * costxf);
		}
	};

	struct CalculationStats {

		int64_t nh_loop_count = 0;
		int64_t nh_loop_occupancy = 0;
		int64_t nh_push_count = 0;
		int64_t nh_push_occupancy = 0;
		int64_t nh_search_count = 0;
		int64_t nh_big_search_count = 0;

		std::chrono::nanoseconds nh_duration{};
		std::chrono::nanoseconds sal_duration{};

		std::chrono::steady_clock::time_point time0;

		template <typename MaskT>
		void nh_record_loop(MaskT occupancy) noexcept {
			using simd::mask2bits;
			const int obits = mask2bits(occupancy);
			nh_loop_count++;
			nh_loop_occupancy += __popcnt(obits);
		}

		template <typename MaskT>
		void nh_record_push(MaskT occupancy) noexcept {
			using simd::mask2bits;
			const int obits = mask2bits(occupancy);
			nh_push_count++;
			nh_push_occupancy += __popcnt(obits);
		}

		template <int Phase>
		void nh_record_search() noexcept {
			if constexpr (Phase == 1) nh_search_count++;
		}

		template <int Phase>
		void nh_record_big_search() noexcept {
			if constexpr (Phase == 1) nh_big_search_count++;
		}

		void merge(const CalculationStats &other) {
			nh_loop_count += other.nh_loop_count;
			nh_loop_occupancy += other.nh_loop_occupancy;
			nh_push_count += other.nh_push_count;
			nh_push_occupancy += other.nh_push_occupancy;
			nh_search_count += other.nh_search_count;
			nh_big_search_count += other.nh_big_search_count;
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
		float noise_height
	);

	float getGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		unsigned rootvdi,
		float radius,
		float noise_height
	);

	void subsampleGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		unsigned rootvdi,
		float radius,
		float noise_height,
		float distribution_radius,
		const std::function<void(unsigned vdi, float r, float s)> &distribute_saliency
	);

}

#endif
