
/*
 * Copyright 2020
 * Computational Media Innovation Centre
 * Victoria University of Wellington
 *
 */

#include "saliency.hpp"

#include <thread>
#include <chrono>
#include <algorithm>
#include <iostream>

#include <omp.h>

#include "plf_colony.h"

#include "model.hpp"
#include "neighborhood.hpp"
#include "curvature.hpp"
#include "meshutils.hpp"

namespace green {

	std::future<bool> compute_saliency_async(Model &model, const saliency_params &params, saliency_progress &progress) {
		return std::async([&, params=params]() { return compute_saliency(model, params, progress); });
	}

	class SaliencyComputation {
	private:
		struct sample_candidate {
			unsigned vdi = -1;
			float summedarea = 0;
			bool deleted = false;
		};

		saliency_params m_params;
		saliency_progress &m_progress;
		TriMesh &m_mesh;
		MeshCache m_meshcache;
		std::vector<sample_candidate> m_candidates0;
		OpenMesh::VPropHandleT<float> m_vertexAreasProperty;
		OpenMesh::EPropHandleT<float> m_edgeLengthProperty;
		OpenMesh::VPropHandleT<float> m_curvatureProperty;
		OpenMesh::VPropHandleT<float> m_saliencyProperty;
		OpenMesh::VPropHandleT<std::vector<float>> m_multiscaleSaliencyProperty;
		float m_surfaceArea = 0;
		float m_hMin = 0;
		float m_hMax = 0;
		std::vector<CalculationStats> m_thread_stats;
		std::chrono::steady_clock::time_point m_time_start = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point m_time_last_percent = std::chrono::steady_clock::now();

	public:
		SaliencyComputation(const saliency_params &params_, saliency_progress &progress_, TriMesh &mesh_)
			: m_params(params_)
			, m_progress(progress_)
			, m_mesh(mesh_)
		{
			std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;
			m_thread_stats.assign(omp_get_max_threads(), {});

			if (m_mesh.get_property_handle(m_saliencyProperty, "quality")) {
				// reuse saliency property
			} else {
				m_mesh.add_property(m_saliencyProperty);
			}

			std::cout << "Computing vertex areas" << std::endl;
			m_vertexAreasProperty = computeVertexAreas(m_mesh);

			std::cout << "Computing edge lengths" << std::endl;
			m_edgeLengthProperty = computeEdgeLengths(m_mesh);

			lce::Histogram curvhist;

			std::cout << "Computing curvature" << std::endl;
			// TODO params and curv selection
			computeDoNMaxDiffs(m_mesh, m_curvatureProperty, curvhist, m_vertexAreasProperty, 1.f);

			m_hMin = curvhist.getMin();
			m_hMax = curvhist.getMax();

			std::cout << "Computing surface area" << std::endl;
			m_surfaceArea = surfaceArea(m_mesh);

			m_mesh.add_property(m_multiscaleSaliencyProperty);
			for (auto vIt = m_mesh.vertices_begin(), vEnd = m_mesh.vertices_end(); vIt != vEnd; ++vIt) {
				m_mesh.property(m_multiscaleSaliencyProperty, *vIt).resize(m_params.levels);
			}

			m_meshcache = MeshCache(m_mesh, m_edgeLengthProperty, m_vertexAreasProperty, m_curvatureProperty);

			std::cout << "Preparing subsampling candidates" << std::endl;
			m_candidates0.reserve(m_meshcache.vdis.size());
			std::transform(m_meshcache.vdis.begin(), m_meshcache.vdis.end(), std::back_inserter(m_candidates0), [](auto &vdi) { return sample_candidate{vdi, 0.f}; });
			{
				float aa = 0;
				for (auto &cand : m_candidates0) {
					const float a = m_meshcache.get_vertex(cand.vdi).props[MeshCache::vertex_prop_area];
					cand.summedarea = aa;
					aa += a;
				}
				for (auto &cand : m_candidates0) {
					TriMesh::VertexHandle v(m_meshcache.get_vertex(cand.vdi).vi);
					// normalize summed areas
					cand.summedarea /= aa;
				}
			}
		}

		void run() {

			// TODO cancellation

			m_progress.completed_levels = 0;
			m_progress.completed_vertices = 0;
			m_progress.total_vertices = m_mesh.n_vertices();

			const float MaxRadius = sqrt(m_surfaceArea * m_params.area / 3.14159265f);

			// TODO params
			const float samples_per_neighbourhood = 100;
			std::cout << "Saliency samples per neighbourhood: " << samples_per_neighbourhood << std::endl;

			// compute saliency at multiple levels
			for (int currentLevel = 0; currentLevel < m_params.levels; currentLevel++)
			{
				m_progress.completed_vertices = 0;

				const float currentRadius = MaxRadius / pow(2.0f, static_cast<float>(currentLevel));
				const float currentArea = currentRadius * currentRadius * 3.14159265f;

				// upper bound because infinities and overflow and we don't support more than INT_MAX vertices atm anyway
				const int subsampling = samples_per_neighbourhood > 0
					? 1.5f + std::min((currentArea / m_surfaceArea) * m_mesh.n_vertices() / samples_per_neighbourhood, 1000000.f)
					: 1;

				const bool normalmap_filter = m_params.normalmap_filter && (currentRadius * currentRadius * 3.14159265f / m_surfaceArea) <= 0.00125f;

				std::cout << "Desired saliency subsampling ~" << (100.f / subsampling) << "% (~" << subsampling << "x)" << std::endl;

				if (subsampling >= 5 && samples_per_neighbourhood > 0) {
					run_level_subsampled(currentLevel, currentRadius, subsampling, normalmap_filter);
				} else {
					run_level_full(currentLevel, currentRadius, normalmap_filter);
				}

				m_progress.completed_levels = currentLevel + 1;
			}

			CalculationStats stats;
			for (auto &ts : m_thread_stats) {
				stats.merge(ts);
			}
			stats.dump_stats(omp_get_max_threads());

			// merge saliency values to a single score / property
			std::cout << "Merging saliency values" << std::endl;
			for (auto vIt = m_mesh.vertices_begin(), vEnd = m_mesh.vertices_end(); vIt != vEnd; ++vIt)
			{
				float s = 0.0f;

				for (unsigned i = 0; i < m_params.levels; ++i)
				{
					const std::vector<float> & saliencyValues = m_mesh.property(m_multiscaleSaliencyProperty, *vIt);
					assert(saliencyValues.size() == m_params.levels);
					//s = std::max(s, saliencyValues[i]);
					s += saliencyValues[i];
				}

				s /= m_params.levels;

				m_mesh.property(m_saliencyProperty, *vIt) = s;
			}

			// normalize saliency and add curvature
			std::cout << "Normalizing saliency values" << std::endl;
			float sMin =  FLT_MAX;
			float sMax = -FLT_MAX;
			for (auto vIt = m_mesh.vertices_begin(), vEnd = m_mesh.vertices_end(); vIt != vEnd; ++vIt)
			{
				float s = m_mesh.property(m_saliencyProperty, *vIt);
				sMin = std::min(s, sMin);
				sMax = std::max(s, sMax);
			}
			std::cout << "Raw saliency range: " << sMin << ' ' << sMax << std::endl;
			float sRange = sMax - sMin;
			for (auto vIt = m_mesh.vertices_begin(), vEnd = m_mesh.vertices_end(); vIt != vEnd; ++vIt)
			{
				float curvature = m_mesh.property(m_curvatureProperty, *vIt);
				float s = m_mesh.property(m_saliencyProperty, *vIt);
				float normalizedSal = (s - sMin) / sRange;
				normalizedSal = (normalizedSal + (m_params.curv_weight * curvature)) < 1 ? normalizedSal + (m_params.curv_weight * curvature) : 1;
				m_mesh.property(m_saliencyProperty, *vIt) = normalizedSal;
			}

			auto now = std::chrono::steady_clock::now();
			std::cout << "Saliency computation finished" << std::endl;
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(now - m_time_start);
			std::cout << "Overall: " << (m_progress.elapsed_time / std::chrono::duration<double>(1.0)) << "s" << std::endl;

			// TODO colorization

			// make sure the property is written to the file
			m_mesh.property(m_saliencyProperty).set_persistent(true);

		}

		void run_level_full(int currentLevel, float currentRadius, bool normalmap_filter) {

			std::atomic<int> completion{0};

			std::cout << "Saliency computation [full] (lv " << (currentLevel + 1) << "/" << m_params.levels << "): 0%";
			std::cout.flush();	

#pragma omp parallel for schedule(dynamic, 2)
			for (int i = 0; i < m_mesh.n_vertices(); i += simd_traits::simd_size)
			{
				auto &stats = m_thread_stats[omp_get_thread_num()];
				stats.timer_begin();

				//static thread_local std::vector<simd_traits::index_t> neighbors;
				//neighbors.clear();

				std::array<unsigned, simd_traits::simd_size> vdis;
				for (int j = 0; j < vdis.size(); j++) {
					vdis[j] = i + j < m_mesh.n_vertices() ? m_meshcache.vdis[i + j] : -1;
				}

				// NOTE: now produces vertex data indices, not ordinary vertex indices
				//getGeodesicNeighborhood(meshcache, stats, vdis, currentRadius, neighbors);
				const auto sal = getGeodesicNeighborhoodSaliency(m_meshcache, stats, vdis, currentRadius, m_hMin, m_hMax, normalmap_filter);

				stats.nh_timer_end();

				for (int j = 0; j < simd_traits::simd_size; j++) {
					if (vdis[j] == -1) break;

					TriMesh::VertexHandle v(m_meshcache.get_vertex(vdis[j]).vi);

					//const float sal = computeSaliency(j, meshcache, neighbors, zero1, zero9);

					using simd::simd_extract;
					m_mesh.property(m_multiscaleSaliencyProperty, v)[currentLevel] = sal[j];

				}

				//stats.sal_timer_end();

				const int completion1 = completion.fetch_add(simd_traits::simd_size);

				if (m_params.show_progress && omp_get_thread_num() == 0 && (stats.time0 - m_time_last_percent) > std::chrono::milliseconds(50)) {
					m_time_last_percent = stats.time0;
					auto pc = 100.0f * (float(completion1) / m_mesh.n_vertices());
					std::printf("\rSaliency computation [full] (lv %u/%d): %7.3f%%", currentLevel + 1, m_params.levels, pc);
					m_progress.completed_vertices = completion1;
				}
			}

			std::printf("\rSaliency computation [full] (lv %u/%d): %7.3f%%\n", currentLevel + 1, m_params.levels, 100.0f * (float(completion) / m_mesh.n_vertices()));
			std::cout << "Actual saliency subsampling: 100% (1x)" << std::endl;

		}

		void run_level_subsampled(int currentLevel, float currentRadius, int subsampling, bool normalmap_filter) {
			
			std::atomic<int> completion{0};

			// should subsample
			std::cout << "Saliency computation [subsampled] (lv " << (currentLevel + 1) << "/" << m_params.levels << "): 0%";
			std::cout.flush();

			struct weighted_saliency {
				std::atomic<float> s{0.f};
				std::atomic<float> w{0.f};
			};

			std::vector<weighted_saliency> tempSaliencyProperty(m_mesh.n_vertices());

			// radius of a circle with 1/nsamples of the surface area
			const float subsampling_radius = sqrt(m_surfaceArea * subsampling / m_mesh.n_vertices() / 3.14159265f);

			// correction factor to (experimentally) try to get closer to the desired number of samples.
			// we would expect this to be < 2, because circles of subsampling_radius would overlap on a flat surface.
			// the actual surface area of the geodesic r-neighbourhood roughly decreases with increasing curvature,
			// so more complex models will generally result in taking more samples than intended.
			// (geodesic r-neighbourhood area probably isn't bounded in reality, and can exceed pi*r^2).
			// TODO this may depend on degree of subsampling
			// TODO this depends on parallelism! (but not necessarily predictably)
			// this is also affected by the neighbourhood distance measurement, which is calculated
			// as a sum of edge lengths and so is an upper bound on the real distance.
			constexpr float sampling_correction = 1.8f;

			// radius within which to prevent future samples
			const float exclusion_radius = subsampling_radius * sampling_correction;

			std::vector<plf::colony<sample_candidate>::iterator> candidateProperty(m_mesh.n_vertices());

			plf::colony<sample_candidate> candidates;
			candidates.reserve(m_meshcache.vdis.size());
			for (auto &cand : m_candidates0) {
				TriMesh::VertexHandle v(m_meshcache.get_vertex(cand.vdi).vi);
				// after insertion, candidates colony starts in the same order
				auto it = candidates.insert(cand);
				// store the colony iterator as a vertex property
				// colony iterators do not get invalidated by erasures
				candidateProperty[v.idx()] = it;
			}

			struct next_candidate {
				plf::colony<sample_candidate>::iterator it;
				// tracking the summed-area start point separately from the iterator
				// should be more robust to area differences between vertices
				float summedarea = 0;
			};

			auto cand_inc_loop = [&](auto &it) {
				if (++it == candidates.end()) {
					it = candidates.begin();
					return true;
				}
				return false;
			};

			auto cand_inc_loop_area = [&](next_candidate &nc, float a) {
				const auto it0 = nc.it;
				float a1 = nc.summedarea + a;
				while (nc.it->summedarea < a1 || nc.it->deleted) {
					if (cand_inc_loop(nc.it)) a1 = std::max(0.f, a1 - 1.f);
					if (nc.it == it0) {
						// no more candidates
						// TODO if the only candidate(s) are between it0 and a1 we won't find them
						// this is not really a problem because the main sample loop will move
						// the candidate iterators to valid candidates until they're exhausted
						return true;
					}
				}
				nc.summedarea = a1;
				return false;
			};

			std::vector<next_candidate> thread_next_cand(omp_get_max_threads());
			for (auto &nc : thread_next_cand) {
				auto i = &nc - thread_next_cand.data();
				nc.it = candidates.begin();
				cand_inc_loop_area(nc, float(i) / omp_get_max_threads());
			}

			// note: this isnt strictly random anymore
			auto random_vdi = [&]() {
				if (candidates.empty()) return unsigned(-1);
				auto &nc = thread_next_cand[omp_get_thread_num()];
				if (nc.it == plf::colony<sample_candidate>::iterator()) {
					// candidates already exhausted
					return unsigned(-1);
				}
				if (!nc.it->deleted) {
					// candidate already usable
					return nc.it->vdi;
				}
				if (cand_inc_loop_area(nc, float(subsampling) / m_mesh.n_vertices())) {
					// search for next candidate (probably) exhausted
					// note: the 'probably' currently prevents us from sentinel-izing the iterator
					return unsigned(-1);
				}
				// note: should not use the deleted flag as a reason to return -1
				// deleted flag may be set by another thread after we selected a valid candidate
				return nc.it->vdi;
			};

			while (true) {

				// specifying dynamic scheduling seems very important for this to perform well
				// note: too few samples per thread per will hurt performance (with openmp overhead)
				// while too many samples per thread will result in worse sampling characteristics
				// due to lazy candidate deletion
#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < 50 * omp_get_max_threads(); i++)
				{
					auto &stats = m_thread_stats[omp_get_thread_num()];
					stats.timer_begin();

					// yes, this may select candidates that have been (or will be) flagged as deleted
					// we have to put up with this in order to parallelize
					const auto rootvdi = random_vdi();
					if (rootvdi == -1) {
						// can't break an openmp loop, so we have to spin
						// this should only happen when we're out of sample candidates
						// (but can currently happen spuriously in rare conditions)
						std::this_thread::yield();
						continue;
					}

					const auto atomic_float_accum = [](std::atomic<float> &a, float b) {
						float a0 = 0;
						do {
							a0 = a.load(std::memory_order_relaxed);
						} while (!a.compare_exchange_weak(a0, a0 + b, std::memory_order_relaxed));
					};

					const auto visitor = [&](unsigned vdi, float r, float s) {
						TriMesh::VertexHandle v(m_meshcache.get_vertex(vdi).vi);
						//const float w = rootvdi == vdi;
						//const float w = r < exclusion_radius * 0.5f;
						// approx inverse distance, preventing inf
						// weight falls to some constant at exclusion radius, and zero at 2x that
						// zero is important to avoid discontinuity at edge of distribution radius
						// weight scale is not important
						const float re = exclusion_radius;
						const float b = 0.001f;
						const float w = std::max(re / (r + re * b) - (1 / (b + 2)), 0.f);
						// TODO try gaussian weighting?
						atomic_float_accum(tempSaliencyProperty[v.idx()].s, s * w);
						atomic_float_accum(tempSaliencyProperty[v.idx()].w, w);
						if (r < exclusion_radius) {
							// exclude from future sampling
							const auto none = plf::colony<sample_candidate>::iterator();
							auto it = candidateProperty[v.idx()];
							if (it != none) {
								candidateProperty[v.idx()] = none;
								it->deleted = true;
								// cannot do parallel erase!
								//candidates.erase(it);
							}
						}
					};

					// could be up to exclusion_radius * 2 distance between samples, so distribute over that to ensure coverage
					// need to go slightly more to avoid problems due to vertex discretization (to allow weight to reach zero)
					subsampleGeodesicNeighborhoodSaliency(m_meshcache, stats, rootvdi, currentRadius, m_hMin, m_hMax, normalmap_filter, exclusion_radius * 2.1f, visitor);

					stats.nh_timer_end();

					// TODO ? completion is approximate when subsampling
					const int completion1 = completion.fetch_add(subsampling);

					if (m_params.show_progress && omp_get_thread_num() == 0 && (stats.time0 - m_time_last_percent) > std::chrono::milliseconds(50)) {
						m_time_last_percent = stats.time0;
						auto pc = 100.0f * (float(completion1) / m_mesh.n_vertices());
						std::printf("\rSaliency computation [subsampled] (lv %u/%d): %7.3f%%", currentLevel + 1, m_params.levels, pc);
						m_progress.completed_vertices = completion1;
					}
				}

				// make thread candidate iterators point to non-deleted candidates
				for (auto &nc : thread_next_cand) {
					auto it0 = nc.it;
					while (nc.it->deleted) {
						cand_inc_loop(nc.it);
						if (nc.it == it0) {
							nc.it = plf::colony<sample_candidate>::iterator();
							break;
						}
					}
				}

				// actually erase candidates after some number of samples
				for (auto it = candidates.begin(); it != candidates.end(); ) {
					if (it->deleted) {
						it = candidates.erase(it);
					} else {
						it++;
					}
				}

				if (candidates.empty()) {
					std::printf("\rSaliency computation [subsampled] (lv %u/%d): %7.3f%% : no more sample candidates\n", currentLevel + 1, m_params.levels, 100.0f * (float(completion) / m_mesh.n_vertices()));
					const float actual_subsampling = m_mesh.n_vertices() / float(completion / subsampling);
					std::cout << "Actual saliency subsampling: " << (100.f / actual_subsampling) << "% (" << actual_subsampling << "x)" << std::endl;
					break;
				}

			}

			for (auto vIt = m_mesh.vertices_begin(), vEnd = m_mesh.vertices_end(); vIt != vEnd; ++vIt) {
				const auto &p = tempSaliencyProperty[vIt->idx()];
				m_mesh.property(m_multiscaleSaliencyProperty, *vIt)[currentLevel] = p.s / (p.w + 0.001f);
			}

		}

		~SaliencyComputation() {
			// TODO property reuse
			m_mesh.remove_property(m_vertexAreasProperty);
			m_mesh.remove_property(m_edgeLengthProperty);
			m_mesh.remove_property(m_curvatureProperty);
			m_mesh.remove_property(m_multiscaleSaliencyProperty);
		}
	};
	
	bool compute_saliency(Model &model, const saliency_params &params, saliency_progress &progress) {
		// TODO thread control
		omp_set_num_threads(2);
		const auto time_start = std::chrono::steady_clock::now();
		auto &mesh = model.trimesh();
		SaliencyComputation s(params, progress, mesh);
		s.run();
		return true;
	}

}
