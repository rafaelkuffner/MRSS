
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
#include <random>
#include <sstream>

#include <omp.h>

#include "plf_colony.h"

#include "model.hpp"
#include "neighborhood.hpp"
#include "curvature.hpp"
#include "meshutils.hpp"

namespace green {

	std::string saliency_user_params::str(bool verbose) const {
		// TODO locale independent formatting - this needs to be parseable
		// we dont set a locale atm though
		char buf[128];
		char *end = buf + sizeof(buf);
		char *p = buf;
		if (verbose) {
			// use extra precision
			p += snprintf(p, end - p, "l=%d,a=%g,c=%g,r=%g", levels, area, curv_weight, normal_power);
			if (normalmap_filter) p += snprintf(p, end - p, ",e=%g", noise_height);
		} else {
			p += snprintf(p, end - p, "l=%d,a=%.3g,c=%.2g,r=%.3g", levels, area, curv_weight, normal_power);
			if (normalmap_filter) p += snprintf(p, end - p, ",e=%.3g", noise_height);
		}
		if (subsample_auto) p += snprintf(p, end - p, ",s=%.0f", samples_per_neighborhood);
		return {buf};
	}

	bool saliency_user_params::parse(std::string_view s) {
		std::cout << "parsing saliency param string " << s << std::endl;
		// start with defaults, and construct a separate instance in case of failure
		saliency_user_params uparams;
		// ensure these switches start off, because they are implied on by value params
		uparams.normalmap_filter = false;
		uparams.subsample_manual = false;
		uparams.subsample_auto = false;
		std::istringstream iss{std::string(s)};
		while (iss) {
			const char sp = iss.get();
			float val = 0;
			bool has_val = false;
			auto need_val = [&]() {
				if (!has_val) {
					std::cout << "value needed for saliency param " << sp << std::endl;
					return false;
				}
				return true;
			};
			if (iss.peek() == '=') {
				iss.get();
				has_val = true;
				// TODO use locale independent parsing
				iss >> val;
				if (!iss) {
					std::cout << "error parsing saliency param " << sp << std::endl;
					return false;
				}
			} 
			if (iss.peek() == ',') {
				iss.get();
			}
			switch (sp) {
			case 'l':
				if (!need_val()) return false;
				uparams.levels = val;
				break;
			case 'a':
				if (!need_val()) return false;
				uparams.area = val;
				break;
			case 'c':
				if (!need_val()) return false;
				uparams.curv_weight = val;
				break;
			case 'r':
				if (!need_val()) return false;
				uparams.normal_power = val;
				break;
			case 'e':
				if (!need_val()) return false;
				uparams.normalmap_filter = true;
				uparams.noise_height = val;
				break;
			case 's':
				if (!need_val()) return false;
				uparams.subsample_auto = true;
				uparams.samples_per_neighborhood = val;
				break;
			default:
				std::cout << "unknown saliency param " << sp << std::endl;
				return false;
			}
		}
		*this = uparams;
		return true;
	}

	void saliency_user_params::sanitize() {
		using std::clamp;
		using std::min;
		using std::max;
		levels = clamp(levels, 1, 10);
		area = clamp(area, 0.f, 1.f);
		curv_weight = max(curv_weight, 0.f);
		normal_power = clamp(normal_power, 0.f, 10.f);
		subsampling_rate = max(subsampling_rate, 1.f);
		samples_per_neighborhood = max(samples_per_neighborhood, 1.f);
		noise_height = clamp(noise_height, 0.f, 1.f);
		thread_count = max(thread_count, 0);
	}
	
	std::future<saliency_result> compute_saliency_async(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress) {
		return std::async([=, &progress]() { return compute_saliency(mparams, uparams, progress); });
	}

	class SaliencyComputation {
	private:
		struct sample_candidate {
			unsigned vdi = -1;
			float area = 0;
			float summedarea = 0;
			bool deleted = false;
		};

		saliency_mesh_params m_mparams;
		saliency_user_params m_uparams;
		saliency_progress &m_progress;
		MeshCache m_meshcache;
		std::vector<sample_candidate> m_candidates0;
		float m_surfaceArea = 0;
		float m_total_vertex_area = 0;
		float m_real_noise_height = 0;
		float m_hMin = 0;
		float m_hMax = 0;
		std::vector<CalculationStats> m_thread_stats;
		std::chrono::steady_clock::time_point m_time_start = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point m_time_last_percent = std::chrono::steady_clock::now();

	public:
		SaliencyComputation(const saliency_mesh_params &mparams_, const saliency_user_params &uparams_, saliency_progress &progress_)
			: m_mparams(mparams_)
			, m_uparams(uparams_)
			, m_progress(progress_)
		{
			// ensure params are valid
			m_uparams.sanitize();
			
			std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;
			m_thread_stats.assign(omp_get_max_threads(), {});

			lce::Histogram curvhist;

			// note: we can't cache the curvature because we couldn't adjust the normal power etc
			m_progress.state = saliency_computation_state::curv;
			std::cout << "Computing curvature" << std::endl;
			// TODO curv selection param
			computeDoNMaxDiffs(*m_mparams.mesh, m_mparams.prop_curvature, curvhist, m_mparams.prop_vertex_area, m_uparams.normal_power);
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(std::chrono::steady_clock::now() - m_time_start);

			m_hMin = curvhist.getMin();
			m_hMax = curvhist.getMax();
			std::cout << "Curv min: " << m_hMin << std::endl;
			std::cout << "Curv max: " << m_hMax << std::endl;

			m_progress.state = saliency_computation_state::area;
			std::cout << "Computing surface area" << std::endl;
			m_surfaceArea = surfaceArea(*m_mparams.mesh);
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(std::chrono::steady_clock::now() - m_time_start);
			std::cout << "Surface area: " << m_surfaceArea << std::endl;

			m_real_noise_height = m_uparams.noise_height * sqrt(m_surfaceArea);
			std::cout << "Real noise height: " << m_real_noise_height << std::endl;

			m_progress.state = saliency_computation_state::nhprep;
			m_meshcache = MeshCache(*m_mparams.mesh, m_mparams.prop_edge_length, m_mparams.prop_vertex_area, m_mparams.prop_curvature);
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(std::chrono::steady_clock::now() - m_time_start);
			//m_meshcache.dump_to_file();

			m_progress.state = saliency_computation_state::cand;
			std::cout << "Preparing subsampling candidates" << std::endl;
			m_candidates0.reserve(m_meshcache.vdis.size());
			std::transform(m_meshcache.vdis.begin(), m_meshcache.vdis.end(), std::back_inserter(m_candidates0), [](auto &vdi) { return sample_candidate{vdi}; });

			// randomize sample candidates
			// makes subsampling more stable and less sensitive to parallelization
			std::minstd_rand rand{std::random_device{}()};
			std::shuffle(m_candidates0.begin(), m_candidates0.end(), rand);

			{
				float aa = 0;
				for (auto &cand : m_candidates0) {
					const float a = m_meshcache.get_vertex(cand.vdi).props[MeshCache::vertex_prop_area];
					cand.area = a;
					cand.summedarea = aa;
					aa += a;
				}
				m_total_vertex_area = aa;
				for (auto &cand : m_candidates0) {
					TriMesh::VertexHandle v(m_meshcache.get_vertex(cand.vdi).vi);
					// normalize areas
					cand.area /= aa;
					cand.summedarea /= aa;
				}
			}

			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(std::chrono::steady_clock::now() - m_time_start);
		}

		bool run() {

			if (m_progress.should_cancel) return false;

			m_progress.completed_levels = 0;
			m_progress.total_vertices = m_mparams.mesh->n_vertices();

			const float MaxRadius = sqrt(m_surfaceArea * m_uparams.area / 3.14159265f);

			const float samples_per_neighbourhood = m_uparams.subsample_manual
				? (m_uparams.area * m_mparams.mesh->n_vertices() / (m_uparams.subsampling_rate - 1))
				: (m_uparams.subsample_auto ? m_uparams.samples_per_neighborhood : 0);
			std::cout << "Saliency samples per neighbourhood: " << samples_per_neighbourhood << std::endl;

			// compute saliency at multiple levels
			for (int currentLevel = 0; currentLevel < m_uparams.levels; currentLevel++)
			{
				const float currentRadius = MaxRadius / pow(2.0f, static_cast<float>(currentLevel));
				const float currentArea = currentRadius * currentRadius * 3.14159265f;

				// upper bound because infinities and overflow and we don't support more than INT_MAX vertices atm anyway
				const int subsampling = samples_per_neighbourhood > 0
					? 1.5f + std::min((currentArea / m_surfaceArea) * m_mparams.mesh->n_vertices() / samples_per_neighbourhood, 1000000.f)
					: 1;

				// TODO 4 is a magic number, is it a good one?
				const bool normalmap_filter = m_uparams.normalmap_filter && currentRadius < m_real_noise_height * 4;
				// && (currentRadius * currentRadius * 3.14159265f / m_surfaceArea) <= 0.00125f;
				m_progress.levels[currentLevel].normalmap_filter = normalmap_filter;
				
				std::cout << "Real neighbourhood radius: " << currentRadius << std::endl;
				std::cout << "Desired saliency subsampling ~" << (100.f / subsampling) << "% (~" << subsampling << "x)" << std::endl;

				if (subsampling >= 5 && samples_per_neighbourhood > 0) {
					m_progress.levels[currentLevel].subsampled = true;
					m_progress.state = saliency_computation_state::run_sub;
					run_level_subsampled(currentLevel, currentRadius, subsampling, normalmap_filter);
				} else {
					m_progress.levels[currentLevel].subsampled = false;
					m_progress.state = saliency_computation_state::run_full;
					run_level_full(currentLevel, currentRadius, normalmap_filter);
				}

				m_progress.completed_levels = currentLevel + 1;
				if (m_progress.should_cancel) return false;
			}

			CalculationStats stats;
			for (auto &ts : m_thread_stats) {
				stats.merge(ts);
			}
			stats.dump_stats(omp_get_max_threads());

			// merge saliency values to a single score / property
			m_progress.state = saliency_computation_state::merge;
			std::cout << "Merging saliency values" << std::endl;
			for (auto vIt = m_mparams.mesh->vertices_begin(), vEnd = m_mparams.mesh->vertices_end(); vIt != vEnd; ++vIt)
			{
				float s = 0.0f;

				for (unsigned i = 0; i < m_uparams.levels; ++i)
				{
					s += m_mparams.mesh->property(m_mparams.prop_saliency_levels[i], *vIt);
				}

				s /= m_uparams.levels;

				m_mparams.mesh->property(m_mparams.prop_saliency, *vIt) = s;
			}

			// normalize saliency and add curvature
			m_progress.state = saliency_computation_state::norm;
			std::cout << "Normalizing saliency values" << std::endl;
			float sMin =  FLT_MAX;
			float sMax = -FLT_MAX;
			for (auto vIt = m_mparams.mesh->vertices_begin(), vEnd = m_mparams.mesh->vertices_end(); vIt != vEnd; ++vIt)
			{
				float s = m_mparams.mesh->property(m_mparams.prop_saliency, *vIt);
				sMin = std::min(s, sMin);
				sMax = std::max(s, sMax);
			}
			std::cout << "Raw saliency range: " << sMin << ' ' << sMax << std::endl;
			float sRange = sMax - sMin;
			for (auto vIt = m_mparams.mesh->vertices_begin(), vEnd = m_mparams.mesh->vertices_end(); vIt != vEnd; ++vIt)
			{
				float curvature = m_mparams.mesh->property(m_mparams.prop_curvature, *vIt);
				float s = m_mparams.mesh->property(m_mparams.prop_saliency, *vIt);
				float normalizedSal = (s - sMin) / sRange;
				normalizedSal = (normalizedSal + (m_uparams.curv_weight * curvature)) < 1 ? normalizedSal + (m_uparams.curv_weight * curvature) : 1;
				m_mparams.mesh->property(m_mparams.prop_saliency, *vIt) = normalizedSal;
			}

			auto now = std::chrono::steady_clock::now();
			std::cout << "Saliency computation finished" << std::endl;
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(now - m_time_start);
			std::cout << "Overall: " << (m_progress.elapsed_time / std::chrono::duration<double>(1.0)) << "s" << std::endl;
			m_progress.state = saliency_computation_state::done;

			return true;
		}

		void run_level_full(int currentLevel, float currentRadius, bool normalmap_filter) {

			std::atomic<int> completion{0};

			std::cout << "Saliency computation [full] (lv " << (currentLevel + 1) << "/" << m_uparams.levels << "): 0%";
			std::cout.flush();	

#pragma omp parallel for schedule(dynamic, 2)
			for (int i = 0; i < m_mparams.mesh->n_vertices(); i += simd_traits::simd_size)
			{
				// can't break openmp loop, spin instead
				if (m_progress.should_cancel) {
					std::this_thread::yield();
					continue;
				}
				
				auto &stats = m_thread_stats[omp_get_thread_num()];
				stats.timer_begin();

				//static thread_local std::vector<simd_traits::index_t> neighbors;
				//neighbors.clear();

				std::array<unsigned, simd_traits::simd_size> vdis;
				for (int j = 0; j < vdis.size(); j++) {
					vdis[j] = i + j < m_mparams.mesh->n_vertices() ? m_meshcache.vdis[i + j] : -1;
				}

				// NOTE: now produces vertex data indices, not ordinary vertex indices
				//getGeodesicNeighborhood(meshcache, stats, vdis, currentRadius, neighbors);
				const auto sal = getGeodesicNeighborhoodSaliency(m_meshcache, stats, vdis, currentRadius, m_hMin, m_hMax, normalmap_filter * m_real_noise_height);

				stats.nh_timer_end();

				for (int j = 0; j < simd_traits::simd_size; j++) {
					if (vdis[j] == -1) break;

					TriMesh::VertexHandle v(m_meshcache.get_vertex(vdis[j]).vi);

					//const float sal = computeSaliency(j, meshcache, neighbors, zero1, zero9);

					using simd::simd_extract;
					m_mparams.mesh->property(m_mparams.prop_saliency_levels[currentLevel], v) = sal[j];

					// record sampled vertices (which is all of them here)
					m_mparams.mesh->property(m_mparams.prop_sampled, v) = true;
				}

				//stats.sal_timer_end();

				const int completion1 = completion.fetch_add(simd_traits::simd_size);

				if (m_uparams.show_progress && omp_get_thread_num() == 0 && (stats.time0 - m_time_last_percent) > std::chrono::milliseconds(50)) {
					m_time_last_percent = stats.time0;
					auto pc = 100.0f * (float(completion1) / m_mparams.mesh->n_vertices());
					std::printf("\rSaliency computation [full] (lv %u/%d): %7.3f%%", currentLevel + 1, m_uparams.levels, pc);
					m_progress.levels[currentLevel].completed_samples = completion1;
					m_progress.levels[currentLevel].completion = float(completion1) / m_mparams.mesh->n_vertices();
					m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(stats.time0 - m_time_start);
				}
			}

			m_progress.levels[currentLevel].completed_samples = completion;
			m_progress.levels[currentLevel].completion = 1;
			std::printf("\rSaliency computation [full] (lv %u/%d): %7.3f%%\n", currentLevel + 1, m_uparams.levels, 100.0f * (float(completion) / m_mparams.mesh->n_vertices()));
			std::cout << "Actual saliency subsampling: 100% (1x)" << std::endl;

		}

		void run_level_subsampled(int currentLevel, float currentRadius, int subsampling, bool normalmap_filter) {
			
			std::atomic<int> completion{0};

			// should subsample
			//std::cout << "Saliency computation [subsampled] @%8.2fs (lv " << (currentLevel + 1) << "/" << m_uparams.levels << "): 0%";
			std::cout.flush();
			std::printf("Saliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, 0);

			struct weighted_saliency {
				std::atomic<float> s{0.f};
				std::atomic<float> w{0.f};
			};

			std::vector<weighted_saliency> tempSaliencyProperty(m_mparams.mesh->n_vertices());

			// radius of a circle with 1/nsamples of the surface area
			const float subsampling_radius = sqrt(m_surfaceArea * subsampling / m_mparams.mesh->n_vertices() / 3.14159265f);

			std::vector<plf::colony<sample_candidate>::iterator> candidateProperty(m_mparams.mesh->n_vertices());

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
				if (cand_inc_loop_area(nc, float(subsampling) / m_mparams.mesh->n_vertices())) {
					// search for next candidate (probably) exhausted
					// note: the 'probably' currently prevents us from sentinel-izing the iterator
					return unsigned(-1);
				}
				// note: should not use the deleted flag as a reason to return -1
				// deleted flag may be set by another thread after we selected a valid candidate
				return nc.it->vdi;
			};

			// correction factor to (experimentally) try to get closer to the desired number of samples.
			// we would expect this to be < 2, because circles of subsampling_radius would overlap on a flat surface.
			// the actual surface area of the geodesic r-neighbourhood roughly decreases with increasing curvature,
			// so more complex models will generally result in taking more samples than intended.
			// (geodesic r-neighbourhood area probably isn't bounded in reality, and can exceed pi*r^2).
			// TODO this may depend on degree of subsampling
			// TODO this depends on parallelism! (less so with candidate randomization, and also smaller value)
			// this is also affected by the neighbourhood distance measurement, which is calculated
			// as a sum of edge lengths and so is an approximation of the real distance.
			// NOTE should not initialize sampling correction too high and adjust down because exclusion is irreversible.
			// NOTE using a constant value is desirable rather than trying to aim for a specific number of samples.
			// we want to control the distance between samples (rather than the area 'owned' by a sample)
			// because that is what should determine the interpolation error.
			// NOTE the number of samples is then incidental and the 'overcompletion' is the algorithm operating as desired.
			// NOTE 1.4 is not scientific and may not be quite low enough in all scenarios but is generally ok.
			const float sampling_correction = 1.4f;

			// area of non-excluded vertices
			// can use this to report actual progress %
			float remaining_area = 1;

			// sample in spits until there are no more candidates
			while (!m_progress.should_cancel) {

				// radius within which to prevent future samples
				const float exclusion_radius = subsampling_radius * sampling_correction;

				// adjust samples per spit based on vertices covered per sample
				const float vertices_per_spit_per_thread = 1000000;
				const float vertices_per_sample = m_mparams.mesh->n_vertices() * currentRadius * currentRadius * 3.14159265f / m_surfaceArea;

				// specifying dynamic scheduling seems very important for this to perform well
				// note: too few samples per thread per will hurt performance (with openmp overhead)
				// while too many samples per thread will also hurt performance due to lazy candidate
				// deletion and may result in worse sampling characteristics.
				const int samples_per_spit = std::min<int>(omp_get_max_threads() * vertices_per_spit_per_thread / vertices_per_sample, m_mparams.mesh->n_vertices() / float(subsampling) * 0.5f) + 1;
#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < samples_per_spit; i++)
				{
					if (m_progress.should_cancel) {
						// can't break an openmp loop, so we have to spin
						std::this_thread::yield();
						continue;
					}
					
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
						auto &vert = m_meshcache.get_vertex(vdi);
						TriMesh::VertexHandle v(vert.vi);
						// record sampled vertices
						if (vdi == rootvdi) m_mparams.mesh->property(m_mparams.prop_sampled, v) = true;
						//const float w = rootvdi == vdi;
						//const float w = r < exclusion_radius * 0.5f;
						// approx inverse distance, preventing inf
						// weight falls to some constant at exclusion radius, and zero at 2x that
						// zero is important to avoid discontinuity at edge of distribution radius
						// weight scale is not important
						const float re = exclusion_radius;
						const float b = 0.001f;
						const float w = std::max(re / (r + re * b) - (1 / (b + 2)), 0.f);
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
					subsampleGeodesicNeighborhoodSaliency(m_meshcache, stats, rootvdi, currentRadius, m_hMin, m_hMax, normalmap_filter * m_real_noise_height, exclusion_radius * 2.1f, visitor);

					stats.nh_timer_end();

					// TODO not actually completion, just num samples
					const int completion1 = completion.fetch_add(1);

					// update progress to some extent to show activity
					if (m_uparams.show_progress && omp_get_thread_num() == 0 && (stats.time0 - m_time_last_percent) > std::chrono::milliseconds(110)) {
						m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(stats.time0 - m_time_start);
						m_time_last_percent = stats.time0 - std::chrono::milliseconds(60);
						m_progress.levels[currentLevel].completed_samples = completion1;
						std::printf("\rSaliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, completion1);
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
				remaining_area = 0;
				for (auto it = candidates.begin(); it != candidates.end(); ) {
					if (it->deleted) {
						it = candidates.erase(it);
					} else {
						remaining_area += it->area;
						it++;
					}
				}

				// need to do progress outside inner loop to access remaining area
				if (m_uparams.show_progress && (m_thread_stats[0].time0 - m_time_last_percent) > std::chrono::milliseconds(50)) {
					m_time_last_percent = m_thread_stats[0].time0;
					const int completion1 = completion;
					m_progress.levels[currentLevel].completion = 1 - remaining_area;
					m_progress.levels[currentLevel].completed_samples = completion1;
					m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(m_thread_stats[0].time0 - m_time_start);
					std::printf("\rSaliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples, %7.3f%%", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, completion1, (1 - remaining_area) * 100);
				}

				if (candidates.empty()) {
					const int completion1 = completion;
					m_progress.levels[currentLevel].completion = 1;
					m_progress.levels[currentLevel].completed_samples = completion1;
					std::printf("\rSaliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples, %7.3f%% : no more sample candidates\n", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, completion1, 100.f);
					const float actual_subsampling = m_mparams.mesh->n_vertices() / float(completion1);
					std::cout << "Actual saliency subsampling: " << (100.f / actual_subsampling) << "% (" << actual_subsampling << "x)" << std::endl;
					break;
				}

			}

			for (auto vIt = m_mparams.mesh->vertices_begin(), vEnd = m_mparams.mesh->vertices_end(); vIt != vEnd; ++vIt) {
				const auto &p = tempSaliencyProperty[vIt->idx()];
				m_mparams.mesh->property(m_mparams.prop_saliency_levels[currentLevel], *vIt) = p.s / (p.w + 0.001f);
			}

		}

		~SaliencyComputation() {

		}
	};
	
	saliency_result compute_saliency(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress) {
		// note: saliency computation should not create/destroy properties,
		// only use them, to minimize problems (potential corruption) from concurrent access
		const int prev_max_threads = omp_get_max_threads();
		if (uparams.thread_count) omp_set_num_threads(uparams.thread_count);
		SaliencyComputation s(mparams, uparams, progress);
		bool r = s.run() && !progress.should_cancel;
		if (!r) progress.state = saliency_computation_state::cancelled;
		omp_set_num_threads(prev_max_threads);
		return saliency_result(mparams.cleanup, r);
	}

}
