
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
	}
	
	std::future<saliency_result> compute_saliency_async(const saliency_mesh_params &mparams, const saliency_user_params &uparams, saliency_progress &progress) {
		return std::async([=, &progress]() { return compute_saliency(mparams, uparams, progress); });
	}

	class SaliencyComputation {
	private:
		struct initial_sample_candidate {
			unsigned vdi = -1;
			float area = 0;
			float summedarea = 0;
			//bool deleted = false;
		};

		struct sample_candidate : initial_sample_candidate {
			std::atomic<bool> deleted2{false};

			sample_candidate(const initial_sample_candidate &other) : initial_sample_candidate{other} {}
		};

		saliency_mesh_params m_mparams;
		saliency_user_params m_uparams;
		saliency_progress &m_progress;
		MeshCache m_meshcache;
		std::vector<initial_sample_candidate> m_candidates0;
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
			std::cout << "Contrast: " << m_uparams.normal_power << std::endl;
			
			std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;

			// TODO curv selection param
			// note: we can't cache the (final) curvature because we couldn't adjust the normal power etc
			m_progress.state = saliency_state::curv;
			std::cout << "Computing curvature" << std::endl;
			float curvmax = -9001e19;
			float curvmin = 9001e19;
			for (auto &v : m_mparams.mesh->vertices()) {
				float c = m_mparams.mesh->property(m_mparams.prop_doncurv_raw, v);
				c = std::pow(c, m_uparams.normal_power);
				curvmax = std::max(curvmax, c);
				curvmin = std::min(curvmin, c);
				m_mparams.mesh->property(m_mparams.prop_curvature, v) = c;
			}

			//m_hMin = curvhist.getMin();
			//m_hMax = curvhist.getMax();
			// force 0-1 histogram range (for entropy calc)
			// lower bound is pretty much always very near zero
			// upper bound is often very near 1 (kinda weird - thats adjacent vertices with anti-parallel normals)
			// using constant histogram range makes things more predictable
			// note: range 0-1 currently hardcoded in MeshCache
			m_hMin = 0;
			m_hMax = 1;
			std::cout << "Curv min: " << curvmin << std::endl;
			std::cout << "Curv max: " << curvmax << std::endl;

			m_progress.state = saliency_state::area;
			std::cout << "Computing surface area" << std::endl;
			m_surfaceArea = surfaceArea(*m_mparams.mesh);
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(std::chrono::steady_clock::now() - m_time_start);
			std::cout << "Surface area: " << m_surfaceArea << std::endl;

			m_real_noise_height = m_uparams.noise_height * sqrt(m_surfaceArea);
			std::cout << "Real noise height: " << m_real_noise_height << std::endl;

			m_progress.state = saliency_state::nhprep;
			const auto time_nhprep_start = std::chrono::steady_clock::now();
			m_meshcache = MeshCache(*m_mparams.mesh, m_mparams.prop_edge_length, m_mparams.prop_vertex_area, m_mparams.prop_curvature);
			const auto time_nhprep_finish = std::chrono::steady_clock::now();
			m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(time_nhprep_finish - m_time_start);
			std::cout << "Neighborhood search prep took " << ((time_nhprep_finish - time_nhprep_start) / std::chrono::duration<double>(1.0)) << "s" << std::endl;
			//m_meshcache.dump_to_file();

			m_progress.state = saliency_state::cand;
			std::cout << "Preparing subsampling candidates" << std::endl;
			m_candidates0.reserve(m_meshcache.vdis.size());
			std::transform(m_meshcache.vdis.begin(), m_meshcache.vdis.end(), std::back_inserter(m_candidates0), [](auto &vdi) { return initial_sample_candidate{vdi}; });

			// randomize sample candidates
			// makes subsampling more stable and less sensitive to parallelization
			std::minstd_rand rand{std::random_device{}()};
			std::shuffle(m_candidates0.begin(), m_candidates0.end(), rand);

			{
				float aa = 0;
				for (auto &cand : m_candidates0) {
					const float a = m_meshcache.decode_area(m_meshcache.get_vertex(cand.vdi).areax);
					cand.area = a;
					cand.summedarea = aa;
					aa += a;
				}
				m_total_vertex_area = aa;
				for (auto &cand : m_candidates0) {
					PolyMesh::VertexHandle v(m_meshcache.get_vertex_aux(cand.vdi).vi);
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

			CalculationStats overall_stats;

			// compute saliency at multiple levels
			for (int currentLevel = 0; currentLevel < m_uparams.levels; currentLevel++)
			{
				// thread stats per level
				m_thread_stats.assign(omp_get_max_threads(), {});
				
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
					m_progress.state = saliency_state::run_sub;
					run_level_subsampled(currentLevel, currentRadius, subsampling, normalmap_filter);
				} else {
					m_progress.levels[currentLevel].subsampled = false;
					m_progress.state = saliency_state::run_full;
					run_level_full(currentLevel, currentRadius, normalmap_filter);
				}

				// merge and dump stats for this level
				CalculationStats level_stats;
				for (auto &ts : m_thread_stats) {
					level_stats.merge(ts);
				}
				std::cout << "Level neighborhood stats:" << std::endl;
				level_stats.dump_stats(omp_get_max_threads());
				// then merge that to the overall stats
				overall_stats.merge(level_stats);

				m_progress.completed_levels = currentLevel + 1;
				if (m_progress.should_cancel) return false;
			}

			std::cout << "Overall neighborhood stats:" << std::endl;
			overall_stats.dump_stats(omp_get_max_threads());

			// merge saliency values to a single score / property
			m_progress.state = saliency_state::merge;
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
			m_progress.state = saliency_state::norm;
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
			std::cout << "Saliency time: " << (m_progress.elapsed_time / std::chrono::duration<double>(1.0)) << "s" << std::endl;
			m_progress.state = saliency_state::done;

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
				const auto sal = getGeodesicNeighborhoodSaliency(m_meshcache, stats, vdis, currentRadius, normalmap_filter * m_real_noise_height);

				stats.nh_timer_end();

				for (int j = 0; j < simd_traits::simd_size; j++) {
					if (vdis[j] == -1) break;

					PolyMesh::VertexHandle v(m_meshcache.get_vertex_aux(vdis[j]).vi);

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

			// candidates vector will not be resized, so iterators will remain valid
			using candidate_iter_t = std::vector<sample_candidate>::iterator;
			std::vector<sample_candidate> candidates(m_candidates0.begin(), m_candidates0.end());
			std::vector<candidate_iter_t> candidateProperty(m_mparams.mesh->n_vertices());

			// number of non-deleted candidates
			std::atomic<int> remaining_candidates{int(candidates.size())};

			// area of non-excluded vertices
			// can use this to report actual progress %
			std::atomic<float> remaining_area{1.f};

			// map vertices to candidates
			for (auto it = candidates.begin(); it != candidates.end(); ++it) {
				PolyMesh::VertexHandle v(m_meshcache.get_vertex_aux(it->vdi).vi);
				candidateProperty[v.idx()] = it;
			}

			struct next_candidate {
				candidate_iter_t it;
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
				const candidate_iter_t none = candidates.end();
				const auto it0 = nc.it;
				candidate_iter_t itb = none;
				float a1 = nc.summedarea + a;
				while (nc.it->summedarea < a1 || nc.it->deleted2.load(std::memory_order_relaxed)) {
					if (cand_inc_loop(nc.it)) a1 = std::max(0.f, a1 - 1.f);
					if (nc.it == it0) {
						// no more candidates
						// if the only candidate(s) are between it0 and a1 we won't find them 'properly'
						// so use any valid candidate we've seen so far
						if (itb == none) {
							// no valid candidates => exhausted
							return true;
						} else {
							nc.it = itb;
							nc.summedarea = a1;
							return false;
						}
					}
					if (itb == none && !nc.it->deleted2.load(std::memory_order_relaxed)) itb = nc.it;
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

			// note: the randomness depends on shuffling the candidates beforehand
			// currently on lucy (default params) about 3.5% execution time goes here, and about 6% at -a 0.01 -s 50
			auto random_vdi = [&](int threadindex) {
				const candidate_iter_t none = candidates.end();
				if (candidates.empty()) return unsigned(-1);
				auto &nc = thread_next_cand[threadindex];
				if (nc.it == none) {
					// candidates already exhausted
					return unsigned(-1);
				}
				if (!nc.it->deleted2.load(std::memory_order_relaxed)) {
					// candidate already usable
					return nc.it->vdi;
				}
				if (cand_inc_loop_area(nc, float(subsampling) / m_mparams.mesh->n_vertices())) {
					// search for next candidate exhausted (actually now)
					// note: the 'probably' used to prevent us from sentinel-izing the iterator
					nc.it = none;
					return unsigned(-1);
				}
				// note: should not use the deleted flag as a reason to return -1
				// deleted flag may be set by another thread after we selected a valid candidate
				return nc.it->vdi;
			};

			const auto atomic_float_accum = [](std::atomic<float> &a, float b) {
				float a0 = a.load(std::memory_order_relaxed);
				while (!a.compare_exchange_weak(a0, a0 + b, std::memory_order_relaxed));
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

			// radius within which to prevent future samples
			const float exclusion_radius = subsampling_radius * sampling_correction;

			std::atomic_flag progress_lock{false};

			auto threadproc = [&](int threadindex) {
				while (!m_progress.should_cancel) {
					
					auto &stats = m_thread_stats[threadindex];
					stats.timer_begin();

					// yes, this may select candidates that have been (or will be) flagged as deleted
					// we have to put up with this in order to parallelize
					const auto rootvdi = random_vdi(threadindex);
					if (rootvdi == -1) {
						// this should only happen when we're out of sample candidates (actually now)
						return;
					}

					const auto visitor = [&](unsigned vdi, float r, float s) {
						auto &vadata = m_meshcache.get_vertex_aux(vdi);
						PolyMesh::VertexHandle v(vadata.vi);
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
							// note: cannot actually erase anything
							auto it = candidateProperty[v.idx()];
							bool deleted0 = it->deleted2.load(std::memory_order_relaxed);
							if (!deleted0 && it->deleted2.compare_exchange_strong(deleted0, true, std::memory_order_relaxed)) {
								remaining_candidates.fetch_sub(1, std::memory_order_relaxed);
								atomic_float_accum(remaining_area, -it->area);
							}
						}
					};

					// could be up to exclusion_radius * 2 distance between samples, so distribute over that to ensure coverage
					// need to go slightly more to avoid problems due to vertex discretization (to allow weight to reach zero)
					subsampleGeodesicNeighborhoodSaliency(m_meshcache, stats, rootvdi, currentRadius, normalmap_filter * m_real_noise_height, exclusion_radius * 2.1f, visitor);

					stats.nh_timer_end();

					// TODO not actually completion, just num samples
					const int completion1 = completion.fetch_add(1, std::memory_order_relaxed);

					// update progress
					if (m_uparams.show_progress && !progress_lock.test_and_set(std::memory_order_acquire)) {
						if ((stats.time0 - m_time_last_percent) > std::chrono::milliseconds(110)) {
							m_progress.elapsed_time = std::chrono::duration_cast<decltype(m_progress.elapsed_time)>(stats.time0 - m_time_start);
							m_time_last_percent = stats.time0 - std::chrono::milliseconds(60);
							const float remarea1 = remaining_area.load(std::memory_order_relaxed);
							m_progress.levels[currentLevel].completion = 1 - remarea1;
							m_progress.levels[currentLevel].completed_samples = completion1;
							std::printf("\rSaliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples, %7.3f%%", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, completion1, (1 - remarea1) * 100);
						}
						progress_lock.clear(std::memory_order_release);
					}
				}
			};

			// just use omp as a threadpool (lol)
#pragma omp parallel for schedule(static)
			for (int threadindex = 0; threadindex < omp_get_max_threads(); threadindex++) {
				threadproc(threadindex);
			}

			{
				const int completion1 = completion;
				m_progress.levels[currentLevel].completion = 1;
				m_progress.levels[currentLevel].completed_samples = completion1;
				std::printf("\rSaliency computation [subsampled] @%8.2fs (lv %u/%d): %9d samples, %7.3f%% : finished\n", m_progress.elapsed_time.count() / 1000.0, currentLevel + 1, m_uparams.levels, completion1, 100.f);
				const float actual_subsampling = m_mparams.mesh->n_vertices() / float(completion1);
				std::cout << "Actual saliency subsampling: " << (100.f / actual_subsampling) << "% (" << actual_subsampling << "x)" << std::endl;
				if (remaining_candidates != 0) {
					std::cout << "WARNING: candidates not actually exhausted! remaining: " << remaining_candidates << std::endl;
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
		SaliencyComputation s(mparams, uparams, progress);
		bool r = s.run() && !progress.should_cancel;
		if (!r) progress.state = saliency_state::cancelled;
		return saliency_result(mparams.cleanup, r);
	}

}
