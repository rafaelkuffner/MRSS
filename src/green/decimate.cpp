
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#include <iostream>
#include <vector>
#include <algorithm>

#include "OpenMesh/Tools/Decimater/DecimaterT.hh"
#include "OpenMesh/Tools/Decimater/ModQuadricT.hh"
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>

#include "decimate.hpp"

namespace {

	using namespace green;

	using sal_prop_t = decltype(decimate_mesh_params::prop_saliency);

	// left bin edges (with trailing upper edge); weights must be normalized
	std::vector<float> saliency_bin_edges(const TriMesh &mesh, sal_prop_t prop, const std::vector<float> &bin_weights, float target_ratio) {
		std::vector<float> saliencies = mesh.property(prop).data_vector();
		std::sort(saliencies.begin(), saliencies.end());
		int verts = 0;
		std::vector<float> sal_bin_edges(bin_weights.size() + 1, 0);
		for (int i = 1; i < bin_weights.size(); i++) {
			// bin sizes have to be determined by decimate-vs-keep weighting
			float w = (1.f - target_ratio) * bin_weights[i - 1] + target_ratio / bin_weights.size();
			verts += mesh.n_vertices() * w;
			sal_bin_edges[i] = saliencies[verts];
			std::cout << "saliency left edge for bin " << i << ": " << sal_bin_edges[i] << std::endl;
		}
		sal_bin_edges.back() = 1;
		return sal_bin_edges;
	}

	std::vector<int> saliency_bin_counts(const TriMesh &mesh, sal_prop_t prop, const std::vector<float> &bin_edges) {
		std::vector<int> counts(bin_edges.size() - 1, 0);
		for (auto &s : mesh.property(prop).data_vector()) {
			int b = 0;
			while (b < counts.size() - 1 && bin_edges[b + 1] <= s) b++;
			counts[b]++;
		}
		return counts;
	}

	template <class MeshT>
	class ModVertexWeightingT : public OpenMesh::Decimater::ModBaseT<MeshT>
	{
	public:

		// Defines the types Self, Handle, Base, Mesh, and CollapseInfo
		// and the memberfunction name()
		DECIMATING_MODULE(ModVertexWeightingT, MeshT, VertexWeighting)

	public:

		/** Constructor
		*  \internal
		*/
		ModVertexWeightingT(MeshT& _mesh)
			: Base(_mesh, true)
		{
		}


		/// Destructor
		virtual ~ModVertexWeightingT() {}


	public: // inherited
		virtual void initialize() {
			if (!prop_sal.is_valid()) throw std::runtime_error("no saliency");
		}

		virtual float collapse_priority(const CollapseInfo& _ci) {
			if (cancel) return Base::ILLEGAL_COLLAPSE;
			float s0 = _ci.mesh.property(prop_sal, _ci.v0);
			float s1 = _ci.mesh.property(prop_sal, _ci.v1);
			// supposedly v0 gets removed
			float s = s0; //(s0 + s1) * 0.5f;
			return min_sal <= s && s <= max_sal ? Base::LEGAL_COLLAPSE : Base::ILLEGAL_COLLAPSE;
		}

		virtual void postprocess_collapse(const CollapseInfo& _ci) {
			collapses++;
			if (collapses - collapses_last_progress > 1000) {
				collapses_last_progress = collapses;
				progress->completed_collapses = collapses;
				progress->elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
				if (progress->should_cancel) cancel = true;
			}
		}


	public: // local
		sal_prop_t prop_sal;
		float min_sal = 0;
		float max_sal = 1;
		decimate_progress *progress = nullptr;
		int collapses = 0;
		int collapses_last_progress = 0;
		std::chrono::steady_clock::time_point time_start;
		bool cancel = false;

	private:

	};
}

namespace green {
	
	void decimate_user_params::sanitize() {
		targetverts = std::max(targetverts, 0);
		nbins = std::max(nbins, 1);
		weight_falloff = std::max(weight_falloff, 0.001f);
	}

	bool decimate(const decimate_mesh_params &mparams, const decimate_user_params &uparams, decimate_progress &progress) {
		progress.state = decimation_state::bins;
		const auto time_start = std::chrono::steady_clock::now();
		std::cout << "decimating to " << uparams.targetverts << " vertices" << std::endl;

		if (uparams.targetverts >= mparams.mesh->n_vertices()) return true;

		const int nbins = uparams.nbins;

		const int target_collapses = mparams.mesh->n_vertices() - uparams.targetverts;
		progress.target_collapses = target_collapses;
		const float target_ratio = float(uparams.targetverts) / mparams.mesh->n_vertices();
		std::cout << "total vertices to remove: " << target_collapses << "; " << (1.f - target_ratio) << std::endl;

		std::vector<float> bin_weights(nbins, 1);
		std::vector<int> bin_target_collapses(nbins + 1, 0);
		float bin_weight_divisor = 0;
		for (int i = 0; i < nbins; i++) {
			// TODO better weight curve?
			float w = std::pow(uparams.weight_falloff, float(i) / nbins);
			bin_weights[i] = w;
			bin_weight_divisor += w;
		}
		for (int i = 0; i < nbins; i++) {
			float w = bin_weights[i] /= bin_weight_divisor;
			int t = target_collapses * w;
			bin_target_collapses[i] = t;
		}

		const auto sal_bin_edges = saliency_bin_edges(*mparams.mesh, mparams.prop_saliency, bin_weights, target_ratio);
		const auto init_bin_counts = saliency_bin_counts(*mparams.mesh, mparams.prop_saliency, sal_bin_edges);

		for (int i = 0; i < nbins; i++) {
			int t = bin_target_collapses[i];
			int c = init_bin_counts[i];
			std::cout << "vertices to remove from bin " << i << ": " << t << "/" << c << "; keep " << (c - t) << "; weight=" << bin_weights[i] << std::endl;
		}

		progress.state = decimation_state::run;
		OpenMesh::Decimater::DecimaterT<TriMesh> decimater(*mparams.mesh);
		
		OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle hModQuadrics;
		decimater.add(hModQuadrics);
		decimater.module(hModQuadrics).unset_max_err();

		ModVertexWeightingT<TriMesh>::Handle hModWeighting;
		decimater.add(hModWeighting);
		decimater.module(hModWeighting).progress = &progress;
		decimater.module(hModWeighting).prop_sal = mparams.prop_saliency;
		decimater.module(hModWeighting).time_start = time_start;

		for (int i = 0; i < nbins; i++) {
			if (progress.should_cancel) {
				progress.state = decimation_state::cancelled;
				return false;
			}
			std::cout << "decimating bin " << i << std::endl;
			decimater.initialize();
			auto &mod = decimater.module(hModWeighting);
			mod.min_sal = sal_bin_edges[i];
			mod.max_sal = sal_bin_edges[i + 1];
			int target = std::min<int>(bin_target_collapses[i], (1.f - target_ratio * 0.1f) * init_bin_counts[i]);
			int collapses = decimater.decimate_to(mparams.mesh->n_vertices() - target);
			std::cout << "vertices removed: " << collapses << std::endl;
			int collapses_remaining = bin_target_collapses[i] - collapses;
			if (collapses_remaining > 0) bin_target_collapses[i + 1] += collapses_remaining;
			mparams.mesh->garbage_collection();
		}

		const auto fin_bin_counts = saliency_bin_counts(*mparams.mesh, mparams.prop_saliency, sal_bin_edges);

		for (int i = 0; i < nbins; i++) {
			int collapses = (init_bin_counts[i] - fin_bin_counts[i]);
			std::cout << "vertices removed from bin " << i << ": " << collapses << "; ratio=" << (float(collapses) / target_collapses) << std::endl;
		}

		std::cout << "final vertices: " << mparams.mesh->n_vertices() << std::endl;

		progress.completed_collapses = decimater.module(hModWeighting).collapses;
		progress.elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
		progress.state = progress.should_cancel ? decimation_state::cancelled : decimation_state::done;
		return !progress.should_cancel;
	}

}
