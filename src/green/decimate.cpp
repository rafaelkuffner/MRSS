
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>

#include "OpenMesh/Tools/Decimater/DecimaterT.hh"
#include "OpenMesh/Tools/Decimater/ModQuadricT.hh"
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>

#include "decimate.hpp"

namespace {

	using namespace green;
	using namespace OpenMesh;

	using sal_prop_t = decltype(decimate_mesh_params::prop_saliency);

	// left bin edges (with trailing upper edge)
	std::vector<float> saliency_bin_edges(PolyMesh &mesh, sal_prop_t prop, VPropHandleT<int> prop_bin, int nbins) {
		struct salvert {
			VertexHandle v;
			float s = 0;
			bool operator<(const salvert &rhs) {
				return s < rhs.s;
			}
		};
		std::vector<salvert> saliencies(mesh.n_vertices());
		for (int i = 0; i < saliencies.size(); i++) {
			VertexHandle v(i);
			auto &sv = saliencies[i];
			sv.v = v;
			sv.s = prop.is_valid() ? mesh.property(prop, v) : 0.f;
		}
		std::sort(saliencies.begin(), saliencies.end());
		// round up so we dont have unbinned vertices at the end
		const int verts_per_bin = 1 + (mesh.n_vertices() - 1) / nbins;
		std::vector<float> sal_bin_edges(nbins + 1, 0);
		int verts = 0;
		for (int i = 0; i < nbins; i++) {
			for (int j = verts; j < verts + verts_per_bin && j < saliencies.size(); j++) {
				mesh.property(prop_bin, saliencies[j].v) = i;
			}
			sal_bin_edges[i] = saliencies[verts].s;
			std::cout << "saliency left edge for bin " << i << ": " << sal_bin_edges[i] << std::endl;
			verts += verts_per_bin;
		}
		sal_bin_edges.back() = saliencies.back().s;
		return sal_bin_edges;
	}

	std::vector<int> saliency_bin_counts(const PolyMesh &mesh, VPropHandleT<int> prop_bin, int nbins) {
		std::vector<int> counts(nbins, 0);
		for (auto &bin : mesh.property(prop_bin).data_vector()) {
			counts[bin]++;
		}
		return counts;
	}

	template <typename MeshT>
	auto seam_props(const MeshT &m, VertexHandle vh, HalfedgeHandle hh) {
		return std::make_tuple(
			m.has_halfedge_normals() ? m.normal(hh) : Vec3f(0),
			m.has_halfedge_texcoords2D() ? m.texcoord2D(hh) : Vec2f(0),
			m.has_halfedge_texcoords3D() ? m.texcoord3D(hh) : Vec3f(0)
		);
	}

	template <typename MeshT>
	bool is_seam(const MeshT &m, HalfedgeHandle hh) {
		// also check connectivity boundary
		if (m.is_boundary(hh)) return true;
		VertexHandle vh = m.to_vertex_handle(hh);
		HalfedgeHandle hb = m.prev_halfedge_handle(m.opposite_halfedge_handle(hh));
		auto p0 = seam_props(m, vh, hh);
		auto pb = seam_props(m, vh, hb);
		return pb != p0;
	}
	
	template <typename MeshT>
	void compute_seams(MeshT &mesh, EPropHandleT<uint8_t> prop_seam, VPropHandleT<uint8_t> prop_seamarity) {
		// edge seam prop is effectively boolean (defended against std::vector<bool>)
		// label seam edges
		for (auto eh : mesh.edges()) {
			const bool b0 = is_seam(mesh, mesh.halfedge_handle(eh, 0));
			const bool b1 = is_seam(mesh, mesh.halfedge_handle(eh, 1));
			mesh.property(prop_seam, eh) = bool(b0 || b1);
		}
		for (auto vh : mesh.vertices()) {
			int boundarity = 0;
			// search incoming halfedges at each vertex for seam edges
			for (auto it = mesh.cvih_begin(vh); it.is_valid(); ++it) {
				boundarity += bool(mesh.property(prop_seam, mesh.edge_handle(*it)));
			}
			mesh.property(prop_seamarity, vh) = uint8_t(std::min<int>(255, boundarity));
		}
	}

	template <class MeshT>
	class ModVertexWeightingT : public Decimater::ModBaseT<MeshT>
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
			std::printf("decimating @%8.2fs : %9d collapses, %5.1f%%", progress->elapsed_time / std::chrono::duration<double>(1.0), collapses, 100.f * collapses / progress->target_collapses);
		}

		virtual float collapse_priority(const CollapseInfo& _ci) {
			// note: v0 gets removed
			auto &mesh = _ci.mesh;
			if (cancel) return Base::ILLEGAL_COLLAPSE;
			if (prop_seam.is_valid() && prop_seamarity.is_valid()) {
				const auto sa0 = mesh.property(prop_seamarity, _ci.v0);
				const auto sa1 = mesh.property(prop_seamarity, _ci.v1);
				// dont remove anything with boundarity 1 (terminus)
				if (sa0 == 1) return Base::ILLEGAL_COLLAPSE;
				// dont remove anything with boundarity > 2 (patch junction)
				if (sa0 > 2) return Base::ILLEGAL_COLLAPSE;
				// dont collapse edge into non-boundary
				if (sa0 == 2 && sa1 == 0) return Base::ILLEGAL_COLLAPSE;
				// 2 vertices on boundary edges may be connected by a non-boundary edge
				if (sa0 > 1 && sa1 > 1 && !mesh.property(prop_seam, mesh.edge_handle(_ci.v0v1))) return Base::ILLEGAL_COLLAPSE;
			}
			if (prop_bin.is_valid()) {
				const auto bin = mesh.property(prop_bin, _ci.v0);
				// only allow collapse within current bin
				if (bin != current_bin) return Base::ILLEGAL_COLLAPSE;
			}
			return Base::LEGAL_COLLAPSE;
		}

		virtual void preprocess_collapse(const CollapseInfo &_ci) override {
			fix_properties(_ci);
		}

		virtual void postprocess_collapse(const CollapseInfo &_ci) override {
			update_progress();
			// TODO turn off check when sure of correctness
			check_boundaries(_ci);
		}

		void check_boundaries(const CollapseInfo &_ci) {
			auto &mesh = _ci.mesh;
			for (auto it = mesh.cvih_begin(_ci.v1); it.is_valid(); ++it) {
				if (is_seam(mesh, *it) && !mesh.property(prop_seam, mesh.edge_handle(*it))) {
					std::cout << "BOUNDARY INCONSISTENCY" << std::endl;
				}
			}
		}

		// fix halfedge properties wrt seams so the result of collapsing is correct
		void fix_properties(const CollapseInfo &_ci) {
			if (!prop_seam.is_valid()) return;
			auto &mesh = _ci.mesh;
			const HalfedgeHandle h0 = _ci.v0v1;
			const HalfedgeHandle h1 = mesh.next_halfedge_handle(h0);
			const HalfedgeHandle h2 = mesh.next_halfedge_handle(h1);
			const HalfedgeHandle o0 = _ci.v1v0;
			const HalfedgeHandle o1 = mesh.next_halfedge_handle(o0);
			const HalfedgeHandle o2 = mesh.next_halfedge_handle(o1);

			// copy properties to other incoming halfedges at v0
			HalfedgeHandle hs = h0;
			for (HalfedgeHandle h = mesh.prev_halfedge_handle(h0); h != o0; h = mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(h))) {
				mesh.copy_all_properties(hs, h, true);
				// switch source halfedge when crossing seam
				// this relies on v0 having boundarity <= 2
				if (mesh.property(prop_seam, mesh.edge_handle(h))) hs = mesh.prev_halfedge_handle(o0);
			}

			// preserve seam flags of edges merged by loop collapse (of 2-gons)
			if (mesh.next_halfedge_handle(h2) == h0 && mesh.property(prop_seam, mesh.edge_handle(h2))) {
				mesh.property(prop_seam, mesh.edge_handle(h1)) = true;
			}
			if (mesh.next_halfedge_handle(o2) == o0 && mesh.property(prop_seam, mesh.edge_handle(o1))) {
				mesh.property(prop_seam, mesh.edge_handle(o2)) = true;
			}
		}

		void update_progress() {
			collapses++;
			if (collapses - collapses_last_progress > 5000) {
				collapses_last_progress = collapses;
				progress->completed_collapses = collapses;
				progress->elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
				if (progress->should_cancel) cancel = true;
				if (show_progress) std::printf("\rdecimating @%8.2fs : %9d collapses, %5.1f%%", progress->elapsed_time / std::chrono::duration<double>(1.0), collapses, 100.f * collapses / progress->target_collapses);
			}
		}

	public: // local
		VPropHandleT<int> prop_bin;
		VPropHandleT<uint8_t> prop_seamarity;
		EPropHandleT<uint8_t> prop_seam;
		int current_bin = 0;
		decimate_progress *progress = nullptr;
		int collapses = 0;
		int collapses_last_progress = 0;
		std::chrono::steady_clock::time_point time_start;
		bool cancel = false;
		bool show_progress = true;

	private:

	};
}

namespace green {
	
	void decimate_user_params::sanitize() {
		targetverts = std::max(targetverts, 0);
		targettris = std::max(targettris, 0);
		nbins = std::max(nbins, 1);
		weight = std::max(weight, 0.f);
		power = std::clamp(power, 0.1f, 10.f);
	}

	bool decimate(const decimate_mesh_params &mparams, const decimate_user_params &uparams, decimate_progress &progress) {
		progress.state = decimation_state::bins;
		const auto time_start = std::chrono::steady_clock::now();

		// NOTE needs triangulated mesh to work
		const float verts_per_tri = float(mparams.mesh->n_vertices()) / float(mparams.mesh->n_faces());
		const int targetverts = uparams.use_tris ? int(uparams.targettris * verts_per_tri) : uparams.targetverts;

		if (uparams.use_tris) std::cout << "decimating to approximately " << uparams.targettris << " triangles" << std::endl;
		std::cout << "decimating to " << targetverts << " vertices" << std::endl;

		if (targetverts >= mparams.mesh->n_vertices()) return true;

		// seam flags based on differences in halfedge properties
		// (normal, texcoord2d, texcoord3d)
		OpenMesh::EPropHandleT<uint8_t> prop_seam;
		mparams.mesh->add_property(prop_seam);
		// number of seam edges
		OpenMesh::VPropHandleT<uint8_t> prop_seamarity;
		mparams.mesh->add_property(prop_seamarity);
		compute_seams(*mparams.mesh, prop_seam, prop_seamarity);

		// explicit bin tags to deal with cases where there are large numbers
		// of vertices with the same saliency values preventing bin discrimination
		// TODO ensure this property is removed when eg cancelled
		// TODO smaller data type?
		OpenMesh::VPropHandleT<int> prop_bin;
		mparams.mesh->add_property(prop_bin);

		const int nbins = uparams.nbins;

		const int target_collapses = mparams.mesh->n_vertices() - targetverts;
		progress.target_collapses = target_collapses;
		const float target_ratio = float(targetverts) / mparams.mesh->n_vertices();
		std::cout << "total vertices to remove: " << target_collapses << "; " << (1.f - target_ratio) << std::endl;

		const int bin_min_keep = targetverts * 0.1f / float(nbins);

		std::vector<float> bin_weights(nbins, 0);
		std::vector<int> bin_keep(nbins, 0);
		float bin_weight_divisor = 0;
		for (int i = 0; i < nbins; i++) {
			// bin weight params: higher bin weight => more vertices _kept_ in that bin
			// (saliency) weight: 0 (even weights across bins) .. 1 (least weight to low bin, most weight to high bin)
			// power: non-linearity of weighting; 1 (linear) .. 2 (quadratic, more extreme)
			const float d = uparams.weight;
			const float x = (i + 0.5f) / float(nbins);
			const float w = std::max(pow(x, uparams.power) * d + (1.f - d) * 0.5f, 0.f);
			bin_weights[i] = w;
			bin_weight_divisor += w;
		}
		for (int i = 0; i < nbins; i++) {
			float w = bin_weights[i] /= bin_weight_divisor;
			bin_keep[i] = std::max<int>(targetverts * w, bin_min_keep);
		}

		const auto sal_bin_edges = saliency_bin_edges(*mparams.mesh, mparams.prop_saliency, prop_bin, nbins);
		const auto init_bin_counts = saliency_bin_counts(*mparams.mesh, prop_bin, nbins);

		for (int i = 0; i < nbins; i++) {
			int k = bin_keep[i];
			int c = init_bin_counts[i];
			std::cout << "bin " << i << ": keep " << k << "/" << c << " (" << bin_weights[i] << "); remove " << (c - k) << std::endl;
			if (k > c) {
				std::cout << "warning: can't keep more vertices than available, result will have fewer vertices than desired" << std::endl;
			}
		}

		progress.state = decimation_state::run;
		OpenMesh::Decimater::DecimaterT<PolyMesh> decimater(*mparams.mesh);
		
		OpenMesh::Decimater::ModQuadricT<PolyMesh>::Handle hModQuadrics;
		decimater.add(hModQuadrics);
		decimater.module(hModQuadrics).unset_max_err();

		ModVertexWeightingT<PolyMesh>::Handle hModWeighting;
		decimater.add(hModWeighting);
		decimater.module(hModWeighting).progress = &progress;
		decimater.module(hModWeighting).time_start = time_start;
		decimater.module(hModWeighting).show_progress = uparams.show_progress;
		decimater.module(hModWeighting).prop_seam = prop_seam;
		decimater.module(hModWeighting).prop_seamarity = prop_seamarity;
		
		if (uparams.use_saliency) {
			decimater.module(hModWeighting).prop_bin = prop_bin;
		} else {
			std::cout << "decimating without saliency, bins will be ignored" << std::endl;
		}

		for (int i = 0; i < nbins; i++) {
			if (progress.should_cancel) {
				progress.state = decimation_state::cancelled;
				return false;
			}
			std::cout << "decimating bin " << i << std::endl;
			decimater.initialize();
			auto &mod = decimater.module(hModWeighting);
			mod.current_bin = i;
			int target = init_bin_counts[i] - bin_keep[i];
			int collapses = decimater.decimate_to(std::max<size_t>(mparams.mesh->n_vertices() - target, targetverts));
			progress.completed_collapses = mod.collapses;
			std::printf("\rdecimating @%8.2fs : %9d collapses, %5.1f%%\n", progress.elapsed_time / std::chrono::duration<double>(1.0), progress.completed_collapses, 100.f * progress.completed_collapses / progress.target_collapses);
			std::cout << "vertices removed: " << collapses << std::endl;
			mparams.mesh->garbage_collection();
		}

		const auto fin_bin_counts = saliency_bin_counts(*mparams.mesh, prop_bin, nbins);

		mparams.mesh->remove_property(prop_bin);
		mparams.mesh->remove_property(prop_seam);
		mparams.mesh->remove_property(prop_seamarity);

		for (int i = 0; i < nbins; i++) {
			int k = fin_bin_counts[i];
			int r = init_bin_counts[i] - k;
			std::cout << "bin " << i << ": kept " << k << " (" << (float(k) / mparams.mesh->n_vertices()) << "); removed " << r << " (" << (float(r) / target_collapses) << ")" << std::endl;
		}

		std::cout << "final vertices: " << mparams.mesh->n_vertices() << std::endl;
		std::cout << "final triangles: " << mparams.mesh->n_faces() << std::endl;
		progress.completed_collapses = decimater.module(hModWeighting).collapses;
		progress.elapsed_time = std::chrono::duration_cast<decltype(decimate_progress::elapsed_time)>(std::chrono::steady_clock::now() - time_start);
		std::cout << "decimate time: " << (progress.elapsed_time / std::chrono::duration<double>(1.0)) << "s" << std::endl;
		progress.state = progress.should_cancel ? decimation_state::cancelled : decimation_state::done;
		return !progress.should_cancel;
	}

}
