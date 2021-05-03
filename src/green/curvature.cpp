
/*
 * Modified from:
 *
 * Example Code for the Paper
 *
 * "Mesh Saliency via Local Curvature Entropy", Eurographics 2016
 *
 * by M. Limper, A. Kuijper and D. Fellner
 *
 */

#include "curvature.hpp"

#include <limits>
#include <cmath>
#include <float.h>
#include <map>
#include <array>
#include <algorithm>
#include <iostream>

#include <omp.h>

#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

#include "meshutils.hpp"
#include "neighborhood.hpp"

namespace {

	template <size_t NBins>
	float histogram_entropy(const std::array<float, NBins> &hist) {
		float atot = 0;
		for (auto &a : hist) atot += a;
		const float iatot = 1.f / atot;
		float s = 0;
		const float nilog2 = -1.f / log(2.f);
		for (auto &a : hist) {
			// limit of p*log(p) as p tends to 0 is 0
			// so we can safely discard anything very small
			const float p = a * iatot;
			if (!(p > 1e-6f)) continue;
			s += p * log(p) * nilog2;
		}
		return s;
	}

	template <size_t NBins>
	void print_histogram(const std::array<float, NBins> &hist) {
		using namespace std;
		float atot = 0;
		for (auto &a : hist) atot += a;
		const float iatot = 1.f / atot;
		for (int i = 0; i < NBins; i++) {
			cout << setw(3) << i << " | ";
			const string_view maxbar = "================================================================================";
			int c = maxbar.size() * hist[i] * iatot;
			cout << maxbar.substr(0, c) << endl;
		}
	}

	// by area, normalized
	template <size_t NBins = 256, typename MonotonicTransform>
	inline std::array<float, NBins> histogram(
		const green::PolyMesh &mesh,
		OpenMesh::VPropHandleT<float> vertex_area_prop,
		OpenMesh::VPropHandleT<float> value_prop,
		MonotonicTransform &&func,
		float minval,
		float maxval,
		const std::vector<int> &sample_indices,
		size_t nsamples
	) {
		using namespace green;
		const float fmin = func(minval);
		const float hist_irange = 1.f / (func(maxval) - fmin);
		std::array<float, NBins> hist{};
		float atot = 0;
		for (size_t i = 0; i < std::min(sample_indices.size(), nsamples); i++) {
			PolyMesh::VertexHandle v(sample_indices[i]);
			const float area = mesh.property(vertex_area_prop, v);
			const float val = mesh.property(value_prop, v);
			const auto bin = intptr_t(float(NBins - 1) * std::clamp(hist_irange * (func(val) - fmin), 0.f, 1.f));
			hist[bin] += area;
			atot += area;
		}
		const float iatot = 1.f / atot;
		for (auto &a : hist) {
			a *= iatot;
		}
		return hist;
	}

}

namespace green {

	inline void update_minmax(curvature_measure &curv, float v) {
		curv.curv_min = std::min(curv.curv_min, v);
		curv.curv_max = std::max(curv.curv_max, v);
	}

	float doubleFaceArea(const OpenMesh::Vec3f & e0, const OpenMesh::Vec3f & e2)
	{
		return OpenMesh::cross(e0, -e2).norm();
	}

	void compute_mean_gaussian_curvature(PolyMesh &mesh, curvature_measure &curv_mean, curvature_measure &curv_gauss) {
		assert(mesh.has_vertex_normals());
		assert(curv_mean.prop_curv.is_valid());
		assert(curv_gauss.prop_curv.is_valid());
		curv_mean.init_ranges(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max());
		curv_gauss.init_ranges(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max());

		PolyMesh::VertexIter vIt, vEnd;

		//initialize properties and helpers
		std::vector<float> areaWeight(mesh.n_vertices());
		std::vector<OpenMesh::Vec3f> accDiffVec(mesh.n_vertices());

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			mesh.property(curv_gauss.prop_curv, *vIt) = M_PI * 2.0f;
			mesh.property(curv_mean.prop_curv, *vIt)     = 0.0f;

			areaWeight[vIt->idx()] = 0.0f;
			accDiffVec[vIt->idx()] = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		}

		//compute vertex weights via face areas

		PolyMesh::FaceIter fIt, fEnd;
		PolyMesh::FaceVertexIter vfIt;

		OpenMesh::Vec3f p0, p1, p2, e0, e1, e2;

		PolyMesh::VertexHandle v0, v1, v2;

		float l0, l1, l2;             //squared edge lengths
		float alpha0, alpha1, alpha2; //angles
		float vA0, vA1, vA2;          //areas

		float oneByTanAlpha0, oneByTanAlpha1, oneByTanAlpha2;
		float obtA;

		float k;

		const float PiHalf   = M_PI * 0.5f;
		const float AngleEps = 1e-15;

		for (fIt = mesh.faces_sbegin(), fEnd = mesh.faces_end(); fIt != fEnd; ++fIt)
		{
			vfIt = mesh.fv_iter(*fIt);

			v0 = *vfIt++;
			v1 = *vfIt++;
			v2 = *vfIt;

			bool nonmanifold = mesh.has_face_status() && mesh.status(*fIt).fixed_nonmanifold();
			if (nonmanifold) {
				int fid = fIt->idx();
				continue;
			}
			p0 = mesh.point(v0);
			p1 = mesh.point(v1);
			p2 = mesh.point(v2);

			e0 = p1 - p0;
			e1 = p2 - p1;
			e2 = p0 - p2;

			alpha0 = angle( e0, -e2);
			alpha1 = angle( e1, -e0);
			alpha2 = M_PI - (alpha0 + alpha1);

			if (alpha0 < AngleEps || alpha1 < AngleEps || alpha2 < AngleEps)
			{
				continue;
			}

			//accumulate vertex weights

			//obtuse cases
			if (alpha0 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA;
				areaWeight[v1.idx()] += obtA * 0.5f;
				areaWeight[v2.idx()] += obtA * 0.5f;
			}
			else if (alpha1 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA * 0.5f;
				areaWeight[v1.idx()] += obtA;
				areaWeight[v2.idx()] += obtA * 0.5f;
			}
			else if (alpha2 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA * 0.5f;
				areaWeight[v1.idx()] += obtA * 0.5f;
				areaWeight[v2.idx()] += obtA;
			}
			//non-obtuse case
			else
			{
				l0 = e0.sqrnorm();
				l1 = e1.sqrnorm();
				l2 = e2.sqrnorm();


				oneByTanAlpha0 = 1.0f / tan(alpha0);
				oneByTanAlpha1 = 1.0f / tan(alpha1);
				oneByTanAlpha2 = 1.0f / tan(alpha2);

				vA0 = (l2*oneByTanAlpha1 + l0*oneByTanAlpha2) / 8.0f;
				vA1 = (l0*oneByTanAlpha2 + l1*oneByTanAlpha0) / 8.0f;
				vA2 = (l1*oneByTanAlpha0 + l2*oneByTanAlpha1) / 8.0f;

				areaWeight[v0.idx()] += vA0;
				areaWeight[v1.idx()] += vA1;
				areaWeight[v2.idx()] += vA2;
			}

			//accumulate difference vectors
			if(alpha0 != 0.0f && alpha1 != 0.0f && alpha2 != 0.0f)
			{
				oneByTanAlpha0 = 1.0f / tan(alpha0);
				oneByTanAlpha1 = 1.0f / tan(alpha1);
				oneByTanAlpha2 = 1.0f / tan(alpha2);

				accDiffVec[v0.idx()] += (e2*oneByTanAlpha1 - e0*oneByTanAlpha2) * 0.25f;
				accDiffVec[v1.idx()] += (e0*oneByTanAlpha2 - e1*oneByTanAlpha0) * 0.25f;
				accDiffVec[v2.idx()] += (e1*oneByTanAlpha0 - e2*oneByTanAlpha1) * 0.25f;


				//accumulate values for gaussian curvature
				mesh.property(curv_gauss.prop_curv, v0) -= alpha0;
				mesh.property(curv_gauss.prop_curv, v1) -= alpha1;
				mesh.property(curv_gauss.prop_curv, v2) -= alpha2;

				//TODO: add case for adjusting gaussian curvature at borders
			}
		}


		//finally, compute final curvature values for vertices
		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			float a = 0;
			if (areaWeight[vIt->idx()] <= std::numeric_limits<float>::epsilon())
			{
				mesh.property(curv_gauss.prop_curv, *vIt) = 0.0f;
				mesh.property(curv_mean.prop_curv,  *vIt) = 0.0f;
				update_minmax(curv_gauss, 0);
				update_minmax(curv_mean, 0);
			}
			else
			{
				float g = mesh.property(curv_gauss.prop_curv, *vIt) /= areaWeight[vIt->idx()];
				update_minmax(curv_gauss, g);

				OpenMesh::Vec3f a = accDiffVec[vIt->idx()];
				float b = areaWeight[vIt->idx()];
				float k1 = (a / b).norm();
				k = (OpenMesh::dot(accDiffVec[vIt->idx()], mesh.normal(*vIt)) <= 0.0f ? -1.0f : 1.0f) * k1;

				mesh.property(curv_mean.prop_curv, *vIt) = k;
				update_minmax(curv_mean, k);
			}
		}
	}

	void compute_mean_curvature(PolyMesh &mesh, float scale, curvature_measure &curv) {
		assert(mesh.has_vertex_normals());
		assert(curv.prop_curv.is_valid());
		curv.init_ranges(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max());

		PolyMesh::VertexIter vIt, vEnd;

		//initialize properties and helpers
		std::vector<float> areaWeight(mesh.n_vertices());
		std::vector<OpenMesh::Vec3f> accDiffVec(mesh.n_vertices());

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			mesh.property(curv.prop_curv, *vIt) = 0.0f;

			areaWeight[vIt->idx()] = 0.0f;
			accDiffVec[vIt->idx()] = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		}

		//compute vertex weights via face areas

		PolyMesh::FaceIter fIt, fEnd;
		PolyMesh::FaceVertexIter vfIt;

		OpenMesh::Vec3f p0, p1, p2,
			e0, e1, e2;

		PolyMesh::VertexHandle v0, v1, v2;

		float l0, l1, l2;             //squared edge lengths
		float alpha0, alpha1, alpha2; //angles
		float vA0, vA1, vA2;          //areas

		float oneByTanAlpha0, oneByTanAlpha1, oneByTanAlpha2;
		float obtA;

		float k;

		const float PiHalf = M_PI * 0.5f;
		const float AngleEps = 1e-15;

		for (fIt = mesh.faces_sbegin(), fEnd = mesh.faces_end(); fIt != fEnd; ++fIt)
		{
			vfIt = mesh.fv_iter(*fIt);

			v0 = *vfIt++;
			v1 = *vfIt++;
			v2 = *vfIt;

			bool nonmanifold = mesh.has_face_status() && mesh.status(*fIt).fixed_nonmanifold();
			if (nonmanifold) {
				int fid = fIt->idx();
				continue;
			}
			p0 = scale * mesh.point(v0);
			p1 = scale * mesh.point(v1);
			p2 = scale * mesh.point(v2);

			e0 = p1 - p0;
			e1 = p2 - p1;
			e2 = p0 - p2;

			alpha0 = angle(e0, -e2);
			alpha1 = angle(e1, -e0);
			alpha2 = M_PI - (alpha0 + alpha1);

			if (alpha0 < AngleEps || alpha1 < AngleEps || alpha2 < AngleEps)
			{
				continue;
			}

			//accumulate vertex weights

			//obtuse cases
			if (alpha0 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA;
				areaWeight[v1.idx()] += obtA * 0.5f;
				areaWeight[v2.idx()] += obtA * 0.5f;
			}
			else if (alpha1 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA * 0.5f;
				areaWeight[v1.idx()] += obtA;
				areaWeight[v2.idx()] += obtA * 0.5f;
			}
			else if (alpha2 >= PiHalf)
			{
				obtA = doubleFaceArea(e0, e2) * 0.25f;

				areaWeight[v0.idx()] += obtA * 0.5f;
				areaWeight[v1.idx()] += obtA * 0.5f;
				areaWeight[v2.idx()] += obtA;
			}
			//non-obtuse case
			else
			{
				l0 = e0.sqrnorm();
				l1 = e1.sqrnorm();
				l2 = e2.sqrnorm();


				oneByTanAlpha0 = 1.0f / tan(alpha0);
				oneByTanAlpha1 = 1.0f / tan(alpha1);
				oneByTanAlpha2 = 1.0f / tan(alpha2);

				vA0 = (l2*oneByTanAlpha1 + l0 * oneByTanAlpha2) / 8.0f;
				vA1 = (l0*oneByTanAlpha2 + l1 * oneByTanAlpha0) / 8.0f;
				vA2 = (l1*oneByTanAlpha0 + l2 * oneByTanAlpha1) / 8.0f;

				areaWeight[v0.idx()] += vA0;
				areaWeight[v1.idx()] += vA1;
				areaWeight[v2.idx()] += vA2;
			}

			//accumulate difference vectors
			if (alpha0 != 0.0f && alpha1 != 0.0f && alpha2 != 0.0f)
			{
				oneByTanAlpha0 = 1.0f / tan(alpha0);
				oneByTanAlpha1 = 1.0f / tan(alpha1);
				oneByTanAlpha2 = 1.0f / tan(alpha2);

				accDiffVec[v0.idx()] += (e2*oneByTanAlpha1 - e0 * oneByTanAlpha2) * 0.25f;
				accDiffVec[v1.idx()] += (e0*oneByTanAlpha2 - e1 * oneByTanAlpha0) * 0.25f;
				accDiffVec[v2.idx()] += (e1*oneByTanAlpha0 - e2 * oneByTanAlpha1) * 0.25f;

			}
		}

		//finally, compute final curvature values for vertices
		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			float a = 0;
			if (areaWeight[vIt->idx()] <= 1e-30f)//std::numeric_limits<float>::epsilon())
			{
				mesh.property(curv.prop_curv, *vIt) = 0.0f;
				update_minmax(curv, 0);
			}
			else
			{
				OpenMesh::Vec3f a = accDiffVec[vIt->idx()];
				float b = areaWeight[vIt->idx()];
				float k1 = (a / b).norm();
				k = (OpenMesh::dot(accDiffVec[vIt->idx()], mesh.normal(*vIt)) <= 0.0f ? -1.0f : 1.0f) * k1;

				mesh.property(curv.prop_curv, *vIt) = k;
				update_minmax(curv, k);
			}
		}

		curv.natural_min = std::numeric_limits<float>::lowest();
		curv.natural_max = std::numeric_limits<float>::max();
	}

	float maxdon(const PolyMesh& mesh, PolyMesh::VertexHandle v, float normalPower) {
		float maxDiff = 1.0f;
		for (PolyMesh::ConstVertexOHalfedgeIter vohIt = mesh.cvoh_iter(v); vohIt.is_valid(); ++vohIt) {
			auto nv = mesh.to_vertex_handle(*vohIt);
			auto norm1 = mesh.normal(nv);
			for (PolyMesh::ConstVertexOHalfedgeIter vohIt2 = mesh.cvoh_iter(v); vohIt2.is_valid(); ++vohIt2) {
				auto nv2 = mesh.to_vertex_handle(*vohIt2);
				auto norm2 = mesh.normal(nv2);
				maxDiff = std::min(OpenMesh::dot(norm1, norm2), maxDiff);
			}
		}
		maxDiff = std::pow(((1 - (maxDiff)) / 2), normalPower); //Divided by two to get 0-1 range so our power parameter works :)
		maxDiff = std::clamp(maxDiff, 0.f, 1.f);
		if (!std::isfinite(maxDiff)) maxDiff = 0;
		return maxDiff;
	}

	void compute_maxdon(PolyMesh &mesh, float normal_power, curvature_measure &curv) {
		assert(mesh.has_vertex_normals());
		assert(curv.prop_curv.is_valid());
		curv.init_ranges(0, 1);
#pragma omp parallel for schedule(dynamic, 10000)
		for (int i = 0; i < mesh.n_vertices(); i++) {
			PolyMesh::VertexHandle vh(i);
			float d = maxdon(mesh, vh, normal_power);
			mesh.property(curv.prop_curv, vh) = d;
		}
		std::tie(curv.curv_min, curv.curv_max) = property_range(mesh, curv.prop_curv);
	}

	void compute_don(PolyMesh &mesh, OpenMesh::VPropHandleT<float> prop_vert_area, float normal_power, curvature_measure &curv) {
		assert(mesh.has_vertex_normals());
		assert(curv.prop_curv.is_valid());
		curv.init_ranges(0, 1);
		PolyMesh::VertexIter vIt, vEnd;
		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt) {
			float k = 0.0f;
			float w = 0;
			auto norm0 = mesh.normal(*vIt);
			for (auto vohIt = mesh.cvoh_iter(PolyMesh::VertexHandle(vIt)); vohIt.is_valid(); ++vohIt) {
				auto nv = mesh.to_vertex_handle(*vohIt);
				float area = mesh.property(prop_vert_area, nv);
				w += area;
				auto norm1 = mesh.normal(nv);
				k += area * std::max(OpenMesh::dot(norm0, norm1), 0.f);
			}
			k = std::pow((1 - (k / w)), normal_power);
			if (!std::isfinite(k)) k = 0;
			mesh.property(curv.prop_curv, *vIt) = k;
			update_minmax(curv, k);
		}
		curv.natural_min = 0;
		curv.natural_max = 1;
	}

	curvature_autocontrast autocontrast(
		const PolyMesh &mesh,
		const curvature_measure &curv,
		OpenMesh::VPropHandleT<float> prop_vert_area,
		bool natural_binning
	){
		float bin_min = natural_binning ? curv.natural_min : curv.curv_min;
		float bin_max = natural_binning ? curv.natural_max : curv.curv_max;
		curvature_autocontrast curvac{curv};
		std::minstd_rand rand{std::random_device{}()};
		std::vector<int> samples(mesh.n_vertices());
		std::iota(samples.begin(), samples.end(), 0);
		// TODO using all vertices for now
		const size_t nsamples = 100000000; //samples.size() / 2;
		if (samples.size() > nsamples) {
			// sort first part descending by area
			std::partial_sort(samples.begin(), samples.begin() + nsamples, samples.end(), [&](int ia, int ib) {
				float aa = mesh.property(prop_vert_area, PolyMesh::VertexHandle(ia));
				float ab = mesh.property(prop_vert_area, PolyMesh::VertexHandle(ib));
				return aa > ab;
			});
			// randomize second part
			//std::shuffle(samples.begin() + nsamples / 3, samples.end(), rand);
		}
		auto eval_entropy = [&](float contrast) {
			auto f = [=](float x) { return curvature_contrast(x, contrast); };
			auto hist = histogram<MeshCache::ncurvbins>(mesh, prop_vert_area, curv.prop_curv, f, bin_min, bin_max, samples, nsamples);
			float s = histogram_entropy(hist);
			return s;
		};
		auto eval_entropy_gradient = [&](float contrast, float delta) {
			float s0 = eval_entropy(contrast);
			float s1 = eval_entropy(contrast + delta);
			return std::pair{s0, (s1 - s0) / delta};
		};
		//for (float contrast = 0.1f; contrast < 1.6f; contrast += 0.1f) {
		//	float s = eval_entropy(contrast);
		//	std::cout << "Global entropy at contrast " << contrast << ": " << s << std::endl;
		//}
		const float target_entropy = autocontrast_target_entropy;
		float best_contrast = 1;
		float contrast_delta = 0.1f;
		for (int i = 0; i < 5; i++) {
			auto [s, dsdc] = eval_entropy_gradient(best_contrast, contrast_delta);
			std::cout << "Global entropy step: " << s << std::endl;
			best_contrast = best_contrast - (s - target_entropy) / dsdc;
			contrast_delta *= 0.7f;
		}
		std::cout << "Auto contrast: " << best_contrast << std::endl;
		if (!std::isfinite(best_contrast) || best_contrast < 0.01f) {
			std::cout << "Auto contrast failed, defaulting to 1.0" << std::endl;
			best_contrast = 1;
		}
		curvac.contrast = best_contrast;
		curvac.bin_min = curvature_contrast(bin_min, best_contrast);
		curvac.bin_max = curvature_contrast(bin_max, best_contrast);
		return curvac;
	}

}
