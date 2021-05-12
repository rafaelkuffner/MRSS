
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

#pragma once

#ifndef GREEN_CURVATURE_HPP
#define GREEN_CURVATURE_HPP

#include "meshutils.hpp"

#include <cmath>
#include <random>
#include <algorithm>
#include <limits>
#include <utility>

namespace green {

	struct curvature_measure {
		OpenMesh::VPropHandleT<float> prop_curv;
		// natural range of values for this type of curvature (eg [0,1] for don)
		float natural_min = 0;
		float natural_max = 0;
		// actual range of values (_without_ contrast)
		float curv_min = 0;
		float curv_max = 0;

		curvature_measure() = default;

		explicit curvature_measure(PolyMesh &mesh) {
			mesh.add_property(prop_curv);
		}

		void init_ranges(float natmin, float natmax) {
			natural_min = natmin;
			natural_max = natmax;
			// set suitable for updating with min/max
			curv_min = std::numeric_limits<float>::max();
			curv_max = std::numeric_limits<float>::lowest();
		}
	};

	inline float curvature_contrast(float x, float e) {
		// ensure we can use negative inputs sensibly
		return powsign(x, e);
	}

	struct curvature_autocontrast : curvature_measure {
		float contrast = 1;
		bool natural_binning = false;
	};

	void compute_mean_gaussian_curvature(PolyMesh &mesh, curvature_measure &curv_mean, curvature_measure &curv_gauss);

	void compute_mean_curvature(PolyMesh &mesh, float scale, curvature_measure &curv);

	float maxdon(const PolyMesh &mesh, PolyMesh::VertexHandle v, float normal_power);

	void compute_maxdon(PolyMesh &mesh, float normal_power, curvature_measure &curv);

	void compute_don(PolyMesh &mesh, OpenMesh::VPropHandleT<float> prop_vert_area, float normal_power, curvature_measure &curv);

	curvature_autocontrast autocontrast(
		const PolyMesh &mesh,
		const curvature_measure &curv,
		float target_entropy_factor,
		OpenMesh::VPropHandleT<float> prop_vert_area,
		bool natural_binning
	);

}

#endif
