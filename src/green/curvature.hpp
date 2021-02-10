
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

#include "model.hpp"
#include "Histogram.h"

#include <random>
#include <algorithm>

namespace green {

	bool computeCurvature(
		PolyMesh & mesh,
		OpenMesh::VPropHandleT<float> &gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty
	);

	bool computeMeanCurvature(
		PolyMesh & mesh,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty,
		lce::Histogram &hCurvature,
		float scale
	);

	bool computeCurvatureSimilarity(
		PolyMesh & mesh,
		OpenMesh::VPropHandleT<float> &gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty,
		OpenMesh::VPropHandleT<float> &similarityCurvatureProperty,
		lce::Histogram &hCurvature
	);

	float maxdon(const PolyMesh& mesh, PolyMesh::VertexHandle v, float normalPower);

	bool computeDoNMaxDiffs(
		PolyMesh &mesh,
		OpenMesh::VPropHandleT<float> &DoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower
	);

	bool computeDoN(
		PolyMesh &mesh,
		OpenMesh::VPropHandleT<float> &DoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower
	);

}

#endif
