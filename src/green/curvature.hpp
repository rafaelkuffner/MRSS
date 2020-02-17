
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

namespace green {

	bool computeCurvature(
		TriMesh & mesh,
		OpenMesh::VPropHandleT<float> &gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty
	);

	bool computeMeanCurvature(
		TriMesh & mesh,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty,
		lce::Histogram &hCurvature,
		float scale
	);

	bool computeCurvatureSimilarity(
		TriMesh & mesh,
		OpenMesh::VPropHandleT<float> &gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> &meanCurvatureProperty,
		OpenMesh::VPropHandleT<float> &similarityCurvatureProperty,
		lce::Histogram &hCurvature
	);

	bool computeDoNMaxDiffs(
		TriMesh &mesh,
		OpenMesh::VPropHandleT<float> &DoN,
		lce::Histogram &hDoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPowe
	);

	bool computeDoN(
		TriMesh &mesh,
		OpenMesh::VPropHandleT<float> &DoN,
		lce::Histogram &hDoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower
	);

}

#endif
