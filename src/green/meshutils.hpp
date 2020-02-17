
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

#ifndef GREEN_MESHUTILS_HPP
#define GREEN_MESHUTILS_HPP

#include "model.hpp"

namespace green {

	enum TransferFunction
	{
		RAINBOW,
		HEAT,
		ZBRUSH
	};

	inline float triangleArea3D(const OpenMesh::Vec3f & p0, const OpenMesh::Vec3f & p1, const OpenMesh::Vec3f & p2)
	{
		return ((p0-p1) % (p0-p2)).norm() * 0.5f;
	}

	inline float triangleArea2D(const OpenMesh::Vec2f & p0, const OpenMesh::Vec2f & p1, const OpenMesh::Vec2f & p2)
	{
		return ((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1])) * 0.5f;
	}

	void mapScalarToColor(OpenMesh::Vec3f & rgbColor, float value, TransferFunction transferFunction = HEAT);

	void intersectLines2D(const OpenMesh::Vec2f & p0, const OpenMesh::Vec2f & d0,
		const OpenMesh::Vec2f & p1, const OpenMesh::Vec2f & d1,
		OpenMesh::Vec2f & intersection);

	float surfaceArea(TriMesh & mesh);

	float faceArea(TriMesh & mesh, TriMesh::FaceHandle & face);

	template<typename T>
	inline T angle(const OpenMesh::VectorT<T, 3> & v, const OpenMesh::VectorT<T, 3> & w)
	{
		float nv = v.norm();
		float nw = w.norm();

		// pow(eps,2) approaches float min
		const float eps = 1e-18f;

		if (nv < eps || nw < eps) return 0.f;

		const float q = std::min(std::max((v|w) / (nv * nw), -1.f), 1.f);
		if (q > 0) {
			// using cross and asin is more accurate at small angles
			// but doesnt work for obtuse angles
			const float r = cross(v, w).norm() / (nv * nw);
			return asin(std::min(std::max(r, 0.f), 1.f));
		} else {
			return acos(q);
		}

	}

	OpenMesh::EPropHandleT<float> computeEdgeLengths(TriMesh & mesh);

	OpenMesh::HPropHandleT<float> computeWedgeVoronoiAreas(TriMesh & mesh);

	OpenMesh::VPropHandleT<float> computeVertexAreas(TriMesh & mesh);

}

#endif
