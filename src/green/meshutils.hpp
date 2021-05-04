
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

#include <limits>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <glm/glm.hpp>

namespace green {

	inline float powsign(float x, float e) {
		return std::copysign(std::pow(std::abs(x), e), x);
	}

	inline OpenMesh::Vec2f glm2om(const glm::vec2 &v) {
		OpenMesh::Vec2f vv;
		vv[0] = v.x;
		vv[1] = v.y;
		return vv;
	}

	inline glm::vec2 om2glm(const OpenMesh::Vec2f &v) {
		glm::vec2 vv;
		vv.x = v[0];
		vv.y = v[1];
		return vv;
	}

	inline OpenMesh::Vec3f glm2om(const glm::vec3 &v) {
		OpenMesh::Vec3f vv;
		vv[0] = v.x;
		vv[1] = v.y;
		vv[2] = v.z;
		return vv;
	}

	inline glm::vec3 om2glm(const OpenMesh::Vec3f &v) {
		glm::vec3 vv;
		vv.x = v[0];
		vv.y = v[1];
		vv.z = v[2];
		return vv;
	}
	
	struct GreenMeshTraits : OpenMesh::DefaultTraits {
		// we always calc face/vertex normals
		DeclareVertexAttributes(OpenMesh::AttributeBits::Normal);
		DeclareFaceAttributes(OpenMesh::AttributeBits::Normal);
	};

	struct PolyMesh : public OpenMesh::PolyMesh_ArrayKernelT<GreenMeshTraits> {};

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

	float surfaceArea(PolyMesh & mesh);

	float faceArea(PolyMesh & mesh, PolyMesh::FaceHandle & face);

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

	template <typename T>
	inline std::pair<T, T> property_range(const green::PolyMesh &mesh, OpenMesh::VPropHandleT<T> prop) {
		T lo = std::numeric_limits<T>::max();
		T hi = std::numeric_limits<T>::lowest();
		for (auto v : mesh.vertices()) {
			auto x = mesh.property(prop, v);
			lo = std::min(lo, x);
			hi = std::max(hi, x);
		}
		return std::pair(lo, hi);
	}

	OpenMesh::EPropHandleT<float> computeEdgeLengths(PolyMesh & mesh);

	OpenMesh::HPropHandleT<float> computeWedgeVoronoiAreas(PolyMesh & mesh);

	OpenMesh::VPropHandleT<float> computeVertexAreas(PolyMesh & mesh);

}

#endif
