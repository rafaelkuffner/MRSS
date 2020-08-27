
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

#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

#include "meshutils.hpp"

namespace green {

	using lce::Histogram;

	float doubleFaceArea(const OpenMesh::Vec3f & e0, const OpenMesh::Vec3f & e2)
	{
		return OpenMesh::cross(e0, -e2).norm();
	}

	bool computeCurvature(TriMesh & mesh,
		OpenMesh::VPropHandleT<float> & gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> & meanCurvatureProperty )
	{
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		mesh.update_face_normals();
		mesh.update_vertex_normals();


		mesh.add_property(gaussianCurvatureProperty);
		mesh.add_property(meanCurvatureProperty);


		TriMesh::VertexIter vIt, vEnd;


		//initialize properties and helpers

		//std::map<TriMesh::VertexHandle, float> areaWeight;
		//std::map<TriMesh::VertexHandle, OpenMesh::Vec3f> accDiffVec;
		std::vector<float> areaWeight(mesh.n_vertices());
		std::vector<OpenMesh::Vec3f> accDiffVec(mesh.n_vertices());

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			mesh.property(gaussianCurvatureProperty, *vIt) = M_PI * 2.0f;
			mesh.property(meanCurvatureProperty, *vIt)     = 0.0f;

			areaWeight[vIt->idx()] = 0.0f;
			accDiffVec[vIt->idx()] = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		}


		//compute vertex weights via face areas

		TriMesh::FaceIter fIt, fEnd;
		TriMesh::FaceVertexIter vfIt;

		OpenMesh::Vec3f p0, p1, p2,
			e0, e1, e2;

		TriMesh::VertexHandle v0, v1, v2;

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

			bool nonmanifold = mesh.status(*fIt).fixed_nonmanifold();
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

				mesh.property(gaussianCurvatureProperty, v0) -= alpha0;
				mesh.property(gaussianCurvatureProperty, v1) -= alpha1;
				mesh.property(gaussianCurvatureProperty, v2) -= alpha2;

				//TODO: add case for adjusting gaussian curvature at borders
			}
		}


		//finally, compute final curvature values for vertices

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			float a = 0;
			if (areaWeight[vIt->idx()] <= std::numeric_limits<float>::epsilon())
			{
				mesh.property(gaussianCurvatureProperty, *vIt) = 0.0f;
				mesh.property(meanCurvatureProperty,     *vIt) = 0.0f;
			}
			else
			{
				mesh.property(gaussianCurvatureProperty, *vIt) /= areaWeight[vIt->idx()];

				OpenMesh::Vec3f a = accDiffVec[vIt->idx()];
				float b = areaWeight[vIt->idx()];
				float k1 = (a / b).norm();
				k = (OpenMesh::dot(accDiffVec[vIt->idx()], mesh.normal(*vIt)) <= 0.0f ? -1.0f : 1.0f) * k1;

				mesh.property(meanCurvatureProperty, *vIt) = k;

			}
		}

		return true;
	}


	bool computeCurvatureSimilarity(TriMesh & mesh,
		OpenMesh::VPropHandleT<float> & gaussianCurvatureProperty,
		OpenMesh::VPropHandleT<float> & meanCurvatureProperty, 
		OpenMesh::VPropHandleT<float> & similarityCurvatureProperty, Histogram &hCurvature)
	{
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		mesh.update_face_normals();
		mesh.update_vertex_normals();


		mesh.add_property(gaussianCurvatureProperty);
		mesh.add_property(meanCurvatureProperty);
		mesh.add_property(similarityCurvatureProperty);

		TriMesh::VertexIter vIt, vEnd;


		//initialize properties and helpers

		//std::map<TriMesh::VertexHandle, float> areaWeight;
		//std::map<TriMesh::VertexHandle, OpenMesh::Vec3f> accDiffVec;
		std::vector<float> areaWeight(mesh.n_vertices());
		std::vector<OpenMesh::Vec3f> accDiffVec(mesh.n_vertices());

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			mesh.property(gaussianCurvatureProperty, *vIt) = M_PI * 2.0f;
			mesh.property(meanCurvatureProperty, *vIt) = 0.0f;

			areaWeight[vIt->idx()] = 0.0f;
			accDiffVec[vIt->idx()] = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		}


		//compute vertex weights via face areas

		TriMesh::FaceIter fIt, fEnd;
		TriMesh::FaceVertexIter vfIt;

		OpenMesh::Vec3f p0, p1, p2,
			e0, e1, e2;

		TriMesh::VertexHandle v0, v1, v2;

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

			bool nonmanifold = mesh.status(*fIt).fixed_nonmanifold();
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


				//accumulate values for gaussian curvature

				mesh.property(gaussianCurvatureProperty, v0) -= alpha0;
				mesh.property(gaussianCurvatureProperty, v1) -= alpha1;
				mesh.property(gaussianCurvatureProperty, v2) -= alpha2;

				//TODO: add case for adjusting gaussian curvature at borders
			}
		}


		//finally, compute final curvature values for vertices

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			float a = 0;
			if (areaWeight[vIt->idx()] <= std::numeric_limits<float>::epsilon())
			{
				mesh.property(gaussianCurvatureProperty, *vIt) = 0.0f;
				mesh.property(meanCurvatureProperty, *vIt) = 0.0f;
			}
			else
			{
				float gaussian = mesh.property(gaussianCurvatureProperty, *vIt) / areaWeight[vIt->idx()];
				mesh.property(gaussianCurvatureProperty, *vIt) = gaussian;

				float val = (accDiffVec[vIt->idx()] / areaWeight[vIt->idx()]).norm();
				k = (OpenMesh::dot(accDiffVec[vIt->idx()], mesh.normal(*vIt)) <= 0.0f ? -1.0f : 1.0f) * val;

				mesh.property(meanCurvatureProperty, *vIt) = k;

				float k1 = k - sqrt(pow(k, 2) - gaussian);
				float k2 = k + sqrt(pow(k, 2) - gaussian);
				float k3 = std::min(abs(k1), abs(k2)) / std::max(abs(k1), abs(k2));
				mesh.property(similarityCurvatureProperty, *vIt) = k3;
				hCurvature.add(mesh.property(similarityCurvatureProperty, *vIt));
			}
		}

		return true;
	}

	bool computeMeanCurvature(TriMesh & mesh,
		OpenMesh::VPropHandleT<float> & meanCurvatureProperty, Histogram &hCurvature,float scale)
	{
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		mesh.update_face_normals();
		mesh.update_vertex_normals();


		mesh.add_property(meanCurvatureProperty);


		TriMesh::VertexIter vIt, vEnd;


		//initialize properties and helpers

		//std::map<TriMesh::VertexHandle, float> areaWeight;
		//std::map<TriMesh::VertexHandle, OpenMesh::Vec3f> accDiffVec;
		std::vector<float> areaWeight(mesh.n_vertices());
		std::vector<OpenMesh::Vec3f> accDiffVec(mesh.n_vertices());

		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			mesh.property(meanCurvatureProperty, *vIt) = 0.0f;

			areaWeight[vIt->idx()] = 0.0f;
			accDiffVec[vIt->idx()] = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		}


		//compute vertex weights via face areas

		TriMesh::FaceIter fIt, fEnd;
		TriMesh::FaceVertexIter vfIt;

		OpenMesh::Vec3f p0, p1, p2,
			e0, e1, e2;

		TriMesh::VertexHandle v0, v1, v2;

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

			bool nonmanifold = mesh.status(*fIt).fixed_nonmanifold();
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
				mesh.property(meanCurvatureProperty, *vIt) = 0.0f;
			}
			else
			{
				OpenMesh::Vec3f a = accDiffVec[vIt->idx()];
				float b = areaWeight[vIt->idx()];
				float k1 = (a / b).norm();
				k = (OpenMesh::dot(accDiffVec[vIt->idx()], mesh.normal(*vIt)) <= 0.0f ? -1.0f : 1.0f) * k1;

				mesh.property(meanCurvatureProperty, *vIt) = k;
				hCurvature.add(k);
			}
		}

		return true;
	}

	float maxdon(const TriMesh& mesh, TriMesh::VertexHandle v, float normalPower) {
		float maxDiff = 1.0f;
		for (TriMesh::ConstVertexOHalfedgeIter vohIt = mesh.cvoh_iter(v); vohIt.is_valid(); ++vohIt) {
			auto nv = mesh.to_vertex_handle(*vohIt);
			auto norm1 = mesh.normal(nv);
			for (TriMesh::ConstVertexOHalfedgeIter vohIt2 = mesh.cvoh_iter(v); vohIt2.is_valid(); ++vohIt2) {
				auto nv2 = mesh.to_vertex_handle(*vohIt2);
				auto norm2 = mesh.normal(nv2);
				maxDiff = std::min(OpenMesh::dot(norm1, norm2), maxDiff);
			}
		}
		maxDiff = pow(((1 - (maxDiff)) / 2), normalPower); //Divided by two to get 0-1 range so our power parameter works :)
		maxDiff = std::clamp(maxDiff, 0.f, 1.f);
		if (!isfinite(maxDiff)) maxDiff = 0;
		return maxDiff;
	}

	bool computeDoNMaxDiffs(
		TriMesh& mesh,
		OpenMesh::VPropHandleT<float>& DoN,
		Histogram& hDoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower
	) {
		assert(mesh.has_vertex_normals());
		mesh.add_property(DoN);
		for (auto vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) {
			float maxDiff = maxdon(mesh, *vIt, normalPower);
			mesh.property(DoN, *vIt) = maxDiff;
			hDoN.add(maxDiff);
		}
		return true;
	}

	bool computeDoN(
		TriMesh & mesh,
		OpenMesh::VPropHandleT<float> & DoN,
		Histogram &hDoN,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower
	) {
		assert(mesh.has_vertex_normals());
		mesh.add_property(DoN);

		TriMesh::VertexIter vIt, vEnd;

		//finally, compute final curvature values for vertices


		for (vIt = mesh.vertices_begin(), vEnd = mesh.vertices_end(); vIt != vEnd; ++vIt)
		{
			float k = 0.0f;
			float w = 0;
			auto norm0 = mesh.normal(*vIt);
			for (TriMesh::ConstVertexOHalfedgeIter vohIt = mesh.cvoh_iter(TriMesh::VertexHandle(vIt)); vohIt.is_valid(); ++vohIt) {

				auto nv = mesh.to_vertex_handle(*vohIt);
				float area = mesh.property(vertexAreasProperty, nv);
				w += area;
				auto norm1 = mesh.normal(nv);
				k += area * std::max(OpenMesh::dot(norm0, norm1), 0.f);

			}

			k = pow((1 - (k / w)),normalPower);
			if (!isfinite(k)) k = 0;
			mesh.property(DoN, *vIt) = k;
			hDoN.add(k);
		}

		return true;
	}

}
