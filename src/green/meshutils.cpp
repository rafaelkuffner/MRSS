
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

#include "meshutils.hpp"

namespace green {

	void hsv2rgb(OpenMesh::Vec3f & rgbColor, const OpenMesh::Vec3f hsvColor)
	{
		unsigned region;
		float remainder;
		float p;
		float q;
		float t;

		float h = hsvColor[0];
		float s = hsvColor[1];
		float v = hsvColor[2];

		if (s == 0)
		{
			rgbColor[0] = v;
			rgbColor[1] = v;
			rgbColor[2] = v;
			return;
		}

		region    = static_cast<unsigned>(floor(h / 0.1666f));
		remainder = (h / 0.1666f) - region;

		p = v * (1.0f - s);
		q = v * (1.0f - (s * remainder));
		t = v * (1.0f - (s * (1.0f - remainder)));

		switch (region)
		{
		case 0:
			rgbColor[0] = v; rgbColor[1] = t; rgbColor[2] = p;
			break;
		case 1:
			rgbColor[0] = q; rgbColor[1] = v; rgbColor[2] = p;
			break;
		case 2:
			rgbColor[0] = p; rgbColor[1] = v; rgbColor[2] = t;
			break;
		case 3:
			rgbColor[0] = p; rgbColor[1] = q; rgbColor[2] = v;
			break;
		case 4:
			rgbColor[0] = t; rgbColor[1] = p; rgbColor[2] = v;
			break;
		default:
			rgbColor[0] = v; rgbColor[1] = p; rgbColor[2] = q;
			break;
		}
	}

	void mapScalarToColor(OpenMesh::Vec3f & rgbColor, float value, TransferFunction transferFunction)
	{
		switch (transferFunction)
		{
		case RAINBOW:
		{
			hsv2rgb(rgbColor, OpenMesh::Vec3f((1.0f - value) * 0.66666f, 1.0f, 1.0f));
			break;
		}
		default:
		case HEAT:
		{
			if (value <= 0.25f)
			{
				value *= 4.0f;
				rgbColor = (1.0f - value) * OpenMesh::Vec3f(0.1f, 0.0f, 0.0f) + //very dark red
					value  * OpenMesh::Vec3f(0.6f, 0.0f, 0.0f);  //a bit dark red
			}
			else if (value <= 0.5f)
			{
				value = 4.0f * (value - 0.25f);
				rgbColor = (1.0f - value) * OpenMesh::Vec3f(0.6f, 0.0f, 0.0f) + //a bit dark red
					value  * OpenMesh::Vec3f(0.9f, 0.2f, 0.0f);  //orange-ish red
			}
			else if (value <= 0.75)
			{
				value = 4.0f * (value - 0.5f);
				rgbColor = (1.0f - value) * OpenMesh::Vec3f(0.9f, 0.2f, 0.0f) + //orange-ish red
					value  * OpenMesh::Vec3f(1.0f, 0.8f, 0.0f);  //orange-ish yellow
			}
			else
			{
				value = 4.0f * (value - 0.75f);
				rgbColor = (1.0f - value) * OpenMesh::Vec3f(1.0f, 0.8f, 0.0f) + //orange-ish yellow
					value  * OpenMesh::Vec3f(1.0f, 1.0f, 0.7f);  //very light yellow
			}
			break;
		}
		case ZBRUSH:
		{
			if (value <= 0.5f)
			{
				value *= 2.0f;
				rgbColor = OpenMesh::Vec3f(value, 1.0f, 1.0f); //cyan to white
			}
			else
			{
				value = 2.0f * (value - 0.5f);
				rgbColor =  OpenMesh::Vec3f(1.0f, (1.0f - value), (1.0f - value));  //white to red
			}

			break;
		}
		}
	}


	OpenMesh::EPropHandleT<float> computeEdgeLengths(TriMesh & mesh)
	{
		OpenMesh::EPropHandleT<float> edgeLengthProperty;
		mesh.add_property(edgeLengthProperty);

		TriMesh::ConstEdgeIter eIt = mesh.edges_begin(), eEnd = mesh.edges_end();
		for (; eIt != eEnd; ++eIt)
		{
			mesh.property(edgeLengthProperty, *eIt) = mesh.calc_edge_length(*eIt);
		}

		return edgeLengthProperty;
	}

	OpenMesh::HPropHandleT<float> computeWedgeVoronoiAreas(TriMesh & mesh)
	{
		OpenMesh::HPropHandleT<float> wedgeVoronoiAreaProperty;

		mesh.add_property(wedgeVoronoiAreaProperty);

		TriMesh::ConstFaceIter fIt = mesh.faces_begin(),
			fEnd = mesh.faces_end();

		TriMesh::HalfedgeHandle halfedge;

		OpenMesh::Vec3f p0, p1, p2,
			e0, e1, e2,
			e0n, e1n,
			n, yT;

		OpenMesh::Vec2f p02D, p12D, p22D,
			e12Dn, e22Dn,
			c0, c1, c2,
			d0, d1, d2,
			centerPoint2D;

		const float SqrLengthEps = 1.0e-8f;

		for (; fIt != fEnd; ++fIt)
		{
			halfedge = *mesh.cfh_iter(*fIt);

			//get points within the triangle
			p0 = mesh.point(mesh.to_vertex_handle(halfedge));
			halfedge = mesh.next_halfedge_handle(halfedge);

			p1 = mesh.point(mesh.to_vertex_handle(halfedge));
			halfedge = mesh.next_halfedge_handle(halfedge);

			p2 = mesh.point(mesh.to_vertex_handle(halfedge));
			halfedge = mesh.next_halfedge_handle(halfedge); //points to initial halfedge again


															//compute mixed voronoi cell centers

															//compute triangle's voronoi cell center, being the midpoint of the enclosing circle
															//this can be done by intersecting any two of the three vectors standing orthogonal
															//on the triangle's edges and originating from their midpoints

			e0 = p1 - p0; //3D edge vector p0 -> p1
			e1 = p2 - p1; //3D edge vector p1 -> p2
			e2 = p0 - p2; //3D edge vector p2 -> p0

			if (e0.sqrnorm() < SqrLengthEps ||
				e1.sqrnorm() < SqrLengthEps ||
				e2.sqrnorm() < SqrLengthEps   )
			{
				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0.0f;
				halfedge = mesh.next_halfedge_handle(halfedge);
				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0.0f;
				halfedge = mesh.next_halfedge_handle(halfedge);
				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0.0f;
				continue;
			}

			e0n = e0.normalized();
			e1n = e1.normalized();

			//project to 2D, using the 3D reference system (e0, yT, n), with origin p0
			//(where the triangle and center both have an always-zero 3rd component)
			n  = (e0n % e1n).normalized(); //3D normal vector
			yT = n   % e0n;                //y axis within the 2D triangle plane

			p02D = OpenMesh::Vec2f(0.0f, 0.0f);
			p12D = OpenMesh::Vec2f(e0.length(), 0.0f);
			p22D = OpenMesh::Vec2f(-(e2 | e0n), -(e2 | yT));

			c0 = 0.5f *         p12D;  //center of edge p02D -> p12D
			c1 = 0.5f * (p12D + p22D); //center of edge p12D -> p22D
			c2 = 0.5f *  p22D;         //center of edge p22D -> p02D

			e12Dn = (p22D - p12D).normalized();
			e22Dn = -c2.normalized(); //minus for consistency of directions

			d0 = OpenMesh::Vec2f( 0.0f,     1.0f    ); //orthogonal vector on e02D
			d1 = OpenMesh::Vec2f(-e12Dn[1], e12Dn[0]); //orthogonal vector on e12D
			d2 = OpenMesh::Vec2f(-e22Dn[1], e22Dn[0]); //orthogonal vector on e22D

			intersectLines2D(c0, d0, c2, d2,
				centerPoint2D);

			if (isinf(centerPoint2D[0]) || isinf(centerPoint2D[1])) {
				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0;
				halfedge = mesh.next_halfedge_handle(halfedge);

				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0;
				halfedge = mesh.next_halfedge_handle(halfedge);

				mesh.property(wedgeVoronoiAreaProperty, halfedge) = 0;

				continue;
			}

			//clamp the cell center against the borders of the triangle

			//we can assume by construction that the point only lies outside w.r.t one particular half-plane
			//if it is outside, clamp against the half-plane's corresponding edge
			if (((centerPoint2D - c0) | d0) < 0.0f)
			{
				centerPoint2D = c0;
			}
			else if (((centerPoint2D - c1) | d1) < 0.0f)
			{
				centerPoint2D = c1;
			}
			else if (((centerPoint2D - c2) | d2) < 0.0f)
			{
				centerPoint2D = c2;
			}


			//compute and store cell areas of each halfedge (corresponding to its vertex' wedge)

			mesh.property(wedgeVoronoiAreaProperty, halfedge) = triangleArea2D(p02D, c0, centerPoint2D) +
				triangleArea2D(p02D, centerPoint2D, c2);

			halfedge = mesh.next_halfedge_handle(halfedge);

			mesh.property(wedgeVoronoiAreaProperty, halfedge) = triangleArea2D(p12D, c1, centerPoint2D) +
				triangleArea2D(p12D, centerPoint2D, c0);

			halfedge = mesh.next_halfedge_handle(halfedge);

			mesh.property(wedgeVoronoiAreaProperty, halfedge) = triangleArea2D(p22D, c2, centerPoint2D) +
				triangleArea2D(p22D, centerPoint2D, c1);

		}


		return wedgeVoronoiAreaProperty;
	}

	OpenMesh::VPropHandleT<float> computeVertexAreas(TriMesh & mesh)
	{
		OpenMesh::VPropHandleT<float> vertexAreaProperty;
		mesh.add_property(vertexAreaProperty);

		//we'll compute vertex areas by summing up the voronoi areas of the wedges of each vertex
		OpenMesh::HPropHandleT<float> wedgeVoronoiAreaProperty = computeWedgeVoronoiAreas(mesh);

		TriMesh::ConstVertexIter vIt = mesh.vertices_begin(),
			vEnd = mesh.vertices_end();

		TriMesh::ConstVertexIHalfedgeIter vihIt;

		float vArea;

		for (; vIt != vEnd; ++vIt)
		{
			//compute and store the area of the vertex by summing up the areas of all neighboring halfedges
			vArea = 0.0f;

			vihIt = mesh.cvih_iter(*vIt);

			while (vihIt.is_valid())
			{
				vArea += mesh.property(wedgeVoronoiAreaProperty, *vihIt);

				++vihIt;
			}

			mesh.property(vertexAreaProperty, *vIt) = vArea;
		}

		return vertexAreaProperty;
	}

	void intersectLines2D(const OpenMesh::Vec2f & p0, const OpenMesh::Vec2f & d0,
		const OpenMesh::Vec2f & p1, const OpenMesh::Vec2f & d1,
		OpenMesh::Vec2f & intersection)
	{
		if (fabs(d0 | d1) == 1.0f)
		{
			intersection = OpenMesh::Vec2f(std::numeric_limits<float>::infinity(),
				std::numeric_limits<float>::infinity());
			return;
		}

		//solve linear system A * x = b by inverting A,
		//using adjugate matrix and inverse determinant

		OpenMesh::Vec2f b = p1 - p0;

		float det = d0[1]*d1[0] - d0[0]*d1[1];

		if (det == 0.0f)
		{
			intersection = OpenMesh::Vec2f(std::numeric_limits<float>::infinity(),
				std::numeric_limits<float>::infinity());
			return;
		}

		//compute intersection point
		//using line parameter for line0 via (line0, line1) = 1/det(A) * adj(A) * b

		det = 1.0f / det;

		intersection = p0 + d0 * (det * (d1[0]*b[1] - d1[1]*b[0]));
	}

	void intersectLines2D(const OpenMesh::Vec2f & p0, const OpenMesh::Vec2f & d0,
		const OpenMesh::Vec2f & p1, const OpenMesh::Vec2f & d1,
		float & lineParam0, float & lineParam1)
	{
		if (fabs(d0 | d1) == 1.0f)
		{
			lineParam0 = lineParam1 = std::numeric_limits<float>::infinity();
			return;
		}

		//solve linear system A * x = b by inverting A,
		//using adjugate matrix and inverse determinant

		OpenMesh::Vec2f b = p1 - p0;

		float det = d0[1]*d1[0] - d0[0]*d1[1];

		if (det == 0.0f)
		{
			lineParam0 = lineParam1 = std::numeric_limits<float>::infinity();
			return;
		}

		//compute x = 1/det(A) * adj(A) * b

		det = 1.0f / det;

		lineParam0 = det * (d1[0]*b[1] - d1[1]*b[0]);
		lineParam1 = det * (d0[0]*b[1] - d0[1]*b[0]);
	}

	float surfaceArea(TriMesh & mesh)
	{
		float area = 0.0f;

		TriMesh::FaceIter fIt, fEnd;

		for (fIt = mesh.faces_begin(), fEnd = mesh.faces_end(); fIt != fEnd; ++fIt)
		{
			//TODO: check why we can't provide *fIt directly as FaceHandle param to faceArea
			TriMesh::FaceHandle face = *fIt;
			area += faceArea(mesh, face);
		}

		return area;
	}

	float faceArea(TriMesh & mesh, TriMesh::FaceHandle & face)
	{
		TriMesh::HalfedgeHandle he0 = mesh.halfedge_handle(face);
		TriMesh::HalfedgeHandle he1 = mesh.next_halfedge_handle(he0);
		TriMesh::HalfedgeHandle he2 = mesh.next_halfedge_handle(he1);

		OpenMesh::Vec3f p0 = mesh.point(mesh.to_vertex_handle(he0));
		OpenMesh::Vec3f p1 = mesh.point(mesh.to_vertex_handle(he1));
		OpenMesh::Vec3f p2 = mesh.point(mesh.to_vertex_handle(he2));

		return triangleArea3D(p0, p1, p2);
	}

}
