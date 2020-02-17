
/*
* Copyright 2020
* Computational Media Innovation Centre
* Victoria University of Wellington
*
*/

#include "neighborhood.hpp"

#include <deque>

namespace green {

	// assumes T is trivially relocatable and nothrow-movable
	template <typename T>
	class simple_queue {
	private:
		struct elem {
			alignas(T) unsigned char data[sizeof(T)];

			T & as_t() {
				return reinterpret_cast<T &>(*this);
			}

			const T & as_t() const {
				return reinterpret_cast<const T &>(*this);
			}
		};

		size_t m_begin = 0;
		size_t m_end = 0;
		size_t m_size = 0;
		size_t m_mask = 0;
		std::vector<elem> m_data;

		// s must be POT
		void resize_impl(size_t s) {
			//if ((s & (s - 1)) != 0) __debugbreak();
			std::vector<elem> data2;
			data2.resize(s);
			if (m_begin < m_end) {
				std::copy(m_data.begin() + m_begin, m_data.begin() + m_end, data2.begin());
			} else {
				std::copy(m_data.begin() + m_begin, m_data.end(), data2.begin());
				std::copy(m_data.begin(), m_data.begin() + m_end, data2.begin() + m_data.size() - m_begin);
			}
			m_data = std::move(data2);
			m_begin = 0;
			m_end = m_size;
			m_mask = s - 1;
		}

	public:
		simple_queue() {}

		simple_queue(simple_queue &&) = default;
		simple_queue & operator=(simple_queue &&) = default;

		simple_queue(const simple_queue &) = delete;
		simple_queue & operator=(const simple_queue &) = delete;

		size_t capacity() const noexcept {
			return m_mask + 1;
		}

		size_t size() const noexcept {
			return m_size;
		}

		bool empty() const noexcept {
			return !m_size;
		}

		void reserve(size_t s) {
			if (s > std::numeric_limits<size_t>::max() / 2) return;
			size_t z = 32;
			while (z < s) z <<= 1;
			resize_impl(z);
		}

		void reserve_additional(size_t s) {
			if (capacity() == 0) m_data.resize(32);
			if (m_size + s > capacity()) resize_impl(capacity() * 2);
		}

		template <typename U>
		void push_back_unchecked(U &&u) {
			::new(static_cast<void *>(&m_data[m_end].as_t())) T(std::forward<U>(u));
			m_end = m_mask & (m_end + 1);
			m_size++;
		}

		template <typename U>
		void push_front_unchecked(U &&u) {
			m_begin = m_mask & (m_begin - 1);
			::new(static_cast<void *>(&m_data[m_begin].as_t())) T(std::forward<U>(u));
			m_size++;
		}

		template <typename U>
		void push_back(U &&u) {
			reserve_additional(1);
			push_back_unchecked(std::forward<U>(u));
		}

		template <typename U>
		void push_front(U &&u) {
			reserve_additional(1);
			push_front_unchecked(std::forward<U>(u));
		}

		T pop_front() {
			assert(m_size);
			T t = std::move(m_data[m_begin].as_t());
			m_data[m_begin].as_t().~T();
			m_begin = m_mask & (m_begin + 1);
			m_size--;
			return t;
		}

		T pop_back() {
			assert(m_size);
			m_end = m_mask & (m_end - 1);
			T t = std::move(m_data[m_end].as_t());
			m_data[m_end].as_t().~T();
			m_size--;
			return t;
		}

		const T & top() const {
			return m_data[m_begin].as_t();
		}

		~simple_queue() {
			while (m_size) {
				pop_front();
			}
		}

	};

	template <typename SimdTraits = simd_traits>
	struct neighborhood_search {

		using mask_t = typename SimdTraits::mask_t;
		using index_t = typename SimdTraits::index_t;
		using dist_t = typename SimdTraits::dist_t;

		static constexpr int NUM_VERTEX_STATES = 2;
		int VERTEX_OPEN = 1;
		int VERTEX_CLOSED = 2;

		struct per_vertex {
			int state = 0;
			int flags = 0;
			dist_t dist{0};
			// TODO maybe split to reduce memory consumption due to alignment
		};

		std::vector<per_vertex> vertexData;
		//ska::flat_hash_map<unsigned, per_vertex> vertexData;

		template <int Phase, typename VisitorF>
		void run(const MeshCache &mesh, CalculationStats &stats, const std::array<unsigned, SimdTraits::simd_size> &rootvdis, float radius, const VisitorF &visitor) {

			using namespace simd;

			constexpr float huge_dist = 9001e19f;
			const index_t invalid_idx{-1};

			//vertexData.clear();
			vertexData.resize(mesh.data.size() >> MeshCache::vdi_shift);
			//vertexData.reserve(1024);

			VERTEX_OPEN += NUM_VERTEX_STATES;
			VERTEX_CLOSED += NUM_VERTEX_STATES;
			// TODO vertex state overflow

			simple_queue<unsigned> openVertices;
			openVertices.reserve(256);

			// initialize distances for root vertices
			for (unsigned i = 0; i < SimdTraits::simd_size; i++) {
				const unsigned vdi = rootvdis[i];
				if (vdi == -1) continue;
				auto &vdata = vertexData[vdi >> MeshCache::vdi_shift];
				if (Phase == 0) vdata.dist = dist_t(huge_dist);
				vdata.dist = simd_insert(vdata.dist, i, 0.f);
			}

			// if scalar, phase 0 only needs to init distances
			if (Phase == 0 && SimdTraits::simd_size == 1) return;

			// add the root vertices to the set of open vertices
			for (unsigned i = 0; i < SimdTraits::simd_size; i++) {
				const unsigned vdi = rootvdis[i];
				if (vdi == -1) continue;
				openVertices.push_back(vdi);
				auto &vdata = vertexData[vdi >> MeshCache::vdi_shift];
				vdata.flags = 0;
				vdata.state = VERTEX_OPEN;
			}

			// while there are open vertices, process them
			int phase_i = 0;
			do {
				if (Phase == 0) {
					phase_i++;
					// phase 0 does a minimal search to try and establish connectivity between root vertices
					// ideally, the root vertices are closely connected (i.e. a surface patch)
					if (phase_i > SimdTraits::simd_size * 3) break;
				}

				const unsigned cvdi = openVertices.pop_front();

				if (openVertices.size()) {
					const unsigned xvdi = openVertices.top();
					_mm_prefetch((const char *)(mesh.data.data() + xvdi), _mm_hint(1));
					_mm_prefetch(64 + (const char *)(mesh.data.data() + xvdi), _mm_hint(1));
					//_mm_prefetch((const char *)(vertexData.data() + (xvdi >> MeshCache::vdi_shift)), 1);
				}

				auto &cData = vertexData[cvdi >> MeshCache::vdi_shift];

				// no longer in queue => closed
				cData.state = VERTEX_CLOSED;

				stats.nh_record_loop(cmplt(cData.dist, dist_t(radius)));

				// add to results
				const mask_t resmask = andnot(cmplt(cData.dist, dist_t(radius)), bits2mask<index_t>(cData.flags));
				if (Phase > 0) if (!none(resmask)) {
					// NOTE: now produces vertex data indices, not ordinary vertex indices
					// mesh.get_vertex(cvdi).vi
					//neighbors.push_back(select(invalid_idx, index_t{int(cvdi)}, resmask));
					if constexpr (Phase == 1) {
						// during search, distance should not be passed to the visitor
						// as the vertex could be revisited later with a smaller distance
						// (visitor is not re-executed in this case)
						visitor(cvdi, resmask);
					} else if constexpr (Phase > 1) {
						// in repeated searches (same or smaller neighbourhood radius)
						// smallest distances are already known
						visitor(cvdi, resmask, cData.dist);
					}
					cData.flags |= mask2bits(resmask);
				}

				// cData ref may now become invalid
				const auto cDist = cData.dist;

				// visit direct neighbors
				const unsigned nedges = mesh.get_vertex(cvdi).nedges;
				openVertices.reserve_additional(nedges);
				for (unsigned i = 0; i < nedges; i++) {
					const unsigned nv = mesh.get_vertex(cvdi).edges[i].ndi;
					const float ecost = mesh.get_vertex(cvdi).edges[i].cost;

					auto &nData = vertexData[nv >> MeshCache::vdi_shift];
					const auto oldState = nData.state;
					const bool unvisited = oldState < VERTEX_OPEN;

					// neighbour dist can end up set arbitrarily large, so we must ensure the thresh is <= radius
					const dist_t costthresh = select(min(nData.dist, dist_t(radius)), dist_t(radius), mask_t(unvisited));
					const dist_t oldcost = select(nData.dist, dist_t(huge_dist), mask_t(unvisited));

					// new best cost to neighbor vertex; repeat searches read cached distances
					const dist_t newCost = Phase < 2 ? min(oldcost, cDist + dist_t(ecost)) : nData.dist;

					// update distances even if we don't visit, so that distances are defined for the neighbours of all
					// the vertices we do visit. this ensures that repeat searches can reliably use the cached neighbour distances.
					// repeat searches don't update distances and will only enqueue vertices once
					if (Phase < 2) nData.dist = newCost;

					const auto m = cmplt(newCost, costthresh);
					if (!none(m)) {
						nData.flags = select(nData.flags, 0, unvisited);
						// only add to queue if not already in it
						//if (oldState != VERTEX_OPEN) {
						//	stats.nh_record_push(m);
						//	openVertices.push_back_unchecked(nv);
						//	nData.state = VERTEX_OPEN;
						//}
						if (oldState < VERTEX_OPEN) {
							stats.nh_record_push(m);
							openVertices.push_back_unchecked(nv);
							nData.state = VERTEX_OPEN;
						} else if (oldState > VERTEX_OPEN) {
							stats.nh_record_push(m);
							openVertices.push_front_unchecked(nv);
							nData.state = VERTEX_OPEN;
						}
					}
				}

			} while (openVertices.size());

		}

	};

	void getGeodesicNeighborhood(
		const MeshCache &mesh,
		CalculationStats &stats,
		const std::array<unsigned, simd_traits::simd_size> &rootvdis,
		float radius,
		std::vector<simd_traits::index_t> & neighbors
	) {
		// so we can use < instead of <=
		// probably doesnt really matter, but should allow more exact testing
		radius = std::nextafter(radius, radius * 2);

		using mask_t = simd_traits::mask_t;
		using index_t = simd_traits::index_t;
		using dist_t = simd_traits::dist_t;

		auto visitor = [&](unsigned vdi, mask_t m) {
			using simd::select;
			neighbors.push_back(select(index_t{-1}, index_t{int(vdi)}, m));
		};

		static thread_local neighborhood_search<> search;
		search.run<0>(mesh, stats, rootvdis, radius, visitor);
		search.run<1>(mesh, stats, rootvdis, radius, visitor);

		// this dumps neighborhoods to file for (manual) verification purposes
		// NOT USABLE WHEN MULTITHREADED
		//std::vector<int> neighbors_scalar;
		//neighbors_scalar.reserve(neighbors.size());
		//for (int i = 0; i < simd_traits::simd_size; i++) {
		//	const unsigned vdi = rootvdis[i];
		//	if (vdi == -1) continue;
		//	neighbors_scalar.clear();
		//	for (auto &x : neighbors) {
		//		using simd::simd_extract;
		//		int y = simd_extract(x, i);
		//		if (y < 0) continue;
		//		// NOTE: converting data index back to vertex index
		//		neighbors_scalar.push_back(mesh.get_vertex(y).vi);
		//	}
		//	std::sort(neighbors_scalar.begin(), neighbors_scalar.end(), [](const auto &a, const auto &b) { return a < b; });
		//	static float file_radius = -9001;
		//	static std::ofstream ofs_neighborhoods;
		//	if (radius != file_radius) {
		//		file_radius = radius;
		//		ofs_neighborhoods.close();
		//		ofs_neighborhoods.open(std::string("./neighborhoods_") + std::to_string(radius) + ".txt", std::ios::trunc);
		//	}
		//	for (auto &a : neighbors_scalar) {
		//		ofs_neighborhoods << a << ' ';
		//	}
		//	ofs_neighborhoods << '\n';
		//}

	}

	template <typename SimdTraits = simd_traits>
	std::array<float, SimdTraits::simd_size> getGeodesicNeighborhoodSaliencyImpl(
		neighborhood_search<SimdTraits> &search,
		const MeshCache &mesh,
		CalculationStats &stats,
		const std::array<unsigned, SimdTraits::simd_size> &rootvdis,
		float radius,
		float curvmin,
		float curvmax,
		bool minimizeSmallChanges
	) {
		// so we can use < instead of <=
		// probably doesnt really matter, but should allow more exact testing
		radius = std::nextafter(radius, radius * 2);

		using mask_t = typename SimdTraits::mask_t;
		using index_t = typename SimdTraits::index_t;
		using dist_t = typename SimdTraits::dist_t;

		const float hist_irange = 255.f / (curvmax - curvmin);

		dist_t hist[256]{};
		dist_t totarea{0};

		dist_t avg_normal_x{};
		dist_t avg_normal_y{};
		dist_t avg_normal_z{};

		auto visitor = [&](unsigned vdi, mask_t m) {
			using simd::select;
#ifdef USE_VERTEX_AREA_WEIGHTING
			const float area = mesh.get_vertex(vdi).props[MeshCache::vertex_prop_area];
#elif
			const float area = 1;
#endif
			totarea = select(totarea, totarea + dist_t{area}, m);
			const float curv = simd::max(simd::min(mesh.get_vertex(vdi).props[MeshCache::vertex_prop_curv], curvmax), curvmin);
			const auto curvbin = intptr_t(hist_irange * (curv - curvmin));
			hist[curvbin] = select(hist[curvbin], hist[curvbin] + dist_t{area}, m);
			const auto norm = mesh.get_vertex(vdi).norm;
			avg_normal_x = select(avg_normal_x, avg_normal_x + norm[0], m);
			avg_normal_y = select(avg_normal_y, avg_normal_y + norm[1], m);
			avg_normal_z = select(avg_normal_z, avg_normal_z + norm[2], m);
		};

		search.run<0>(mesh, stats, rootvdis, radius, visitor);
		search.run<1>(mesh, stats, rootvdis, radius, visitor);

		OpenMesh::Vec3f avg_normal[SimdTraits::simd_size];
		dist_t plane_min{+9001e19}, plane_max{-9001e19};

		if (minimizeSmallChanges) {

			for (int k = 0; k < SimdTraits::simd_size; k++) {
				using simd::simd_extract;
				using simd::simd_insert;
				avg_normal[k][0] = simd_extract(avg_normal_x, k);
				avg_normal[k][1] = simd_extract(avg_normal_y, k);
				avg_normal[k][2] = simd_extract(avg_normal_z, k);
				avg_normal[k].normalize();
				avg_normal_x = simd_insert(avg_normal_x, k, avg_normal[k][0]);
				avg_normal_y = simd_insert(avg_normal_y, k, avg_normal[k][1]);
				avg_normal_z = simd_insert(avg_normal_z, k, avg_normal[k][2]);
			}

			auto visitor2 = [&](unsigned vdi, mask_t m, dist_t r) {
				using simd::select;
				using simd::min;
				using simd::max;
				const auto pos = mesh.get_vertex(vdi).pos;
				dist_t dotprods = avg_normal_x * pos[0] + avg_normal_y * pos[1] + avg_normal_z * pos[2];
				//dist_t dotprods = abs(avg_normal_x * pos[0] + avg_normal_y * pos[1] + avg_normal_z * pos[2] - avg_normal_x * x0 - avg_normal_y * y0 - avg_normal_z * z0);

				plane_min = select(plane_min, min(plane_min, dotprods), m);
				plane_max = select(plane_max, max(plane_max, dotprods), m);
			};

			// TODO don't really have to do a search again
			// could cache the positions from this first search
			search.run<2>(mesh, stats, rootvdis, radius, visitor2);

		}

		auto plane_range = plane_max - plane_min;



		const float nilog2 = -1.f / log(2.f);
		const auto itotarea = dist_t{1.f} / totarea;
		std::array<float, SimdTraits::simd_size> score{};

		for (auto a : hist)
		{
			// if a == 0 this bin of the histogram is empty
			const auto pp = a * itotarea;
			// limit of p*log(p) as p tends to 0 is 0
			// so we can safely discard anything very small
			for (int i = 0; i < SimdTraits::simd_size; i++) {
				using simd::simd_extract;
				const auto p = simd_extract(pp, i);
				if (!(p > 1e-6f)) continue;
				score[i] += p * log(p) * nilog2;
			}
		}

		if (minimizeSmallChanges) {
			for (int i = 0; i < SimdTraits::simd_size; i++) {
				using simd::simd_extract;
				float nmapWeight = std::min(pow(simd_extract(plane_range, i) / radius + 0.65f, 10.f), 1.f);
				score[i] *= nmapWeight;
			}
		}

		return score;
	}

	std::array<float, simd_traits::simd_size> getGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		const std::array<unsigned, simd_traits::simd_size> &rootvdis,
		float radius,
		float curvmin,
		float curvmax,
		bool minimizeSmallChanges
	) {
		static thread_local neighborhood_search<> search;
		return getGeodesicNeighborhoodSaliencyImpl(search, mesh, stats, rootvdis, radius, curvmin, curvmax, minimizeSmallChanges);
	}

	void subsampleGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		unsigned rootvdi,
		float radius,
		float curvmin,
		float curvmax,
		bool minimizeSmallChanges,
		float distribution_radius,
		const std::function<void(unsigned vdi, float r, float s)> &distribute_saliency
	) {
		static thread_local neighborhood_search<simd1_traits> search;
		const float s = getGeodesicNeighborhoodSaliencyImpl<simd1_traits>(search, mesh, stats, {rootvdi}, radius, curvmin, curvmax, minimizeSmallChanges)[0];
		using mask_t = simd1_traits::mask_t;
		using dist_t = simd1_traits::dist_t;
		if (distribution_radius <= radius) {
			// if distribution will visit a sub-neighbourhood, can use optimized search
			auto visitor = [&](unsigned vdi, mask_t m, dist_t r) {
				distribute_saliency(vdi, r, s);
			};
			search.run<2>(mesh, stats, {rootvdi}, distribution_radius, visitor);
		} else {
			// otherwise, have to use ordinary search
			// (shouldn't happen normally, but can when testing extreme subsampling)
			static thread_local std::vector<int> nvdis;
			nvdis.clear();
			auto visitor = [&](unsigned vdi, mask_t m) {
				nvdis.push_back(vdi);
			};
			search.run<1>(mesh, stats, {rootvdi}, distribution_radius, visitor);
			for (auto vdi : nvdis) {
				auto &vdata = search.vertexData[vdi >> MeshCache::vdi_shift];
				distribute_saliency(vdi, vdata.dist, s);
			}
		}
	}


	MeshCache::MeshCache(
		const TriMesh & mesh,
		OpenMesh::EPropHandleT<float> edgeLengthProperty,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		OpenMesh::VPropHandleT<float> curvatureMeasure
	) {

		std::cout << "Preparing mesh data for neighborhood search" << std::endl;

		vi2di.resize(mesh.n_vertices());

		std::vector<unsigned char> visited(mesh.n_vertices(), false);
		std::deque<int> open_vertices_outer;
		std::deque<int> open_vertices_inner;

		unsigned prev_di = -1;
		int vertices_added = 0;
		int inner_loop_count = 0;
		size_t total_edges = 0;

		auto alloc_data_index = [&]() {
			const unsigned di = data.size();
			vertices_added++;
			if ((prev_di >> vdi_shift) == (di >> vdi_shift)) std::abort();
			prev_di = di;
			vdis.push_back(di);
			return di;
		};

		auto dovertex = [&](int vi0) {
			open_vertices_outer.clear();
			open_vertices_outer.push_back(vi0);
			for (int i = 0; i < 16 && open_vertices_outer.size(); i++) {
				const int ovi = open_vertices_outer.front();
				open_vertices_outer.pop_front();
				if (visited[ovi]) continue;
				inner_loop_count++;
				open_vertices_inner.clear();
				open_vertices_inner.push_back(ovi);
				do {
					const int vi = open_vertices_inner.front();
					open_vertices_inner.pop_front();
					if (visited[vi]) continue;
					visited[vi] = true;
					const unsigned di = alloc_data_index();
					vi2di[vi] = di;
					static_assert(sizeof(vertex) == 16 + 12 + 12, "");
					static_assert(sizeof(edge) == 8, "");
					for (int k = 0; k < sizeof(vertex) / sizeof(unsigned); k++) {
						data.push_back(0);
					}
					reinterpret_cast<vertex &>(data[di]).vi = vi;
					reinterpret_cast<vertex &>(data[di]).props[vertex_prop_area] = mesh.property(vertexAreasProperty, TriMesh::VertexHandle(vi));
					reinterpret_cast<vertex &>(data[di]).props[vertex_prop_curv] = mesh.property(curvatureMeasure, TriMesh::VertexHandle(vi));
					reinterpret_cast<vertex &>(data[di]).pos = mesh.point(TriMesh::VertexHandle(vi));
					reinterpret_cast<vertex &>(data[di]).norm = mesh.normal(TriMesh::VertexHandle(vi));
					for (TriMesh::ConstVertexOHalfedgeIter vohIt = mesh.cvoh_iter(TriMesh::VertexHandle(vi)); vohIt.is_valid(); ++vohIt) {
						const int nvi = mesh.to_vertex_handle(*vohIt).idx();
						if (!visited[nvi]) open_vertices_inner.push_back(nvi);
						data.push_back(0);
						data.push_back(0);
						// can't resolve neighbour data index on first pass, so use a second pass
						// careful, as push_back will invalidate vertex refs
						reinterpret_cast<vertex &>(data[di]).nedges++;
						total_edges++;
					}
					// pad data to min 6 edges (so entire vertex is at least 64 bytes; ref vdi_shift)
					// 6 edges seems about average for triangulated meshes
					for (int k = reinterpret_cast<vertex &>(data[di]).nedges; k < 6; k++) {
						data.push_back(0);
						data.push_back(0);
					}
					// loop until we have a simd-aligned set of near neighbors
				} while (open_vertices_inner.size() && (vertices_added % simd_traits::simd_size != 0));
				// carry over inner queue
				open_vertices_outer.insert(open_vertices_outer.begin(), open_vertices_inner.begin(), open_vertices_inner.end());
			}
		};

		for (int vi = 0; vi < mesh.n_vertices(); vi++) {
			if (visited[vi]) continue;
			dovertex(vi);
		}

		// resolve vertex neighbor indices
		for (int vi = 0; vi < mesh.n_vertices(); vi++) {
			const unsigned di = vi2di[vi];
			auto &v = reinterpret_cast<vertex &>(data[di]);
			int i = 0;
			for (TriMesh::ConstVertexOHalfedgeIter vohIt = mesh.cvoh_iter(TriMesh::VertexHandle(vi)); vohIt.is_valid(); ++vohIt, ++i) {
				const int nV = mesh.to_vertex_handle(*vohIt).idx();
				auto &e = v.edges[i];
				e.cost = mesh.property(edgeLengthProperty, mesh.edge_handle(*vohIt));
				e.ndi = vi2di[nV];
			}
		}

		std::cout << "Edges per vertex: " << (double(total_edges) / mesh.n_vertices()) << std::endl;
		std::cout << "Inner layout grouping: " << (float(mesh.n_vertices()) / inner_loop_count)
			<< " / " << simd_traits::simd_size << std::endl;

		const size_t bytes_per_thread = (data.size() >> vdi_shift) * sizeof(neighborhood_search<>::per_vertex);
		std::cout << "Neighborhood search memory per thread: " << (bytes_per_thread / (1024 * 1024)) << "MiB" << std::endl;


	}

	void CalculationStats::dump_stats(int threadcount) const {

		using namespace std;

		cout << "Neighborhood loop count: " << nh_loop_count
			<< ", loop occupancy: " << (float(nh_loop_occupancy) / nh_loop_count)
			<< ", push occupancy: " << (float(nh_push_occupancy) / nh_push_count)
			<< endl;

		cout << "Neighborhood: " << (nh_duration / (threadcount * 1.0s)) << "s (total " << (nh_duration / 1.0s) << "s)" << endl;
		//cout << "Saliency    : " << (sal_duration / (threadcount * 1.0s)) << "s (total " << (sal_duration / 1.0s) << "s)" << endl;

	}

}
