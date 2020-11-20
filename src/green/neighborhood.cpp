
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
		void push_back_unchecked(U &&u) noexcept {
			::new(static_cast<void *>(&m_data[m_end].as_t())) T(std::forward<U>(u));
			m_end = m_mask & (m_end + 1);
			m_size++;
		}

		template <typename U>
		void push_front_unchecked(U &&u) noexcept {
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

		const T & front() const noexcept {
			return m_data[m_begin].as_t();
		}

		~simple_queue() {
			while (m_size) {
				pop_front();
			}
		}

	};

	struct simple_neighborhood_queue {
		simple_queue<unsigned> queue;

		simple_neighborhood_queue() {
			queue.reserve(256);
		}

		size_t size() const noexcept {
			return queue.size();
		}

		unsigned pop() noexcept {
			return queue.pop_front();
		}

		unsigned top() const noexcept {
			return queue.front();
		}

		void reserve_neighbors(size_t s) {
			queue.reserve_additional(s);
		}

		bool push(unsigned x) {
			queue.push_back(x);
			return true;
		}

		bool push_neighbor(unsigned x, int statediff) {
			// only add to queue if not already in it
			if (statediff < 0) {
				queue.push_back_unchecked(x);
				return true;
			} else if (statediff > 0) {
				queue.push_front_unchecked(x);
				return true;
			}
			return false;
		}
	};

	template <typename SimdTraits = simd_traits>
	struct neighborhood_search {

		using mask_t = typename SimdTraits::mask_t;
		using index_t = typename SimdTraits::index_t;
		using dist_t = typename SimdTraits::dist_t;
		using neighborhood_queue_t = typename SimdTraits::neighborhood_queue_t;

		static constexpr auto null_visitor = [](unsigned vdi, mask_t m) {};

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

		template <int Phase, typename VisitorF = decltype(null_visitor)>
		void run(const MeshCache &mesh, CalculationStats &stats, const std::array<unsigned, SimdTraits::simd_size> &rootvdis, float radius, const VisitorF &visitor) {

			using namespace simd;

			constexpr float huge_dist = 9001e19f;
			const index_t invalid_idx{-1};

			//vertexData.clear();
			vertexData.resize(mesh.vdi_uid(mesh.data.size()) + 1);
			//vertexData.reserve(1024);

			VERTEX_OPEN += NUM_VERTEX_STATES;
			VERTEX_CLOSED += NUM_VERTEX_STATES;
			// TODO vertex state overflow

			neighborhood_queue_t open_vertices;

			// initialize distances for root vertices
			for (unsigned i = 0; i < SimdTraits::simd_size; i++) {
				const unsigned vdi = rootvdis[i];
				if (vdi == -1) continue;
				auto &vdata = vertexData[mesh.vdi_uid(vdi)];
				if (Phase == 0) vdata.dist = dist_t(huge_dist);
				vdata.dist = simd_insert(vdata.dist, i, 0.f);
			}

			// if scalar, phase 0 only needs to init distances
			if (Phase == 0 && SimdTraits::simd_size == 1) return;

			// add the root vertices to the set of open vertices
			for (unsigned i = 0; i < SimdTraits::simd_size; i++) {
				const unsigned vdi = rootvdis[i];
				if (vdi == -1) continue;
				open_vertices.push(vdi);
				auto &vdata = vertexData[mesh.vdi_uid(vdi)];
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

				const unsigned cvdi = open_vertices.pop();

				if (open_vertices.size()) {
					const unsigned xvdi = open_vertices.top();
					_mm_prefetch((const char *)(mesh.data.data() + xvdi), _mm_hint(1));
					_mm_prefetch(64 + (const char *)(mesh.data.data() + xvdi), _mm_hint(1));
					//_mm_prefetch((const char *)(vertexData.data() + (xvdi >> MeshCache::vdi_shift)), 1);
				}

				auto &cData = vertexData[mesh.vdi_uid(cvdi)];
				auto &vert = mesh.get_vertex(cvdi);

				// no longer in queue => closed
				cData.state = VERTEX_CLOSED;

				stats.nh_record_loop(cmplt(cData.dist, dist_t(radius)));

				// add to results
				const mask_t resmask = andnot(cmplt(cData.dist, dist_t(radius)), bits2mask<index_t>(cData.flags));
				if (Phase > 0) if (!none(resmask)) {
					// NOTE: now produces vertex data indices, not ordinary vertex indices
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
				const unsigned nedges = vert.nedges;
				open_vertices.reserve_neighbors(nedges);
				for (unsigned i = 0; i < nedges; i++) {
					const auto &edge = vert.edges[i];
					const float ecost = mesh.decode_cost(edge.costx);

					auto &nData = vertexData[mesh.vdi_uid(edge.ndi)];
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
						if (open_vertices.push_neighbor(edge.ndi, oldState - VERTEX_OPEN)) {
							stats.nh_record_push(m);
							nData.state = VERTEX_OPEN;
						}
					}
				}

			} while (open_vertices.size());

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
		float noise_height
	) {
		// so we can use < instead of <=
		// probably doesnt really matter, but should allow more exact testing
		radius = std::nextafter(radius, radius * 2);

		using mask_t = typename SimdTraits::mask_t;
		using index_t = typename SimdTraits::index_t;
		using dist_t = typename SimdTraits::dist_t;

		dist_t hist[256]{};
		dist_t totarea{0};

		dist_t avg_normal_x{};
		dist_t avg_normal_y{};
		dist_t avg_normal_z{};

		auto visitor = [&](unsigned vdi, mask_t m) {
			using simd::select;
#ifdef USE_VERTEX_AREA_WEIGHTING
			const float area = mesh.decode_area(mesh.get_vertex(vdi).areax);
#elif
			const float area = 1;
#endif
			totarea = select(totarea, totarea + dist_t{area}, m);
			const intptr_t curvbin = mesh.get_vertex(vdi).curvbin;
			hist[curvbin] = select(hist[curvbin], hist[curvbin] + dist_t{area}, m);
			const auto norm = mesh.get_vertex_aux(vdi).norm;
			avg_normal_x = select(avg_normal_x, avg_normal_x + norm[0], m);
			avg_normal_y = select(avg_normal_y, avg_normal_y + norm[1], m);
			avg_normal_z = select(avg_normal_z, avg_normal_z + norm[2], m);
		};

		search.run<0>(mesh, stats, rootvdis, radius, visitor);
		search.run<1>(mesh, stats, rootvdis, radius, visitor);

		OpenMesh::Vec3f avg_normal[SimdTraits::simd_size];
		dist_t plane_min{+9001e19}, plane_max{-9001e19};

		if (noise_height > 0) {

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
				const auto pos = mesh.get_vertex_aux(vdi).pos;
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

		if (noise_height > 0) {
			for (int i = 0; i < SimdTraits::simd_size; i++) {
				using simd::simd_extract;
                //float nmapWeight = std::min(pow(simd_extract(plane_range, i) / radius + 0.65f, 10.f), 1.f);
				// FIXME smoothstep ?
				float nmapWeight = std::min(simd_extract(plane_range, i) / noise_height, 1.f);
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
		float noise_height
	) {
		static thread_local neighborhood_search<> search;
		return getGeodesicNeighborhoodSaliencyImpl(search, mesh, stats, rootvdis, radius, noise_height);
	}

	void subsampleGeodesicNeighborhoodSaliency(
		const MeshCache &mesh,
		CalculationStats &stats,
		unsigned rootvdi,
		float radius,
		float noise_height,
		float distribution_radius,
		const std::function<void(unsigned vdi, float r, float s)> &distribute_saliency
	) {
		static thread_local neighborhood_search<simd1_traits> search;
		const float s = getGeodesicNeighborhoodSaliencyImpl<simd1_traits>(search, mesh, stats, {rootvdi}, radius, noise_height)[0];
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
				auto &vdata = search.vertexData[mesh.vdi_uid(vdi)];
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
		vdis.reserve(mesh.n_vertices());

		// in data units (unsigned)
		auto vertex_size = [&](unsigned nedges) {
			// pad data so entire vertex is at least the specified size
			// 6 edges seems about average for triangulated meshes
			// note: rounding up to alignment
			const unsigned zvert = (sizeof(vertex) + sizeof(edge) * nedges + sizeof(unsigned) - 1) / sizeof(unsigned);
			return std::max(zvert, vert_min_size);
		};

		// work out precisely how big the data vectors actually need to be
		size_t datasize = 0;
		for (auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
			unsigned nedges = 0;
			for (auto vohIt = mesh.cvoh_iter(vit); vohIt.is_valid(); ++vohIt) {
				nedges++;
				// get max cost for scaling
				const float cost = mesh.property(edgeLengthProperty, mesh.edge_handle(*vohIt));
				cost_scale = std::max(cost_scale, cost);
			}
			datasize += vertex_size(nedges);
			// get max vertex area for scaling
			const float area = mesh.property(vertexAreasProperty, TriMesh::VertexHandle(vit));
			area_scale = std::max(area_scale, area);
		}
		
		// alloc data vectors
		data.resize(datasize);
		dataaux.resize((vdi_uid(datasize) + 1) * (sizeof(vertex_aux) / sizeof(unsigned)));

		std::vector<unsigned char> visited(mesh.n_vertices(), false);
		std::deque<int> open_vertices_outer;
		std::deque<int> open_vertices_inner;

		unsigned prev_di = -1;
		int vertices_added = 0;
		int inner_loop_count = 0;
		size_t total_edges = 0;

		auto alloc_data_index = [&]() {
			unsigned di = 0;
			if (prev_di != -1) di = prev_di + vertex_size(get_vertex(prev_di).nedges);
			if (vdi_uid(prev_di) == vdi_uid(di)) std::abort();
			prev_di = di;
			vdis.push_back(di);
			vertices_added++;
			return di;
		};

		auto dovertex = [&](int vi0) {
			open_vertices_outer.clear();
			open_vertices_outer.push_back(vi0);
			for (int i = 0; i < 128 && open_vertices_outer.size(); ) {
				const int ovi = open_vertices_outer.front();
				open_vertices_outer.pop_front();
				if (visited[ovi]) continue;
				// only count unvisited to group size
				i++;
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
					// note: push_back will invalidate vertex refs (although we're not doing that anymore)
					{
						vertex &vdata = get_vertex(di);
						vertex_aux &vadata = get_vertex_aux(di);
						const float area = std::max(mesh.property(vertexAreasProperty, TriMesh::VertexHandle(vi)), 0.f);
						// encode area
						vdata.areax = encode_area(area);
						const float curv = std::clamp(mesh.property(curvatureMeasure, TriMesh::VertexHandle(vi)), 0.f, 1.f);
						// bin curvature based on 0-1 range
						vdata.curvbin = unsigned(255.f * curv);
						vadata.vi = vi;
						vadata.pos = mesh.point(TriMesh::VertexHandle(vi));
						vadata.norm = mesh.normal(TriMesh::VertexHandle(vi));
					}
					unsigned nedges = 0;
					for (auto vohIt = mesh.cvoh_iter(TriMesh::VertexHandle(vi)); vohIt.is_valid(); ++vohIt) {
						const int nvi = mesh.to_vertex_handle(*vohIt).idx();
						if (!visited[nvi]) open_vertices_inner.push_back(nvi);
						// can't resolve neighbour data index on first pass, so use a second pass
						nedges++;
						total_edges++;
					}
					// FIXME what do about this?
					// possibly just drop the excess edges and emit a warning
					if (nedges > 255) std::abort();
					get_vertex(di).nedges = nedges;
					// loop until we have a simd-aligned set of near neighbors
				} while (open_vertices_inner.size() && (vertices_added % simd_traits::simd_size != 0));
				// carry over inner queue
				open_vertices_outer.insert(open_vertices_outer.begin(), open_vertices_inner.begin(), open_vertices_inner.end());
			}
		};

		for (auto v : mesh.vertices()) {
			const auto vi = v.idx();
			if (visited[vi]) continue;
			dovertex(vi);
		}

		// resolve vertex neighbor indices
		for (int vi = 0; vi < mesh.n_vertices(); vi++) {
			const unsigned vdi = vi2di[vi];
			auto &v = get_vertex(vdi);
			int i = 0;
			for (auto vohIt = mesh.cvoh_iter(TriMesh::VertexHandle(vi)); vohIt.is_valid(); ++vohIt, ++i) {
				const int nV = mesh.to_vertex_handle(*vohIt).idx();
				edge &e = v.edges[i];
				e.ndi = vi2di[nV];
				e.costx = encode_cost(mesh.property(edgeLengthProperty, mesh.edge_handle(*vohIt)));
			}
		}

		std::cout << "Edges per vertex: " << (double(total_edges) / mesh.n_vertices()) << std::endl;
		std::cout << "Inner layout grouping: " << (float(mesh.n_vertices()) / inner_loop_count)
			<< " / " << simd_traits::simd_size << std::endl;

		const size_t bytes_mesh_data_ideal = mesh.n_vertices() * sizeof(vertex) + mesh.n_edges() * 2 * sizeof(edge);
		const size_t bytes_mesh_data_real = data.size() * sizeof(unsigned);
		const size_t bytes_mesh_dataaux_real = dataaux.size() * sizeof(unsigned);
		std::cout << "Mesh data size: " << (bytes_mesh_data_real / (1024 * 1024)) << "MiB" << std::endl;
		std::cout << "Mesh data padding: " << ((bytes_mesh_data_real - bytes_mesh_data_ideal) / (1024 * 1024)) << "MiB" << std::endl;
		std::cout << "Mesh aux data size: " << (bytes_mesh_dataaux_real / (1024 * 1024)) << "MiB" << std::endl;

		const size_t bytes_per_thread_full = (vdi_uid(data.size()) + 1) * sizeof(neighborhood_search<>::per_vertex);
		const size_t bytes_per_thread_subsampled = (vdi_uid(data.size()) + 1) * sizeof(neighborhood_search<simd1_traits>::per_vertex);
		std::cout << "Neighborhood search memory per thread (full): " << (bytes_per_thread_full / (1024 * 1024)) << "MiB" << std::endl;
		std::cout << "Neighborhood search memory per thread (subsampled): " << (bytes_per_thread_subsampled / (1024 * 1024)) << "MiB" << std::endl;

	}

	void MeshCache::dump_to_file() const {
		std::cout << "Dumping meshcache" << std::endl;
		// TODO custom filename
		std::ofstream out("./meshcache.txt");
		for (int i = 0; i < vi2di.size(); i++) {
			int vdi = vi2di[i];
			out << i << ":";
			auto &v = get_vertex(vdi);
			auto &va = get_vertex_aux(vdi);
			out << " a=" << decode_area(v.areax);
			out << " c=" << v.curvbin;
			out << " p=(" << va.pos[0] << "," << va.pos[1] << "," << va.pos[2] << ")";
			out << " n=(" << va.norm[0] << "," << va.norm[1] << "," << va.norm[2] << ")";
			out << " [";
			for (int j = 0; j < v.nedges; j++) {
				int ndi = v.edges[j].ndi;
				int ni = get_vertex_aux(ndi).vi;
				out << ni << ":";
				out << " l=" << decode_cost(v.edges[i].costx);
				out << "; ";
			}
			out << "]\n";
		}
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
