//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Copyright (C) 2017 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include "mesh_processing.h"
#include <set>
#include <cmath>
namespace mesh_processing {

	using surface_mesh::Point;
	using surface_mesh::Scalar;
	using surface_mesh::Color;
	using std::min;
	using std::max;
	using std::cout;
	using std::endl;

	MeshProcessing::MeshProcessing(const string& filename) {
		load_mesh(filename);
	}

	MeshProcessing::~MeshProcessing() {
		// TODO
	}

	void MeshProcessing::remesh(const REMESHING_TYPE &remeshing_type,
		const int &num_iterations) {
		calc_weights();
		calc_mean_curvature();
		calc_uniform_mean_curvature();
		calc_gauss_curvature();
		calc_target_length(remeshing_type);

		// main remeshing loop
		for (int i = 0; i < num_iterations; ++i)
		{
			split_long_edges();
			collapse_short_edges();
			equalize_valences();
			tangential_relaxation();
		}
	}

	void MeshProcessing::calc_target_length(const REMESHING_TYPE &remeshing_type) {
		Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
		Mesh::Vertex_around_vertex_circulator  vv_c, vv_end;
		Scalar                   length;
		Scalar                   mean_length;
		Scalar                   H;
		Scalar                   K;

		Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:meancurvature", 0);
		Mesh::Vertex_property<Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
		Mesh::Vertex_property<Scalar> target_new_length = mesh_.vertex_property<Scalar>("v:newlength", 0);

		const float TARGET_LENGTH = 2.0;

		if (remeshing_type == AVERAGE)
		{
			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
				target_length[*v_it] = TARGET_LENGTH;

		}
		else if (remeshing_type == CURV)
		{
			// ------------- IMPLEMENT HERE ---------
			// Get the maximal curvature at each vertex (use the precomputed mean and gaussian curvature)
			// Calculate the desired edge length as the TARGET_LENGTH divided by the maximal curvature at each vertex, and assign it to the property target_length
			// Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
			// Rescale the property target_new_length such that it's mean equals the user specified TARGET_LENGTH
			// ------------- IMPLEMENT HERE ---------

			Mesh::Vertex_property<Scalar> curvature1 = mesh_.vertex_property<Scalar>("v:curvature", 0);
			Mesh::Vertex_property<Scalar> gauss_curvature1 = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0);

			calc_mean_curvature();
			calc_gauss_curvature();

			// calculate desired length
			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
				length = 1.0;
				if (!mesh_.is_boundary(*v_it)) {
					length = 1 / (curvature1[*v_it] + sqrt((curvature1[*v_it] * curvature1[*v_it]) - gauss_curvature1[*v_it]));
					//std::cout << "length: " << length << std::endl;
				}
				target_length[*v_it] = length;
			}

			Scalar laplace;
			int numVertices;
			// smooth desired length
			for (int i = 0; i < 5; i++) {
				for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
					// Initialize variables
					vv_c = mesh_.vertices(*v_it);
					if (!vv_c) {
						continue;
					}
					vv_end = vv_c;
					const Scalar& refPoint = target_length[*v_it];
					laplace = 0.0f;
					numVertices = 0;

					// Iterate over adjacent vertices    
					do {
						++numVertices;
						const Scalar& vi = target_length[*vv_c];
						laplace += vi - refPoint;
					} while (++vv_c != vv_end);

					// Average and normalize the Laplacian
					//std::cout << "laplace: " << laplace << std::endl;
					laplace /= numVertices;
					//std::cout << "laplace 2 : " << laplace << std::endl;
					target_length[*v_it] = target_length[*v_it] + laplace;
				}

			}

			// rescale desired length
			Scalar totalLength = 0;

			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
				totalLength += target_length[*v_it];
			}

			//std::cout << "total length: " << totalLength << std::endl;

			Scalar avrgLength = totalLength / mesh_.n_vertices();

			std::cout << "avrg length: " << avrgLength << std::endl;
			std::cout << "tagret length: " << TARGET_LENGTH << std::endl;

			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
				target_length[*v_it] = target_length[*v_it] * (TARGET_LENGTH / avrgLength);
			}
		}
	}

	void MeshProcessing::split_long_edges() {
		

		/* ===============
			TODO: update the edge set while iterating ok ? If no store all modifs in a vector
			and update everything at once at the end of the current iteration
		   =============== */

		// Maximum number of iterations
		const unsigned int MAX_IT = 100;

		// Get properties
		Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
		
		unsigned int i = 0;
		bool finished = false;
		for (; !finished && i < MAX_IT; ++i) {

			// Finished is true by default
			finished = true;

			// Iterate over edges
			Mesh::Edge_iterator e_end{ mesh_.edges_end() };
			for (Mesh::Edge_iterator e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {

				Mesh::Edge cEdge = *e_it;

				// Get both endpoints, edge length, and target length
				Mesh::Vertex v0 = mesh_.vertex(cEdge, 0);
				Mesh::Vertex v1 = mesh_.vertex(cEdge, 1);
				Scalar targetLength = (target_length[v0] + target_length[v1]) / 2;
				Scalar edgeLength = mesh_.edge_length(cEdge);

				// Check splitting condition
				if (edgeLength > ((4.f / 3.f) * targetLength)) {

					// Finished is set to false when an edge is splitted
					finished = false;

					// Create edge midpoint
					Point p = (mesh_.position(v0) + mesh_.position(v1)) / 2;
					Mesh::Vertex v = mesh_.add_vertex(p);

					// Interpolate normal and target length
					normals[v] = (normals[v0] + normals[v1]) / 2;
					target_length[v] = (target_length[v0] + target_length[v1]) / 2;

					// Split the edge and update the mesh datastructure
					mesh_.split(cEdge, v);
				}
			}
		}
	}

	void MeshProcessing::collapse_short_edges() {

		/* ===============
			TODO: update the edge set while iterating ok ? If no store all modifs in a vector
			and update everything at once at the end of the current iteration
			TODO: in each way does mesh_.collapse() work ?
		   =============== */

		   // Maximum number of iterations
		const unsigned int MAX_IT = 100;

		// Get properties
		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

		unsigned int i = 0;
		bool finished = false;
		for (; !finished && i < MAX_IT; ++i) {

			// Finished is true by default
			finished = true;

			// Iterate over edges
			Mesh::Edge_iterator e_end{ mesh_.edges_end() };
			for (Mesh::Edge_iterator e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {

				if (!mesh_.is_deleted(*e_it)) { // The edge might already have been deleted

					// ------------- IMPLEMENT HERE ---------
					// INSERT CODE:
					// Compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
					// If the edge is shorter than 4/5 of the desired length
					//		Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse.
					//		Check if halfedges collapsible
					//		Select the halfedge to be collapsed if at least one halfedge can be collapsed
					//		Collapse the halfedge
					// Leave the loop running until no collapse has been done (use the finished variable)
					// ------------- IMPLEMENT HERE ---------	

					Mesh::Edge cEdge = *e_it;

					// Get both endpoints, edge length, and target length
					Mesh::Vertex v0 = mesh_.vertex(cEdge, 0);
					Mesh::Vertex v1 = mesh_.vertex(cEdge, 1);
					Scalar targetLength = (target_length[v0] + target_length[v1]) / 2;
					Scalar edgeLength = mesh_.edge_length(cEdge);

					// Check collapsing condition
					if (edgeLength < ((4.f / 5.f) * targetLength)) {

						// Check which halfedge is collapsable
						Mesh::Halfedge h0 = mesh_.halfedge(cEdge, 0);
						Mesh::Halfedge h1 = mesh_.halfedge(cEdge, 1);
						bool collapse0 = mesh_.is_collapse_ok(h0);
						bool collapse1 = mesh_.is_collapse_ok(h1);

						// Finished is set to false when at least one of the halfedge is collapsable
						finished = !(collapse0 || collapse1);

						// Collapse the correct halfedge
						if (collapse0 && collapse1) {
							if (mesh_.valence(v0) > mesh_.valence(v1))	mesh_.collapse(h0);
							else										mesh_.collapse(h1);
						}
						else {
							if (collapse0)								mesh_.collapse(h0);
							else if (collapse1)							mesh_.collapse(h1);
						}
					}
				}
			}
		}

		mesh_.garbage_collection();

		if (i == MAX_IT && !finished) std::cerr << "collapse break\n";
	}

	void MeshProcessing::equalize_valences()
	{
		Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
		Mesh::Vertex   v0, v1, v2, v3;
		Mesh::Halfedge   h;
		int             val0, val1, val2, val3;
		int             val_opt0, val_opt1, val_opt2, val_opt3;
		int             ve0, ve1, ve2, ve3, ve_before, ve_after;
		bool            finished;
		int             i;


		// flip all edges
		for (finished = false, i = 0; !finished && i < 100; ++i)
		{
			finished = true;

			for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it)
			{
				if (!mesh_.is_boundary(*e_it))
				{
					// ------------- IMPLEMENT HERE ---------
					//  Extract valences of the four vertices involved to an eventual flip.
					//  Compute the sum of the squared valence deviances before flip
					//  Compute the sum of the squared valence deviances after an eventual flip
					//  If valence deviance is decreased and flip is possible, flip the vertex
					//  Leave the loop running until no collapse has been done (use the finished variable)
					// ------------- IMPLEMENT HERE ---------

					if (mesh_.is_flip_ok(*e_it)) {

						h = mesh_.halfedge(*e_it, 0);
						v2 = mesh_.to_vertex(mesh_.next_halfedge(h));
						//assert(h == mesh_.next_halfedge(mesh_.next_halfedge(mesh_.next_halfedge(h))));

						h = mesh_.halfedge(*e_it, 1);
						v3 = mesh_.to_vertex(mesh_.next_halfedge(h));
						//assert(h == mesh_.next_halfedge(mesh_.next_halfedge(mesh_.next_halfedge(h))));

						v0 = mesh_.vertex(*e_it, 0);
						v1 = mesh_.vertex(*e_it, 1);

						val0 = mesh_.valence(v0);
						val1 = mesh_.valence(v1);
						val2 = mesh_.valence(v2);
						val3 = mesh_.valence(v3);

						ve_before = (val0 - 6)*(val0 - 6) +
							(val1 - 6)*(val1 - 6) +
							(val2 - 6)*(val2 - 6) +
							(val3 - 6)*(val3 - 6);

						ve_after = (val0 - 7)*(val0 - 7) +
							(val1 - 7)*(val1 - 7) +
							(val2 - 5)*(val2 - 5) +
							(val3 - 5)*(val3 - 5);

						if (ve_after < ve_before) {
							finished = false;
							mesh_.flip(*e_it);
						}
					}
				}
			}
			//std::cout << "flip iterations: " << i << std::endl;
		}

		if (i == 100) std::cerr << "flip break\n";
	}

	void MeshProcessing::tangential_relaxation()
	{
		Mesh::Vertex_iterator     v_it, v_end(mesh_.vertices_end());
		Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
		int    valence;
		float dotProduct;
		Point     u, n;
		Point     laplace;
		int numVertices;
		int numBig100;
		int numSmall10;
		int totalVertice;

		Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
		Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");

		Mesh::Vertex_property<Scalar> v_unicurvature;

		// smooth
		for (int iters = 0; iters < 10; ++iters)
		{
			totalVertice = 0;
			numBig100 = 0;
			numSmall10 = 0;
			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
			{
				if (!mesh_.is_boundary(*v_it))
				{
					++totalVertice;
					// ------------- IMPLEMENT HERE ---------
					//  Compute uniform laplacian curvature approximation vector
					//  Compute the tangential component of the laplacian vector and move the vertex
					//  Store smoothed vertex location in the update vertex property.
					// ------------- IMPLEMENT HERE ---------

					// Initialize variables
					vv_c = mesh_.vertices(*v_it);
					if (!vv_c) {
						continue;
					}
					vv_end = vv_c;
					const Point& refPoint = mesh_.position(*v_it);
					laplace = 0.0f;
					numVertices = 0;

					// Iterate over adjacent vertices    
					do {
						++numVertices;
						const Point& vi = mesh_.position(*vv_c);
						laplace += vi - refPoint;
					} while (++vv_c != vv_end);

					// Average and normalize the Laplacian
					laplace /= numVertices;

					//std::cout << "lapl" << laplace << std::endl;
					//std::cout << "norm" << normals[*v_it] << std::endl;
					n = normals[*v_it];
					dotProduct = laplace.x * n.x + laplace.y * n.y + laplace.z * n.z;
					//std::cout << "u" << u << std::endl;

					n = n * dotProduct / norm(n);
					laplace = laplace - n;
					if (norm(laplace) > 100) {
						++numBig100;
					}
					else if (norm(laplace) < 10) {
						++numSmall10;
					}

					//std::cout << "refoint" << refPoint << std::endl;
					//std::cout << "laplace: " << norm(laplace) << std::endl;
					update[*v_it] = laplace;
					//break;
				}

			}
			//std::cout << "smooth iterations: " << iters << std::endl;
			//std::cout << "big100: " << numBig100 << std::endl;
			//std::cout << "smol10: " << numSmall10 << std::endl;
			//std::cout << "total: " << norm(laplace) << std::endl;

			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
				if (!mesh_.is_boundary(*v_it))
					mesh_.position(*v_it) += update[*v_it];
		}
	}

	Point MeshProcessing::calculateUniformDiscreteLaplacian(Mesh::Vertex v) {

		// Initialize variables
		Point acc_laplace = Point(0.0f);
		Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
		if (!vh_c) {
			return Point(0.0f);
		}
		Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
		const Point& refPoint = mesh_.position(v);
		int num_vertices = 0;
		bool hasBoundaryEdge = false;

		do {
			// Increment number of vertices
			num_vertices++;

			Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
			Mesh::Edge e = mesh_.edge(*vh_c);

			// Check for boundary
			if (mesh_.is_boundary(e)) {
				hasBoundaryEdge = true;
				break;
			}

			const Point& vi = mesh_.position(neighbor_v);
			acc_laplace += (vi - refPoint);

		} while (++vh_c != vh_end);

		if (hasBoundaryEdge) {
			// Curvature is 0 on boundary
			return Point(0.0f);
		}
		else {
			acc_laplace = acc_laplace / num_vertices;
			return acc_laplace;
		}
	}

	Point MeshProcessing::calculateCotanDiscreteLaplacian(Mesh::Vertex v, bool norm_total_weights) {

		Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0);

		// Initialize variables
		Point acc_laplace(0.0);
		Scalar total_weights = 0.f;
		Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
		if (!vh_c) {
			return Point(0.0f);
		}
		Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
		const Point& refPoint = mesh_.position(v);
		bool hasBoundaryEdge = false;

		// Iterate over adjacent vertices
		do {
			Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
			Mesh::Edge e = mesh_.edge(*vh_c);

			// Check for boundary
			if (mesh_.is_boundary(e)) {
				hasBoundaryEdge = true;
				break;
			}

			const Point& vi = mesh_.position(neighbor_v);
			acc_laplace += e_weight[e] * (vi - refPoint) / 2.0;
			total_weights += e_weight[e];

		} while (++vh_c != vh_end);

		if (hasBoundaryEdge) {
			// Curvature is 0 on boundary
			return Point(0.0f);
		}
		else {
			if (norm_total_weights) {
				return acc_laplace / total_weights;
			}
			else {
				return acc_laplace;
			}
		}

	}

	void MeshProcessing::calc_uniform_mean_curvature() {

		Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {
			v_unicurvature[v] = norm(calculateUniformDiscreteLaplacian(v));
		}
	}

	void MeshProcessing::calc_mean_curvature() {

		Mesh::Vertex_property<Scalar>  v_curvature = mesh_.vertex_property<Scalar>("v:curvature", 0);

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {
			v_curvature[v] = norm(calculateCotanDiscreteLaplacian(v, false));
		}
	}

	void MeshProcessing::calc_gauss_curvature() {

		Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0);
		Mesh::Vertex_property<Scalar> v_weight = mesh_.vertex_property<Scalar>("v:weight", 0);

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {

			// Initialize variables
			Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
			if (!vh_c) {
				continue;
			}
			Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
			const Point& refPoint = mesh_.position(v);
			Scalar theta = 0.f;
			bool hasBoundaryEdge = false;

			// Get the first vertex
			Point neighbor_p_before = mesh_.position(mesh_.to_vertex(*vh_c));
			++vh_c;

			// Iterate over adjacent vertices
			do {

				// Check for boundary
				if (mesh_.is_boundary(*vh_c)) {
					hasBoundaryEdge = true;
					break;
				}

				// Get next vertex
				Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));

				// Normalize the two triange edges
				Point d0 = normalize(neighbor_p_before - refPoint);
				Point d1 = normalize(neighbor_p_after - refPoint);

				// Accumulate the angle
				theta += acos(min(0.99f, max(-0.99f, dot(d0, d1))));

				// The "after" vertex become the "before" vertex
				neighbor_p_before = neighbor_p_after;

			} while (++vh_c != vh_end);

			if (hasBoundaryEdge) {
				// Curvature is 0 on boundary
				v_gauss_curvature[v] = 0.0f;
			}
			else {
				// Get last angle
				Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));
				Point d0 = normalize(neighbor_p_before - refPoint);
				Point d1 = normalize(neighbor_p_after - refPoint);
				theta += acos(min(0.99f, max(-0.99f, dot(d0, d1))));

				// Normalize
				v_gauss_curvature[v] = (2 * M_PI - theta) * 2.0f * v_weight[v];
			}
		}

	}

	void MeshProcessing::calc_weights() {
		calc_edges_weights();
		calc_vertices_weights();
	}

	void MeshProcessing::calc_edges_weights() {
		auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
		auto points = mesh_.vertex_property<Point>("v:point");

		Mesh::Halfedge h0, h1, h2;
		Point p0, p1, p2, d0, d1;

		for (auto e : mesh_.edges())
		{
			e_weight[e] = 0.0;

			h0 = mesh_.halfedge(e, 0);
			p0 = points[mesh_.to_vertex(h0)];

			h1 = mesh_.halfedge(e, 1);
			p1 = points[mesh_.to_vertex(h1)];

			if (!mesh_.is_boundary(h0))
			{
				h2 = mesh_.next_halfedge(h0);
				p2 = points[mesh_.to_vertex(h2)];
				d0 = p0 - p2;
				d1 = p1 - p2;
				e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
			}

			if (!mesh_.is_boundary(h1))
			{
				h2 = mesh_.next_halfedge(h1);
				p2 = points[mesh_.to_vertex(h2)];
				d0 = p0 - p2;
				d1 = p1 - p2;
				e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
			}
		}
	}

	void MeshProcessing::calc_vertices_weights() {
		Mesh::Face_around_vertex_circulator vf_c, vf_end;
		Mesh::Vertex_around_face_circulator fv_c;
		Scalar area;
		auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

		for (auto v : mesh_.vertices()) {
			area = 0.0;
			vf_c = mesh_.faces(v);

			if (!vf_c) {
				continue;
			}

			vf_end = vf_c;

			do {
				fv_c = mesh_.vertices(*vf_c);

				const Point& P = mesh_.position(*fv_c);  ++fv_c;
				const Point& Q = mesh_.position(*fv_c);  ++fv_c;
				const Point& R = mesh_.position(*fv_c);

				area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

			} while (++vf_c != vf_end);

			v_weight[v] = 0.5 / area;
		}
	}

	void MeshProcessing::load_mesh(const string &filename) {
		if (!mesh_.read(filename)) {
			std::cerr << "Mesh not found, exiting." << std::endl;
			exit(-1);
		}

		cout << "Mesh " << filename << " loaded." << endl;
		cout << "# of vertices : " << mesh_.n_vertices() << endl;
		cout << "# of faces : " << mesh_.n_faces() << endl;
		cout << "# of edges : " << mesh_.n_edges() << endl;

		// Compute the center of the mesh
		mesh_center_ = Point(0.0f, 0.0f, 0.0f);
		for (auto v : mesh_.vertices()) {
			mesh_center_ += mesh_.position(v);
		}
		mesh_center_ /= mesh_.n_vertices();

		// Compute the maximum distance from all points in the mesh and the center
		dist_max_ = 0.0f;
		for (auto v : mesh_.vertices()) {
			if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
				dist_max_ = distance(mesh_center_, mesh_.position(v));
			}
		}

		compute_mesh_properties();

		// Store the original mesh, this might be useful for some computations
		mesh_init_ = mesh_;
	}

	void MeshProcessing::compute_mesh_properties() {
		Mesh::Vertex_property<Point> vertex_normal =
			mesh_.vertex_property<Point>("v:normal");
		mesh_.update_face_normals();
		mesh_.update_vertex_normals();
		Mesh::Vertex_property<Color> v_color_valence =
			mesh_.vertex_property<Color>("v:color_valence",
				Color(1.0f, 1.0f, 1.0f));
		Mesh::Vertex_property<Color> v_color_unicurvature =
			mesh_.vertex_property<Color>("v:color_unicurvature",
				Color(1.0f, 1.0f, 1.0f));
		Mesh::Vertex_property<Color> v_color_curvature =
			mesh_.vertex_property<Color>("v:color_curvature",
				Color(1.0f, 1.0f, 1.0f));
		Mesh::Vertex_property<Color> v_color_gaussian_curv =
			mesh_.vertex_property<Color>("v:color_gaussian_curv",
				Color(1.0f, 1.0f, 1.0f));

		Mesh::Vertex_property<Scalar> vertex_valence =
			mesh_.vertex_property<Scalar>("v:valence", 0.0f);
		for (auto v : mesh_.vertices()) {
			vertex_valence[v] = mesh_.valence(v);
		}

		Mesh::Vertex_property<Scalar> v_unicurvature =
			mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
		Mesh::Vertex_property<Scalar> v_curvature =
			mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
		Mesh::Vertex_property<Scalar> v_gauss_curvature =
			mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

		calc_weights();
		calc_uniform_mean_curvature();
		calc_mean_curvature();
		calc_gauss_curvature();
		color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
			8 /* max */);
		color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
		color_coding(v_curvature, &mesh_, v_color_curvature);
		color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

		// get the mesh attributes and upload them to the GPU
		int j = 0;
		unsigned int n_vertices(mesh_.n_vertices());

		// Create big matrices to send the data to the GPU with the required
		// format
		color_valence_ = Eigen::MatrixXf(3, n_vertices);
		color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
		color_curvature_ = Eigen::MatrixXf(3, n_vertices);
		color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
		normals_ = Eigen::MatrixXf(3, n_vertices);
		points_ = Eigen::MatrixXf(3, n_vertices);
		indices_ = MatrixXu(3, mesh_.n_faces());

		for (auto f : mesh_.faces()) {
			std::vector<float> vv(3);
			int k = 0;
			for (auto v : mesh_.vertices(f)) {
				vv[k] = v.idx();
				++k;
			}
			indices_.col(j) << vv[0], vv[1], vv[2];
			++j;
		}

		j = 0;
		for (auto v : mesh_.vertices()) {
			points_.col(j) << mesh_.position(v).x,
				mesh_.position(v).y,
				mesh_.position(v).z;

			normals_.col(j) << vertex_normal[v].x,
				vertex_normal[v].y,
				vertex_normal[v].z;

			color_valence_.col(j) << v_color_valence[v].x,
				v_color_valence[v].y,
				v_color_valence[v].z;

			color_unicurvature_.col(j) << v_color_unicurvature[v].x,
				v_color_unicurvature[v].y,
				v_color_unicurvature[v].z;

			color_curvature_.col(j) << v_color_curvature[v].x,
				v_color_curvature[v].y,
				v_color_curvature[v].z;

			color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
				v_color_gaussian_curv[v].y,
				v_color_gaussian_curv[v].z;
			++j;
		}
	}

	void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
		Mesh::Vertex_property<Color> color_prop, Scalar min_value,
		Scalar max_value, int bound) {
		// Get the value array
		std::vector<Scalar> values = prop.vector();

		if (min_value == 0.0 && max_value == 0.0) {
			// discard upper and lower bound
			unsigned int n = values.size() - 1;
			unsigned int i = n / bound;
			std::sort(values.begin(), values.end());
			min_value = values[i];
			max_value = values[n - 1 - i];
		}

		// map values to colors
		for (auto v : mesh->vertices())
		{
			set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
		}
	}

	void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
		Mesh::Vertex_property<Color> color_prop)
	{
		color_prop[v] = col;
	}

	Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
		Scalar v0, v1, v2, v3, v4;
		v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
		v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
		v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
		v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
		v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

		Color col(1.0f, 1.0f, 1.0f);

		if (value < v0) {
			col = Color(0, 0, 1);
		}
		else if (value > v4) {
			col = Color(1, 0, 0);
		}
		else if (value <= v2) {
			if (value <= v1) { // [v0, v1]
				Scalar u = (value - v0) / (v1 - v0);
				col = Color(0, u, 1);
			}
			else { // ]v1, v2]
				Scalar u = (value - v1) / (v2 - v1);
				col = Color(0, 1, 1 - u);
			}
		}
		else {
			if (value <= v3) { // ]v2, v3]
				Scalar u = (value - v2) / (v3 - v2);
				col = Color(u, 1, 0);
			}
			else { // ]v3, v4]
				Scalar u = (value - v3) / (v4 - v3);
				col = Color(1, 1 - u, 0);
			}
		}
		return col;
	}


}


