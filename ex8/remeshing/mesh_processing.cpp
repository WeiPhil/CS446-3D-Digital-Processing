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
					target_new_length[*v_it] = target_length[*v_it] + laplace;
				}
				
			}

			// rescale desired length
			Scalar totalLength = 0;

			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
				totalLength += target_new_length[*v_it];
			}

			std::cout << "total length: " << totalLength << std::endl;

			Scalar avrgLength = totalLength / mesh_.n_vertices();

			std::cout << "avrg length: " << avrgLength << std::endl;
			std::cout << "tagret length: " << TARGET_LENGTH << std::endl;

			for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
				target_new_length[*v_it] = target_new_length[*v_it] * (TARGET_LENGTH / avrgLength);
			}
		}
	}

	void MeshProcessing::split_long_edges()
	{
		Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
		Mesh::Vertex   v0, v1, v;
		bool            finished;
		int             i;
		Point p0, p1, p;
		float avrgLength, avrgDesired;
		int nmbrPoints, nmbrSplits, nmbrEgdeUnder;

		Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);


		for (finished = false, i = 0; !finished && i < 100; ++i)
		{
			finished = true;
			// ------------- IMPLEMENT HERE ---------
			// INSERT CODE:
			//  Compute the desired length as the mean between the property target_length of two vertices of the edge
			//  If the edge is longer than 4/3 * desired length
			//		add the midpoint to the mesh
			//		set the interpolated normal and interpolated vtargetlength_ property to the vertex
			//		split the edge with this vertex (use openMesh function split)
			// Leave the loop running until no splits are done (use the finished variable)
			// ------------- IMPLEMENT HERE ---------

			Scalar edgeLength;
			Scalar targetLength;
			avrgDesired = 0;
			avrgLength = 0;
			nmbrPoints = 0;
			nmbrSplits = 0;
			nmbrEgdeUnder = 0;

			for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it)
			{


				edgeLength = mesh_.edge_length(*e_it);
				v0 = mesh_.vertex(*e_it, 0);
				v1 = mesh_.vertex(*e_it, 1);
				targetLength = (target_length[v0] + target_length[v1]) / 2;

				
					avrgLength += edgeLength;
					avrgDesired += targetLength;
					++nmbrPoints;
					if (edgeLength < (04.0f / 3.0f)) {
						++nmbrEgdeUnder;
					}
		

				if (edgeLength > ((4.f / 3.f) * targetLength)) {
					//std::cout << edgeLength << " > " << targetLength << std::endl;
					++nmbrSplits;

					finished = false;

					p0 = mesh_.position(v0);
					p1 = mesh_.position(v1);
					p = (p0 + p1) / 2;

					v = mesh_.add_vertex(p);
					
					//std::cout << "added vertex" << std::endl;

					normals[v] = (normals[v0] + normals[v1]) / 2;
					target_length[v] = (target_length[v0] + target_length[v1]) / 2;

					//std::cout << "computed normals" << std::endl;

					//mesh_.insert_vertex(*e_it, v);
					//std::cout << "insert" << std::endl;
					mesh_.split(*e_it, v);
					//std::cout << "split" << std::endl;
				}

			}
			//std::cout << "split iterations: " << i << std::endl;
			//std::cout << "avrg length: " << avrgLength/nmbrPoints << std::endl;
			//std::cout << "avrg desired: " << avrgDesired/nmbrPoints << std::endl;
			//std::cout << "avrg cutoff: " << (4.f/3.f)*avrgDesired / nmbrPoints << std::endl;
			//std::cout << "nmbr splits: " << nmbrSplits << std::endl;
			//std::cout << "nmbr egde under 4/3: " << nmbrEgdeUnder << std::endl;
			//std::cout << "nmbr egdes: " << nmbrPoints << std::endl;

		}
	}

	void MeshProcessing::collapse_short_edges()
	{
		Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
		Mesh::Vertex   v0, v1;
		Mesh::Halfedge  h01, h10;
		bool            finished, b0, b1;
		int             i;
		bool            hcol01, hcol10;

		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

		Scalar edgeLength;
		Scalar targetLength;
		//100 de base
		for (finished = false, i = 0; !finished && i < 100; ++i)
		{
			finished = true;

			for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it)
			{
				if (!mesh_.is_deleted(*e_it)) // might already be deleted
				{
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

					edgeLength = mesh_.edge_length(*e_it);
					v0 = mesh_.vertex(*e_it, 0);
					v1 = mesh_.vertex(*e_it, 1);
					targetLength = (target_length[v0] + target_length[v1]) / 2;

					//std::cout << edgeLength << " < " << targetLength << std::endl;

					if (edgeLength < ((4.f / 5.f) * targetLength)) {
						if (!mesh_.is_boundary(*e_it)) {
							h01 = mesh_.halfedge(*e_it, 0);
							h10 = mesh_.halfedge(*e_it, 1);
							b0 = mesh_.is_collapse_ok(h01);
							b1 = mesh_.is_collapse_ok(h10);
							if (b0 && b1) {
								finished = false;
								if (mesh_.valence(v0) > mesh_.valence(v1)) {
									mesh_.collapse(h01);
								}else {
									mesh_.collapse(h10);
								}
							}else if(b0) {
								finished = false;
								mesh_.collapse(h01);
							}
							else if (b1) {
								finished = false;
								mesh_.collapse(h10);
							}
						}
					}
				}
			}
			//std::cout << "collapse iterations: " << i << std::endl;
		}

		mesh_.garbage_collection();

		if (i == 100) std::cerr << "collapse break\n";
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

					n = n * dotProduct/norm(n);
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

	void MeshProcessing::calc_uniform_mean_curvature() {
		Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);

		// ------------- IMPLEMENT HERE ---------
		// For each non-boundary vertex, approximate mean curvature using
		// the length of the uniform Laplacian approximation
		// Save your approximation in unicurvature vertex property of the mesh.
		// ------------- IMPLEMENT HERE ---------

		Mesh::Vertex_around_vertex_circulator vv_c, vv_end;

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {

			// Initialize variables
			Point acc_laplace(0.0);
			vv_c = mesh_.vertices(v);
			if (!vv_c) {
				continue;
			}
			vv_end = vv_c;
			const Point& refPoint = mesh_.position(v);
			int numVertices = 0;

			// Iterate over adjacent vertices    
			do {
				++numVertices;
				const Point& vi = mesh_.position(*vv_c);
				acc_laplace += vi - refPoint;
			} while (++vv_c != vv_end);

			// Average and normalize the Laplacian
			acc_laplace /= numVertices;
			v_unicurvature[v] = norm(acc_laplace);
		}
	}

	void MeshProcessing::calc_mean_curvature() {
		Mesh::Vertex_property<Scalar>  v_curvature = mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
		Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
		Mesh::Vertex_property<Scalar>  v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

		// ------------- IMPLEMENT HERE ---------
		// For all non-boundary vertices, approximate the mean curvature using
		// the length of the Laplace-Beltrami approximation.
		// Save your approximation in v_curvature vertex property of the mesh.
		// Use the weights from calc_weights(): e_weight and v_weight
		// ------------- IMPLEMENT HERE ---------

		Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {

			// Initialize variables
			Point acc_laplace(0.0);
			vh_c = mesh_.halfedges(v);
			if (!vh_c) {
				continue;
			}
			vh_end = vh_c;
			const Point& refPoint = mesh_.position(v);

			// Iterate over adjacent vertices
			do {
				Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
				Mesh::Edge e = mesh_.edge(*vh_c);
				const Point& vi = mesh_.position(neighbor_v);
				acc_laplace += e_weight[e] * 2.0f * (vi - refPoint);

			} while (++vh_c != vh_end);

			// Multiply by vertex's weight and normalize
			acc_laplace *= v_weight[v];
			v_curvature[v] = norm(acc_laplace);
		}

	}

	void MeshProcessing::calc_gauss_curvature() {
		Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
		Mesh::Vertex_property<Scalar> v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

		// ------------- IMPLEMENT HERE ---------
		// For each non-boundary vertex, approximate Gaussian curvature,
		// and store it in the vertex property v_gauss_curvature.
		// Hint: When calculating angles out of cross products make sure the value
		// you pass to the acos function is between -1.0 and 1.0.
		// Use the v_weight property for the area weight.
		// ------------- IMPLEMENT HERE ---------

		Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;

		// Iterate over vertices
		for (auto v : mesh_.vertices()) {

			// Initialize variables
			vh_c = mesh_.halfedges(v);
			if (!vh_c) {
				continue;
			}
			vh_end = vh_c;
			const Point& refPoint = mesh_.position(v);
			Scalar theta = 0.f;

			// Get the first vertex
			Point neighbor_p_before = mesh_.position(mesh_.to_vertex(*vh_c));
			++vh_c;

			// Iterate over adjacent vertices
			do {

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

			// Get last angle
			Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));
			Point d0 = normalize(neighbor_p_before - refPoint);
			Point d1 = normalize(neighbor_p_after - refPoint);
			theta += acos(min(0.99f, max(-0.99f, dot(d0, d1))));

			// Normalize
			v_gauss_curvature[v] = (2 * M_PI - theta) * 2.0f * v_weight[v];
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


