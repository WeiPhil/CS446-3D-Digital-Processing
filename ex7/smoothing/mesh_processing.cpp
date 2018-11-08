//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#include "mesh_processing.h"
#include <set>

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

Point MeshProcessing::calclateUniformDiscreteLaplacian(Mesh::Vertex v){

    // Initialize variables
    Point acc_laplace = Point(0.0f);
    Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
    if(!vh_c) {
        return Point(0.0f);
    }
    Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
    const Point& refPoint = mesh_.position(v);
    int num_vertices = 0;
    bool hasBoundaryEdge = false;

    do {
        // Increment number of vertices
        num_vertices ++;
        
        Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
        Mesh::Edge e = mesh_.edge(*vh_c);

        // Check for boundary
        if(mesh_.is_boundary(e)){ 
            hasBoundaryEdge = true;
            break;
        }

        const Point& vi = mesh_.position(neighbor_v);
        acc_laplace += (vi-refPoint);

    } while(++vh_c != vh_end);

    if(hasBoundaryEdge){
        // Curvature is 0 on boundary
        return Point(0.0f);
    } else {
        acc_laplace = acc_laplace/num_vertices;
        return acc_laplace;
    }
}

Point MeshProcessing::calculateCotanDiscreteLaplacian(Mesh::Vertex v, bool norm_total_weights){

    Mesh::Edge_property<Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight", 0);

    // Initialize variables
    Point acc_laplace(0.0); 
    Scalar total_weights = 0.f;
    Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
    if(!vh_c) {
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
        if(mesh_.is_boundary(e)){ 
            hasBoundaryEdge = true;
            break;
        }

        const Point& vi = mesh_.position(neighbor_v);
        acc_laplace += e_weight[e] * (vi-refPoint);
        total_weights += e_weight[e];

    } while(++vh_c != vh_end);

    if(hasBoundaryEdge){
        // Curvature is 0 on boundary
        return Point(0.0f);
    } else {
        if (norm_total_weights) {
            return acc_laplace / total_weights;
        } else {
            return acc_laplace;
        }
    }
        
}
// ======================================================================
// EXERCISE 1.1
// ========================================================================
void MeshProcessing::uniform_smooth(const unsigned int iterations) {

    Scalar delta_t_lambda = 0.2f; // this is arbitrary and could be changed

    // We need to implement p'i = pi + delta_t_lambda * L(pi)
    // Where L(pi) = 1 / sum(wj) * sum(wj * (pj - pi) )

    for (unsigned int iter=0; iter<iterations; ++iter) {

         // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        // ------------- IMPLEMENT HERE ---------    

        // NON BOUNDARY CASE NOT HANDLED!

        for(auto v : mesh_.vertices()){

            const Point& refPoint = mesh_.position(v);

            mesh_.position(v) = refPoint + delta_t_lambda * calculateUniformDiscreteLaplacian(v);

        }
           

    }
}

// ======================================================================
// EXERCISE 1.2
// ========================================================================
void MeshProcessing::smooth(const unsigned int iterations) {

    Scalar delta_t_lambda = 0.5f; // this is arbitrary and could be changed

    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // Perform Cotan Laplacian smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        // 2) for each non-boundary vertex, update its position using the normalized cotan Laplacian operator
        //    (Hint: use the precomputed edge weights in the edge property "e:weight")
        // ------------- IMPLEMENT HERE ---------


         // NON BOUNDARY CASE NOT HANDLED!!

        calc_edges_weights(); // I suppose we must recall this so that it updates

        for(auto v : mesh_.vertices()){

            const Point& refPoint = mesh_.position(v);

            mesh_.position(v) = refPoint + delta_t_lambda * calculateCotanDiscreteLaplacian(v, true);

        }

    }
}

// ======================================================================
// EXERCISE 2
// ========================================================================
void MeshProcessing::implicit_smoothing(const double timestep) {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd B(n,3);

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets;

    Scalar diffusion = 1.37f/(M_PI*M_PI);
    Scalar delta_t_lambda = timestep * diffusion;

    // Calculate D then D^-1 
    // Calculate M
    for (auto v: mesh_.vertices()){

        // Right part
        Scalar area = 1.0f/area_inv[v];

        B(v.idx(),0) = area * points[v].x;
        B(v.idx(),1) = area * points[v].y;
        B(v.idx(),2) = area * points[v].z;

        // Left part

        Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
        Mesh::Vertex neighbor_v;
        Mesh::Edge e;

        Scalar total_cotan_weight = 0.0f;

        vh_c = mesh_.halfedges(v);

        vh_end = vh_c;

        do {
            neighbor_v = mesh_.to_vertex(*vh_c);
            e = mesh_.edge(*vh_c);

            total_cotan_weight += cotan[e];
            
            triplets.push_back(Eigen::Triplet<double>(v.idx(),neighbor_v.idx(),- delta_t_lambda * cotan[e]));

        } while(++vh_c != vh_end);

        triplets.push_back(Eigen::Triplet<double>(v.idx(),v.idx(),area+ delta_t_lambda * total_cotan_weight));

        
    }
    
    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 2 HERE
    // ========================================================================

    // build sparse matrix from triplets
    A.setFromTriplets(triplets.begin(), triplets.end());


    // solve A*X = B
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(A);
    Eigen::MatrixXd X = solver.solve(B);

    // copy solution
    for (int i = 0; i < n; ++i)
    {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = X(i, dim);
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

// ======================================================================
// EXERCISE 3.1
// ========================================================================
void MeshProcessing::uniform_laplacian_enhance_feature(const unsigned int iterations,
                                                       const unsigned int coefficient) {

    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the uniform Laplacian operator:
    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------

    Mesh::Vertex_property<Point> points_in = mesh_.add_vertex_property<Point>("points_in");
    Mesh::Vertex_property<Point> points_out = mesh_.add_vertex_property<Point>("points_out");

    for(auto v : mesh_.vertices()){
        points_in[v] = mesh_.position(v);

    }

    uniform_smooth(iterations);

    for(auto v : mesh_.vertices()){
        points_out[v] = mesh_.position(v);

    }

    Eigen::MatrixXf pout = *get_points();

    int n = mesh_.n_vertices();

    auto points = mesh_.vertex_property<Point>("v:point");

    for(auto v : mesh_.vertices()){ 
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = points_in[v][dim] + coefficient * (points_in[v][dim] - points_out[v][dim]);

    }

    mesh_.remove_vertex_property(points_in);
    mesh_.remove_vertex_property(points_out);
    
}

// ======================================================================
// EXERCISE 3.2
// ========================================================================
void MeshProcessing::cotan_laplacian_enhance_feature(const unsigned int iterations,
                                                      const unsigned int coefficient) {

    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the normalized cotan Laplacian operator:
    // 1) perform cotan Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------
    Mesh::Vertex_property<Point> points_in = mesh_.add_vertex_property<Point>("points_in");
    Mesh::Vertex_property<Point> points_out = mesh_.add_vertex_property<Point>("points_out");

    for(auto v : mesh_.vertices()){
        points_in[v] = mesh_.position(v);

    }

    smooth(iterations);

    for(auto v : mesh_.vertices()){
        points_out[v] = mesh_.position(v);

    }

    int n = mesh_.n_vertices();

    auto points = mesh_.vertex_property<Point>("v:point");

    for(auto v : mesh_.vertices()){
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = points_in[v][dim] + coefficient * (points_in[v][dim] - points_out[v][dim]);

    }

    mesh_.remove_vertex_property(points_in);
    mesh_.remove_vertex_property(points_out);


}

void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_uniform_mean_curvature() {

    Mesh::Vertex_property<Scalar> v_unicurvature = mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);

    // Iterate over vertices
    for (auto v: mesh_.vertices()){
        v_unicurvature[v] = norm(calculateUniformDiscreteLaplacian(v));
    }
}

void MeshProcessing::calc_mean_curvature() {
    
    Mesh::Vertex_property<Scalar>  v_curvature = mesh_.vertex_property<Scalar>("v:curvature", 0);

    // Iterate over vertices
    for (auto v: mesh_.vertices()){
        v_curvature[v] = norm(calculateCotanDiscreteLaplacian(v, false));
    }
}

void MeshProcessing::calc_gauss_curvature() {
    
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh_.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh_.vertex_property<Scalar>("v:weight", 0);

    // Iterate over vertices
	for (auto v: mesh_.vertices()){

        // Initialize variables
        Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
        if(!vh_c) {
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
            if(mesh_.is_boundary(*vh_c)){
                hasBoundaryEdge = true;
                break;
            }

            // Get next vertex
            Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));
            
            // Normalize the two triange edges
            Point d0 = normalize(neighbor_p_before - refPoint);
            Point d1 = normalize(neighbor_p_after - refPoint);

            // Accumulate the angle
            theta += acos(min(0.99f, max(-0.99f, dot(d0, d1)))) ;

            // The "after" vertex become the "before" vertex
            neighbor_p_before = neighbor_p_after;

        } while(++vh_c != vh_end);

        if(hasBoundaryEdge){
            // Curvature is 0 on boundary
            v_gauss_curvature[v] = 0.0f;
        } else {
            // Get last angle
            Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));
            Point d0 = normalize(neighbor_p_before - refPoint);
            Point d1 = normalize(neighbor_p_after - refPoint);
            theta += acos(min(0.99f, max(-0.99f, dot(d0, d1)))) ;

            // Normalize
            v_gauss_curvature[v] = (2*M_PI - theta) * 2.0f * v_weight[v];
        }
    }
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges())
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
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
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
    for (auto v: mesh_.vertices()) {
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
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
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

    for(auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
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
                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
               Mesh::Vertex_property<Color> color_prop) {
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}

MeshProcessing::~MeshProcessing() {}
}
