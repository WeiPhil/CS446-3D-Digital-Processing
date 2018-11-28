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
#include <Eigen/Dense>
#include <set>
#include <unordered_map>
#include <cmath>
namespace mesh_processing
{

using surface_mesh::Color;
using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Vec2d;
using surface_mesh::Vec3d;

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pair;
using std::unordered_map;
using std::vector;

MeshProcessing::MeshProcessing(const string &filename)
{
    load_mesh(filename);
}

// ======================================================================
// EXERCISE Begin
// ======================================================================

// ======================================================================
// EXERCISE 1.1 Mapping the Surface Boundary to 2D Circle
// ========================================================================
void MeshProcessing::map_suface_boundary_to_circle()
{

    Mesh::Vertex_property<Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.0));
    int n_vertices = mesh_.n_vertices();

    //Homework starting from here

    // Accumulate total boundary edge length here
    double acc_length = 0.f;
    std::vector<double> boundary_edge_lengths;
    bool has_mapped_boundary = false;
    Vec2d circle_center{0.5, 0.5};

    for (auto v : mesh_.vertices())
    {

        if (mesh_.is_boundary(v))
        { // Found a boundary vertex
            if (!has_mapped_boundary)
            { // Execute this part only once

                // starting_halfedge is guaranteed to be a boundary halfedge
                Mesh::Halfedge starting_halfedge = mesh_.halfedge(v);
                Mesh::Halfedge current_halfedge = starting_halfedge;

                // Now iterates over the surface boundary
                do
                {
                    // Get current edge and its endpoints
                    Mesh::Edge current_edge = mesh_.edge(current_halfedge);
                    Mesh::Vertex start_point = mesh_.vertex(current_edge, 0);
                    Mesh::Vertex end_point = mesh_.vertex(current_edge, 1);

                    // Accumulate edge length
                    Point pos = mesh_.position(end_point) - mesh_.position(start_point);
                    acc_length += sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
                    boundary_edge_lengths.push_back(acc_length);

                    // Get the next halfedge
                    current_halfedge = mesh_.halfedge(end_point);

                } while (current_halfedge != starting_halfedge);

                // We now have the total boundary length, re-iterate over the boundary
                size_t i = 0;
                do
                {
                    // Get current edge and its endpoints
                    Mesh::Edge current_edge = mesh_.edge(current_halfedge);
                    Mesh::Vertex start_point = mesh_.vertex(current_edge, 0);
                    Mesh::Vertex end_point = mesh_.vertex(current_edge, 1);

                    // angle is in [0, 2 * PI] (radian angle on circle)
                    double angle = (2 * EIGEN_PI) * (boundary_edge_lengths[i] / acc_length);

                    // Get corresponding point on cirle
                    v_texture[start_point] = circle_center + 0.5 * Vec2d{cos(angle), sin(angle)};

                    // Get the next halfedge
                    current_halfedge = mesh_.halfedge(end_point);
                    i++;

                } while (current_halfedge != starting_halfedge);

                // Map the boundary once
                has_mapped_boundary = true;
            }
        }
        else
        {
            // Every interior vertix goes to circle center
            v_texture[v] = circle_center;
        }
    }

    //Homework stopping from here
    //Update the texture matrix
    texture_ = Eigen::MatrixXf(2, n_vertices);
    int j = 0;
    for (auto v : mesh_.vertices())
    {
        texture_.col(j) << v_texture[v][0], v_texture[v][1];
        j++;
    }
}

// ======================================================================
// EXERCISE 1.2 Iterative Solve Textures
// ========================================================================
void MeshProcessing::iterative_solve_textures(int item_times)
{
    int n_vertices = mesh_.n_vertices();
    Mesh::Vertex_property<Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.5, 0.5));

    //Homework starting from here

    Mesh::Edge_property<Scalar> e_weight =
        mesh_.edge_property<Scalar>("e:weight", 0.0f);

    Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;

    for (auto v : mesh_.vertices())
    {
        //for non-boundary vertices
        if (!mesh_.is_boundary(v))
        {

            vh_c = mesh_.halfedges(v);
            if (!vh_c)
            {
                continue;
            }
            vh_end = vh_c;

            std::vector<std::pair<Scalar, Vec2d>> cotan_weightsTex;

            Scalar tot_cotan_weight = 0.f;
            // we calculate the non-normalized cotan weight and store them with the corresponding neighboor texture
            do
            {
                Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
                Mesh::Edge e = mesh_.edge(*vh_c);
                tot_cotan_weight += e_weight[e];
                cotan_weightsTex.push_back(std::pair<Scalar, Vec2d>(e_weight[e], v_texture[neighbor_v]));

            } while (++vh_c != vh_end);

            vh_c = mesh_.halfedges(v);

            Scalar cotan_weight = 0.f;
            Vec2d newTexCoord(0.0f);
            // Now we noramlize the cotan weights and calculate the new texture coordinates
            for (auto pair : cotan_weightsTex)
            {
                Scalar normalizedWeight = pair.first / tot_cotan_weight;
                newTexCoord += normalizedWeight * pair.second;
                cotan_weight += normalizedWeight;
            }

            v_texture[v] = newTexCoord / cotan_weight;
        }
    }

    //Homework stopping from here
    //Update the texture matrix
    texture_ = Eigen::MatrixXf(2, n_vertices);
    int j = 0;
    for (auto v : mesh_.vertices())
    {
        texture_.col(j) << v_texture[v][0], v_texture[v][1];
        j++;
    }
}
// ======================================================================
// EXERCISE 1.3 Direct Solve Textures
// ========================================================================
void MeshProcessing::direct_solve_textures()
{
    Mesh::Vertex_property<Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.0));
    int n_vertices = mesh_.n_vertices();
    //Homework starting from here

	// setup
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

	// Variable to store vertices and their position in the matrix
    std::map<Mesh::Vertex, int> boundaries;
    std::map<Mesh::Vertex, int> non_boundaries;

	// Storing vertices and their positions in the matrix
    int n = 0, m = 0;
    for (auto v : mesh_.vertices())
    {
        if (!mesh_.is_boundary(v))
        {
            non_boundaries.emplace(v, n);
            n++;
        }
        else
        {
            boundaries.emplace(v, m);
            m++;
        }
    }

	// Setting up matrices to compute solution
    Eigen::SparseMatrix<double> minusM(n+m, n+m);
    Eigen::MatrixXd minusB(Eigen::MatrixXd::Zero(n+m, 2));
    std::vector<Eigen::Triplet<double>> triplets_minusM;

	// Filling matrices with boundary vertices
	for (auto v : boundaries) {
		minusB(v.second, 0) = v_texture[v.first][0];
		minusB(v.second, 1) = v_texture[v.first][1];
		triplets_minusM.push_back(Eigen::Triplet<double>(v.second, v.second, 1.0));
	}

	// Filling matrices with non boundary vetices
    for (auto v : non_boundaries)
    {
		// Storing cotan weights
        Vec2d zero_vector_boundary_condition(0.0f);
        Scalar total_cotan_weight = 0.0f;

        // Get circulator
        Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v.first);
        Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;

        // Iterate over neighbors for each non-boundary vertex
        do
        {
            Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
            Mesh::Edge e = mesh_.edge(*vh_c);
            // i != j , j belongs to neighborhood and vj is non-boundary
            if (!mesh_.is_boundary(neighbor_v))
            {
                total_cotan_weight += cotan[e];
                triplets_minusM.push_back(Eigen::Triplet<double>(v.second + m, non_boundaries.at(neighbor_v) + m, cotan[e]));
            }
            else // i != j , j belongs to neighborhood and vj is boundary
            {
				total_cotan_weight += cotan[e];
				triplets_minusM.push_back(Eigen::Triplet<double>(v.second + m, boundaries.at(neighbor_v), cotan[e]));
            }
        } while (++vh_c != vh_end);

        triplets_minusM.push_back(Eigen::Triplet<double>(v.second+m, v.second+m, -total_cotan_weight));

        // Setting - bi
        minusB(v.second + m, 0) = zero_vector_boundary_condition[0];
        minusB(v.second + m, 1) = zero_vector_boundary_condition[1];
    }
	// Setting up L
    minusM.setFromTriplets(triplets_minusM.begin(), triplets_minusM.end());

	// Solving for solution
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(minusM);
    if (solver.info() != Eigen::Success)
    {
        printf("linear solver init failed.\n");
    }
    Eigen::MatrixXd X = solver.solve(minusB);
    if (solver.info() != Eigen::Success)
    {
        printf("linear solver failed.\n");
    }

	// Copying solutions into the texture
    for (auto v : non_boundaries)
    {
        v_texture[v.first][0] = X(v.second + m, 0);
        v_texture[v.first][1] = X(v.second + m, 1);
    }

    //Homework stopping from here
    //Update the texture matrix
    texture_ = Eigen::MatrixXf(2, n_vertices);
    int j = 0;
    for (auto v : mesh_.vertices())
    {
        texture_.col(j) << v_texture[v][0], v_texture[v][1];
        j++;
    }
}

// ======================================================================
// EXERCISE 2 Minimal Surfaces
// ======================================================================
void MeshProcessing::minimal_surface()
{
}

// ======================================================================
// EXERCISE End
// ======================================================================

void MeshProcessing::calc_weights()
{
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_uniform_mean_curvature()
{
    Mesh::Vertex_property<Scalar> v_unicurvature =
        mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------

    for (auto v : mesh_.vertices())
    {

        // Initialize variables
        Point acc_laplace = Point(0.0f);
        Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
        if (!vh_c)
        {
            v_unicurvature[v] = 0.f;
            continue;
        }
        Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
        const Point &refPoint = mesh_.position(v);
        int num_vertices = 0;
        bool hasBoundaryEdge = false;

        do
        {
            // Increment number of vertices
            num_vertices++;

            Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
            Mesh::Edge e = mesh_.edge(*vh_c);

            // Check for boundary
            if (mesh_.is_boundary(e))
            {
                hasBoundaryEdge = true;
                break;
            }

            const Point &vi = mesh_.position(neighbor_v);
            acc_laplace += (vi - refPoint);

        } while (++vh_c != vh_end);

        if (hasBoundaryEdge)
        {
            v_unicurvature[v] = 0.f;
        }
        else
        {
            v_unicurvature[v] = norm(acc_laplace / num_vertices);
        }
    }
}

void MeshProcessing::calc_mean_curvature()
{
    Mesh::Vertex_property<Scalar> v_curvature =
        mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
        mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
        mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    //@TODO Our implementation doesn't use v_weight (the line is there originally), wtf ?

    for (auto v : mesh_.vertices())
    {

        // Initialize variables
        Point acc_laplace(0.0);
        Scalar total_weights = 0.f;
        Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
        if (!vh_c)
        {
            v_curvature[v] = 0.f;
            continue;
        }
        Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
        const Point &refPoint = mesh_.position(v);
        bool hasBoundaryEdge = false;

        // Iterate over adjacent vertices
        do
        {
            Mesh::Vertex neighbor_v = mesh_.to_vertex(*vh_c);
            Mesh::Edge e = mesh_.edge(*vh_c);

            // Check for boundary
            if (mesh_.is_boundary(e))
            {
                hasBoundaryEdge = true;
                break;
            }

            const Point &vi = mesh_.position(neighbor_v);
            acc_laplace += e_weight[e] * (vi - refPoint) / 2.0;
            total_weights += e_weight[e];

        } while (++vh_c != vh_end);

        // Update curvature
        if (hasBoundaryEdge)
        {
            v_curvature[v] = 0.f;
        }
        else
        {
            v_curvature[v] = norm(acc_laplace / total_weights);
        }
    }
}

void MeshProcessing::calc_gauss_curvature()
{
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
        mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
        mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------

    // Iterate over vertices
    for (auto v : mesh_.vertices())
    {

        // Initialize variables
        Mesh::Halfedge_around_vertex_circulator vh_c = mesh_.halfedges(v);
        if (!vh_c)
        {
            continue;
        }
        Mesh::Halfedge_around_vertex_circulator vh_end = vh_c;
        const Point &refPoint = mesh_.position(v);
        Scalar theta = 0.f;
        bool hasBoundaryEdge = false;

        // Get the first vertex
        Point neighbor_p_before = mesh_.position(mesh_.to_vertex(*vh_c));
        ++vh_c;

        // Iterate over adjacent vertices
        do
        {

            // Check for boundary
            if (mesh_.is_boundary(*vh_c))
            {
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

        if (hasBoundaryEdge)
        {
            // Curvature is 0 on boundary
            v_gauss_curvature[v] = 0.0f;
        }
        else
        {
            // Get last angle
            Point neighbor_p_after = mesh_.position(mesh_.to_vertex(*vh_c));
            Point d0 = normalize(neighbor_p_before - refPoint);
            Point d1 = normalize(neighbor_p_after - refPoint);
            theta += acos(min(0.99f, max(-0.99f, dot(d0, d1))));

            // Normalize
            v_gauss_curvature[v] = (2 * EIGEN_PI - theta) * 2.0f * v_weight[v];
        }
    }
}

void MeshProcessing::calc_edges_weights()
{
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e : mesh_.edges())
    {
        double w = 0;
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
            w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));
        }

        w = w < 0 ? 0 : w;
        e_weight[e] = w * 0.5;
    }
}

void MeshProcessing::calc_vertices_weights()
{
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v : mesh_.vertices())
    {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if (!vf_c)
        {
            continue;
        }

        vf_end = vf_c;

        do
        {
            fv_c = mesh_.vertices(*vf_c);

            const Point &P = mesh_.position(*fv_c);
            ++fv_c;
            const Point &Q = mesh_.position(*fv_c);
            ++fv_c;
            const Point &R = mesh_.position(*fv_c);

            area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

        } while (++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename)
{
    if (!mesh_.read(filename))
    {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh " << filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v : mesh_.vertices())
    {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v : mesh_.vertices())
    {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_)
        {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    init_textures();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::init_textures()
{
    Mesh::Vertex_property<Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0, 0));
    int n_vertices = mesh_.n_vertices();
    texture_ = Eigen::MatrixXf(2, n_vertices);
    int j = 0;

    double min[3] = {1e10, 1e10, 1e10};
    double max[3] = {-1e10, -1e10, -1e10};
    for (auto v : mesh_.vertices())
    {
        Point p = mesh_.position(v);
        for (int kd = 0; kd < 3; kd++)
        {
            if (p[kd] < min[kd])
                min[kd] = p[kd];
            if (p[kd] > max[kd])
                max[kd] = p[kd];
        }
    }
    for (auto v : mesh_.vertices())
    {
        Point p = mesh_.position(v);
        v_texture[v][0] = (p[0] - min[0]) / (max[0] - min[0]);
        v_texture[v][1] = (p[1] - min[1]) / (max[1] - min[1]);
        texture_.col(j) << v_texture[v][0],
            v_texture[v][1];
        j++;
    }
}

void MeshProcessing::compute_mesh_properties()
{
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
    for (auto v : mesh_.vertices())
    {
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

    for (auto f : mesh_.faces())
    {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v : mesh_.vertices(f))
        {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v : mesh_.vertices())
    {
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
                                  Mesh::Vertex_property<Color> color_prop, int bound)
{
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size() - 1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n - 1 - i];

    // map values to colors
    for (auto v : mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color &col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value)
{
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
    v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
    v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
    v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
    v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0)
    {
        col = Color(0, 0, 1);
    }
    else if (value > v4)
    {
        col = Color(1, 0, 0);
    }
    else if (value <= v2)
    {
        if (value <= v1)
        { // [v0, v1]
            Scalar u = (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        }
        else
        { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1 - u);
        }
    }
    else
    {
        if (value <= v3)
        { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        }
        else
        { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1 - u, 0);
        }
    }
    return col;
}

MeshProcessing::~MeshProcessing() {}
} // namespace mesh_processing
