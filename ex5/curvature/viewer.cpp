//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

using std::min;
using std::max;
using namespace surface_mesh;

typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

void Viewer::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void Viewer::calc_edges_weights() {
    Mesh::Halfedge h0, h1, h2;
    Mesh::Vertex v0, v1;
    Point p0, p1, p2, d0, d1;
    Scalar w;
    auto eweight = mesh.edge_property<Scalar>("e:weight", 0);
    for (auto e: mesh.edges()) {
        w = 0.0;

        h0 = mesh.halfedge(e, 0);
        v0 = mesh.to_vertex(h0);
        p0 = mesh.position(v0);

        h1 = mesh.halfedge(e, 1);
        v1 = mesh.to_vertex(h1);
        p1 = mesh.position(v1);

        h2 = mesh.next_halfedge(h0);
        p2 = mesh.position(mesh.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(min(0.99f, max(-0.99f, dot(d0, d1)))));

        h2 = mesh.next_halfedge(h1);
        p2 = mesh.position(mesh.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(min(0.99f, max(-0.99f, dot(d0, d1)))));

        w = max(0.0f, w);
        eweight[e] = w * 0.5;
    }
}

void Viewer::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto vweight = mesh.vertex_property<Scalar>("v:weight", 0);

    for (auto v: mesh.vertices()) {
        area = 0.0;
        vf_c = mesh.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh.vertices(*vf_c);

            const Point& P = mesh.position(*fv_c);  ++fv_c;
            const Point& Q = mesh.position(*fv_c);  ++fv_c;
            const Point& R = mesh.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        vweight[v] = 0.5 / area;
    }
}

void Viewer::computeValence() {
   Mesh::Vertex_property<Scalar> vertex_valence =
            mesh.vertex_property<Scalar>("v:valence", 0);
    for (auto v: mesh.vertices()) {
        vertex_valence[v] = mesh.valence(v);
    }
}

std::tuple<Normal,Scalar> Viewer::computeNormalOfFaceAndArea(Mesh::Vertex_around_face_circulator fv_c){

    const Point& xi = mesh.position(*fv_c);  ++fv_c;
    const Point& xj = mesh.position(*fv_c);  ++fv_c;
    const Point& xk = mesh.position(*fv_c);

    Point crossProduct = cross(xj-xi,xk-xi);
    Scalar area = norm(crossProduct) * 0.5f;

    return std::make_tuple(crossProduct.normalize(), area);
}

std::tuple<Normal,Scalar> Viewer::computeNormalOfFaceAndAngle(Mesh::Vertex_around_face_circulator fv_c){

    const Point& xi = mesh.position(*fv_c);  ++fv_c;
    const Point& xj = mesh.position(*fv_c);  ++fv_c;
    const Point& xk = mesh.position(*fv_c);

    return std::make_tuple(cross(xj-xi,xk-xi).normalize(), std::acos(dot((xj-xi).normalize(), (xk-xi).normalize()) ) );
}

Normal Viewer::computeNormalOfFace(Mesh::Vertex_around_face_circulator fv_c){

    const Point& xi = mesh.position(*fv_c);  ++fv_c;
    const Point& xj = mesh.position(*fv_c);  ++fv_c;
    const Point& xk = mesh.position(*fv_c);

    return cross(xj-xi,xk-xi).normalize();
}
// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Viewer::computeNormalsWithConstantWeights() {
    Point default_normal(0.0, 1.0, 0.0);
	Point base_point(0.0, 0.0, 0.0);
    Mesh::Vertex_property<Point> v_cste_weights_n =
            mesh.vertex_property<Point>("v:cste_weights_n", default_normal);

    Mesh::Vertex_around_face_circulator fv_c;
    Mesh::Face_around_vertex_circulator vf_c, vf_end;

    for (auto v: mesh.vertices()){
        
        vf_c = mesh.faces(v);
        
        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        std::vector<Normal> face_normals;

        do {
            fv_c = mesh.vertices(*vf_c);

            face_normals.push_back(computeNormalOfFace(fv_c));

        } while(++vf_c != vf_end);

		Normal vertexNormal = base_point;
        for(auto n : face_normals){
            vertexNormal += n;
        }

        v_cste_weights_n[v] = vertexNormal.normalize();
       
    }

}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Viewer::computeNormalsByAreaWeights() {
    Point default_normal(0.0, 1.0, 0.0);
	Point base_point(0.0, 0.0, 0.0);
    Mesh::Vertex_property<Point> v_area_weights_n =
            mesh.vertex_property<Point>("v:area_weight_n", default_normal);

    Mesh::Vertex_around_face_circulator fv_c;
    Mesh::Face_around_vertex_circulator vf_c, vf_end;

    for (auto v: mesh.vertices()){
        
        vf_c = mesh.faces(v);
        
        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        std::vector<std::tuple<Normal,Scalar>> face_normals_area;

        do {
            fv_c = mesh.vertices(*vf_c);

            face_normals_area.push_back(computeNormalOfFaceAndArea(fv_c));

        } while(++vf_c != vf_end);

        Normal vertexNormal = base_point;
        
        for(auto n : face_normals_area){
            vertexNormal += (std::get<0>(n) * std::get<1>(n) );
        }
        v_area_weights_n[v] = vertexNormal.normalize(); 
       
    }
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Viewer::computeNormalsWithAngleWeights() {
    Point default_normal(0.0, 1.0, 0.0);
	Point base_point(0.0, 0.0, 0.0);
    Mesh::Vertex_property<Point> v_angle_weights_n =
            mesh.vertex_property<Point>("v:angle_weight_n", default_normal);

    Mesh::Vertex_around_face_circulator fv_c;
    Mesh::Face_around_vertex_circulator vf_c, vf_end;

    for (auto v: mesh.vertices()){
        
        vf_c = mesh.faces(v);
        
        if(!vf_c) {
            cout << "No faces?" << endl;
            continue;
        }

        vf_end = vf_c;

        std::vector<std::tuple<Normal,Scalar>> face_normals_angle;

        do {
            fv_c = mesh.vertices(*vf_c);

            face_normals_angle.push_back(computeNormalOfFaceAndAngle(fv_c));

        } while(++vf_c != vf_end);

		Normal vertexNormal = base_point;
        for(auto n : face_normals_angle){
            vertexNormal += std::get<0>(n)*std::get<1>(n);
        }

        v_angle_weights_n[v] = vertexNormal.normalize(); 
       
    }
    // ------------- IMPLEMENT HERE ---------
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the angles technique (see .pdf) and store it inside v_angle_weights_n[v]
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Viewer::calc_uniform_laplacian() {
    Mesh::Vertex_property<Scalar> v_uniLaplace = mesh.vertex_property<Scalar>("v:uniLaplace", 0);
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    Point laplace(0.0);
    min_uniLaplace = 1000;
    max_uniLaplace = -1000;

    // Need to test if a vertex is a boundary vertex, if it is not we can calculate the uniform laplacian

    for (auto v: mesh.vertices()){
        
        vv_c = mesh.vertices(v);
        
        if(!vv_c) {
            continue;
        }

        vv_end = vv_c;

        const Point& refPoint = mesh.position(v);

        int numVertices = 0;
    
        do {
            ++ numVertices;
            const Point& vi = mesh.position(*vv_c);
            laplace += refPoint - vi ;
        } while(++vv_c != vv_end);

        laplace = laplace/numVertices;

        v_uniLaplace[v] = norm(laplace);

        if(min_uniLaplace > v_uniLaplace[v])
            min_uniLaplace = v_uniLaplace[v];

        if(max_uniLaplace < v_uniLaplace[v])
            max_uniLaplace = v_uniLaplace[v];

    }
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, compute uniform Laplacian operator vector
    // and store the vector length in the vertex property of the
    // mesh called v_uniLaplace[v].
    // Store min and max values of v_uniLaplace[v] in min_uniLaplace and max_uniLaplace.
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Viewer::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
    Mesh::Vertex_property<Scalar>  v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
    Mesh::Vertex neighbor_v;
    Mesh::Edge e;
    Point laplace(0.0);
    min_mean_curvature = 1000;
    max_mean_curvature = -1000;



     for (auto v: mesh.vertices()){
        
        vh_c = mesh.halfedges(v);
        
        if(!vh_c) {
            continue;
        }

        vh_end = vh_c;

        const Point& refPoint = mesh.position(v);
    
        do {
            neighbor_v = mesh.to_vertex(*vh_c);
            e = mesh.edge(*vh_c);

            const Point& vi = mesh.position(neighbor_v);
            laplace += e_weight[e] * (refPoint - vi);
        } while(++vh_c != vh_c);

        laplace = laplace/v_weight[v];

        v_curvature[v] = norm(laplace);

        if(min_mean_curvature > v_curvature[v])
            min_mean_curvature = v_curvature[v];

        if(max_mean_curvature < v_curvature[v])
            max_mean_curvature = v_curvature[v];

    }
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // ------------- IMPLEMENT HERE ---------
}
// ========================================================================
// EXERCISE 2.3
// ========================================================================
void Viewer::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    Mesh::Vertex_around_vertex_circulator vv_c, vv_c2, vv_end;
    Point d0, d1;
    Scalar angles, cos_angle;
    Scalar lb(-1.0f), ub(1.0f);
    min_gauss_curvature = 1000;
    max_gauss_curvature = -1000;

	for (auto v : mesh.vertices()) {

		vv_c = mesh.vertices(v);
		vv_c2 = mesh.vertices(v);

		if (!vv_c) {
			continue;
		}
		if (!vv_c2) {
			continue;
		}

		vv_end = vv_c;
		++vv_c2;
		cos_angle = 0.f;

		const Point& refPoint = mesh.position(v);

		do {
			const Point& vi = mesh.position(*vv_c);
			const Point& vi2 = mesh.position(*vv_c2);
			d0 = norm(refPoint - vi);
			d1 = norm(refPoint - vi2);
			cos_angle += acos(min(ub, max(lb, dot(d0, d1))));

			++vv_c;
		} while (++vv_c2 != vv_end);

		angles =  cos_angle;

		v_gauss_curvature[v] = (2*M_PI - angles) / v_weight[v];

		if (min_gauss_curvature > v_gauss_curvature[v])
			min_gauss_curvature = v_gauss_curvature[v];

		if (max_gauss_curvature < v_gauss_curvature[v])
			max_gauss_curvature = v_gauss_curvature[v];

	}


    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // ------------- IMPLEMENT HERE ---------
}
