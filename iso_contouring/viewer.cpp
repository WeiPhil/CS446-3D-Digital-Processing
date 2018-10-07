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

// ========================================================================
// EXERCISE 1
// ========================================================================
Scalar Viewer::iso_value(Point v_pos)
{
    float x,y;
    x = v_pos.x;
    y = v_pos.y;

    // ----- (un)comment a line to change the function you are testing
    //Scalar iso = sqrt(x*x + y*y) - 1;
    Scalar iso = sin(2*x+2*y) - cos(4*x*y) +1;
    //Scalar iso = y*y - sin(x*x);

    return iso;
}

/*	@brief	Computes the coordinate of an interpolated point on a given segment [AB]

	The implementation assumes that point A is taken for reference, meaning that the function
	will return A when p=0 and B when p=1. In our case A should always be the "outlier" point 
	in a given triangle, that is the point with a differently signed ISO value than the others.

	@param	a	A point (reference, beginning of segment)
	@param	b	Another point (end of segment)
	@param	p	The proportion characterizing the interpolation (in [0,1])

	@return	An interpolated point on segment [AB], at proportion p
*/
Point interpolatePoint(Point a, Point b, Scalar p) {
	return (1 - p) * a + p * b;
}

/*	@brief	Computes a simple proportion

	The implementation assumes that isoA is taken for reference, meaning that the function
	will return 0 isoA=0, whatever isoB's value. In our case isoA should always be the ISO 
	of the "outlier" point in a given triangle, that is the point with a differently signed 
	ISO value than the others.

	@param	isoA	An ISO value (reference)
	@parma	isoB	Another ISO value

	@return	The "proportion" of isoA in (isoA + isoB)
*/
Scalar computeProportion(Scalar isoA, Scalar isoB) {
	return abs(isoA) / (abs(isoA) + abs(isoB));
}

void Viewer::calc_iso_contouring() {
    Mesh::Vertex_property<Scalar> v_iso = mesh.vertex_property<Scalar>("v:iso", 0);
    segment_points.clear();
    std::vector<Point> v_positions(mesh.n_vertices());
    std::vector<std::vector<int> > triangle_ids;

    for (auto v: mesh.vertices()) {
        Point v_pos = mesh.position(v);
        v_positions[v.idx()] = v_pos;
        Scalar iso = 0;

        iso = iso_value(v_pos);

        v_iso[v] = iso; //this variable is for coloring the density; do not change this variable
    }

    for(auto f: mesh.faces()) {
        std::vector<int> vv(3);
        int k = 0;
        for (auto v: mesh.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        triangle_ids.push_back(vv);
    }

	// To store the ISO value for each vertex
    std::vector<Scalar> iso_values(mesh.n_vertices());

	// Fill in the iso_values array
    for(auto v: mesh.vertices())
        iso_values[v.idx()] = iso_value(v_positions[v.idx()]);
    
    for(auto triangle : triangle_ids){ // For each triangle

        // Get the 3 vertices
		Point vertex0 = v_positions[triangle[0]];
		Point vertex1 = v_positions[triangle[1]];
		Point vertex2 = v_positions[triangle[2]];

		// Get the 3 ISO values
        Scalar iso0 = iso_values[triangle[0]];
        Scalar iso1 = iso_values[triangle[1]];
        Scalar iso2 = iso_values[triangle[2]];

		if ((iso0 > 0 && iso1 > 0 && iso2 > 0) || (iso0 <= 0 && iso1 <= 0 && iso2 <= 0)) {
			// Nothing to do
		}
		else if ((iso0 <= 0 && iso1 > 0 && iso2 > 0) || (iso0 > 0 && iso1 <= 0 && iso2 <= 0)) {
			// Point 0 is the "outlier"
			segment_points.push_back(interpolatePoint(vertex0, vertex1, computeProportion(iso0, iso1)));
			segment_points.push_back(interpolatePoint(vertex0, vertex2, computeProportion(iso0, iso2)));
		}
		else if ((iso0 > 0 && iso1 <= 0 && iso2 > 0) || (iso0 <= 0 && iso1 > 0 && iso2 <= 0)) {
			// Point 1 is the "outlier"
			segment_points.push_back(interpolatePoint(vertex1, vertex0, computeProportion(iso1, iso0)));
			segment_points.push_back(interpolatePoint(vertex1, vertex2, computeProportion(iso1, iso2)));
		}
		else if ((iso0 > 0 && iso1 > 0 && iso2 <= 0) || (iso0 <= 0 && iso1 <= 0 && iso2 > 0)) {
			// Point 2 is the "outlier"
			segment_points.push_back(interpolatePoint(vertex2, vertex0, computeProportion(iso2, iso0)));
			segment_points.push_back(interpolatePoint(vertex2, vertex1, computeProportion(iso2, iso1)));
		}
		else {
			// Should never happen
		}
    }
}
