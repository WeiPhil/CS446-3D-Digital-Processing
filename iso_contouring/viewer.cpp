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
    //Scalar iso = sin(2*x+2*y) - cos(4*x*y) +1;
    Scalar iso = y*y - sin(x*x);

    return iso;
}

Scalar lerp(Scalar A,Scalar B,Scalar a,Scalar b){
    return  ( -A * (b-a)/(B-A) ) - a; 
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

    //segment_points is defined in viewer.h as std::vector<Point> segment_points;
    //add points in segment_points forming an edge one after the other;
    //for example segment_points[0] and segment_points[1] are two points forming the first edge
    //and segment_points[2] and segment_points[3] are two points forming the second edge



    /*
    Classify grid nodes as inside/outside
        Is below or above iso-value?
    Classify cell: configura ons
        in/out for each corner
    Determine contour edges
        look-up table for edge configura on
    Determine vertex posi ons
        linear interpola on of grid values along edges

    */
    // ----- add your code here -----
    std::vector<Scalar> iso_values(mesh.n_vertices());

    for(auto v: mesh.vertices())
        iso_values[v.idx()] = iso_value(v_positions[v.idx()]);
    
    for(auto triangle : triangle_ids){
        // triangle is a set of three indices representing the 3 vertices it composes
        Scalar iso_1 = iso_values[triangle[0]];
        Scalar iso_2 = iso_values[triangle[1]];
        Scalar iso_3 = iso_values[triangle[2]];
        if (!((iso_1 < 0 && iso_2 < 0  && iso_3 < 0 ) || (iso_1 >= 0 && iso_2 >= 0  && iso_3 >= 0))){

            if(iso_1 < 0 && iso_2 >= 0){ // case - + - and - + +
                if(iso_3 < 0){ //case - + -
                    Scalar xPos1,yPos1,xPos2,yPos2;
                    if(v_positions[triangle[0]].x < v_positions[triangle[1]].x){
                        xPos1 = lerp(iso_1,iso_2,v_positions[triangle[0]].x,v_positions[triangle[1]].x);
                    }else{
                        xPos1 = lerp(iso_1,iso_2,v_positions[triangle[1]].x,v_positions[triangle[0]].x);
                    }

                    if(v_positions[triangle[0]].y < v_positions[triangle[1]].y){
                        yPos1 = lerp(iso_1,iso_2,v_positions[triangle[0]].y,v_positions[triangle[1]].y);
                    }else{
                        yPos1 = lerp(iso_1,iso_2,v_positions[triangle[1]].y,v_positions[triangle[0]].y);
                    }

                    if(v_positions[triangle[1]].x < v_positions[triangle[2]].x){
                        xPos2 = lerp(iso_1,iso_2,v_positions[triangle[1]].x,v_positions[triangle[2]].x);
                    }else{
                        xPos2 = lerp(iso_1,iso_2,v_positions[triangle[2]].x,v_positions[triangle[1]].x);
                    }

                    if(v_positions[triangle[1]].y < v_positions[triangle[2]].y){
                        yPos2 = lerp(iso_1,iso_2,v_positions[triangle[1]].y,v_positions[triangle[2]].y);
                    }else{
                        yPos2 = lerp(iso_1,iso_2,v_positions[triangle[2]].y,v_positions[triangle[1]].y);
                    }

                    segment_points.push_back(Point(xPos1,yPos1,0.0));
                    segment_points.push_back(Point(xPos2,yPos2,0.0));



                }else{ //case - + +

                }
            
            }
            else if(iso_1 >= 0 && iso_2 < 0){ // case - + - and - + +

            }


        }
    }



    // ------------------------------
}
