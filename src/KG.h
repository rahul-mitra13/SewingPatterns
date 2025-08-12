//geometry-central includes 
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/surface_point.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

//file includes 
#include "helpers.h"

#pragma once 
#include <vector> 
#include <optional>
#include <iostream>
#include <fstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//redefining for consistency
struct KGVertex{
    int id = -1;
    Vector3 position;//position embedded in R^3 (if position information is available)
    Vector3 baryCoords;//barycentric coordinates of the vertex if position information is not available
    int row_in = -1;
    int row_out = -1;
    int col_in[2] = {-1, -1};
    int col_out[2] = {-1, -1};
    double alpha_tag = -1;
    double beta_tag = -1;
    Face face;//associated face of a vertex (on the underlying mesh)
    std::optional<Edge> edge;//associated edge of a vertex
    std::optional<Halfedge> halfedge;//associated halfedge of a vertex
    SurfacePoint surfacePoint;//used for intersections on singular faces
    bool isBaryCenter = false;//if this is a vertex at the barycenter
    bool hasBeenHandled = false;//if this vertex has been handled by a merge
    bool isVirtual = false;//if this is a virtual vertex
    bool isAlphaVirtual = false;//is a virtual vertex in the course direction
    bool isBetaVirtual = false;//is a virtual vertex in the wale direction 

    void printVertexInfo(){
          
        std::cout << "vertex id = " << id << std::endl;
        std::cout << "row in = " << row_in << std::endl;
        std::cout << "row out = " << row_out << std::endl;
        std::cout << "col_in[0] = " << col_in[0] << std::endl;
        std::cout << "col_in[1] = " << col_in[1] << std::endl;
        std::cout << "col_out[0] = " << col_out[0] << std::endl;
        std::cout << "col_out[1] = " << col_out[1] << std::endl;
        std::cout << "isVirtual = " << isVirtual << std::endl;
    }
};