#pragma once

//C++ includes
#include <cfloat>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
#include <map>
#include <set>

//Eigen includes
#include <Eigen/Sparse>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"

//json include
#include "nlohmann/json.hpp"

//polyscope includes
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this file contains even more auxilliary functions used for experimentation, accessing basic mesh data etc
//1. getBoundaryVertices()
//2. getBoundaryEdges()
//3. getVertexPositionsandFaceLists()
//4. shortestEdgePathOnBoundary()

//global boundary conditions
struct globalBoundaryConditions{
    std::vector<int> courseStartBoundaryVertices;//vertices where t = 0
    std::vector<int> courseEndBoundaryVertices;//vertices where t = 1

    std::vector<int> courseBdyEdges;//boundary edges in the course direction where sigma = 0
    std::vector<std::vector<double>> waleBdyPathConstraints;//vector of vectors (each of size nEdges()) where each entry stores the weight of the edge in the integration with \sigma
};


//
//get a pair of vertex lists that are at the extremes of a panel (this probably needs user input eventually)
//
//@param[in]    geometry                    VertexPositionGeometry                  input geometry
//@param[in]    axis                        int                                     axis whose extremes we're interested in (0 -> x, 1 -> y, 2 -> z)
//
//@return       pair of two vertex lists    pair<vector<Vertex>, vector<Vertex>>    pair of vertex lists that represent the vertices at the extremes
std::pair<std::vector<Vertex>, std::vector<Vertex>> getBoundaryVertices(VertexPositionGeometry& geometry, int axis);

//
//get a pair of boundary edges that are at the extremes of a panel (this probably needs user input eventually)
//@param[in]    geometry                VertexPositionGeometry              input geometry 
//@param[in]    axis                    int                                 axis whose extremes we're interested in (0 -> x, 1 -> y, 2 -> z)
//
//@return       pair of 2 edge lists    pair<vector<Edge, vector<Edge>>     pair of edge lists that represent the edges at the extremes 
std::pair<std::vector<Edge>, std::vector<Edge>> getBoundaryEdges(VertexPositionGeometry& geometry, int axis);

//
// Project a vector onto a given plane.
//
// @param[in]  vec     The vector to project onto the plane.
// @param[in]  normal  Normal to the plane to project onto.
// @param[in]  axis    Axis vector whose projection onto the plane defines the local X axis.
//
// @return     Projected vector.
//
Vector2 projectOntoPlane(const Eigen::Vector3d &vec, const Eigen::Vector3d &normal, const Eigen::Vector3d &axis);

// Get vertex positions and face lists from an input geometry.
//
// @param[in]   geometry                VertexPositionGeometry                          input geometry  
//
//@return       pair of 2 matrices      pair(Eigen::MatrixXd V, Eigen::MatrixXi F)      V -> vertex positions, F -> face index list
//
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> getVertexPositionsandFaceLists(VertexPositionGeometry& geometry);

//get shortest edge path between two vertices by taking only boundary edges of the mesh
//
//@param[in]    geom            VertexPositionGeometry      input geometry
//@param[in]    startVert       Vertex                      start vertex on the boundary 
//@param[in]    endVert         Vertex                      end vertex on the boundary
//
//@return       halfedge list   std::vector<Haledge>        returns a list of halfedges on the shortest boundary path between the two vertices    
std::vector<Halfedge> shortestEdgePathOnBoundary(VertexPositionGeometry& geom, Vertex startVert, Vertex endVert);

//get the shortest dijkstra edge path on the boundary
//
//@parma[in]     geom               VertexPositionGeometry                                          input geometry 
//@param[in]     startVert          Vertex                                                          start vertex on the boundary 
//@param[in]     endVert            Vertex                                                          end vertex on the boundary 
//
//@return        tuple              std::tuple<vector<Vertex>, vector<Edge>, vector<double>>        returns a list of vertices and edges on the boundary between the two vertices and edge weights on this he path
std::tuple<std::vector<Vertex>, std::vector<Edge>, std::vector<double>> getVerticesAndEdgesInShortestEdgePathOnBoundary(VertexPositionGeometry& geom, Vertex startVert, Vertex endVert);


//build a vertex mapping map from an input txt file (for a "local" mapping across different models)
//
//@param[in]    filename    const std::string&                                                  name of the file that stores the mappings 
//@return       map         std::map<std::pair<std::string, int>, std::pair<std::string, int>>  this is a map (panel1name, panel1Vertex -> panel2name, panel2Vertex)
std::map<std::pair<std::string, int>, std::pair<std::string, int>> buildLocalVertexMappingMapFromFile(const std::string& filename);

//builds a vector of "stitched together" vertices
std::vector<std::pair<int, int>> buildPairOfStitchedVerticesFromFile(const std::string& filename);

//builds a vector of "stitched together" edges from a vector of "stitched" together vertices 
std::vector<std::pair<int, int>> buildPairOfStitchedEdges(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs);

//takes in a vector of vertex indices and returns a list of edge indices making up that vertex list
std::vector<int> creatEdgeListFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList);

//takes in a vector of vertex indices in order and returns edge weights corresponoding to those vertices
std::vector<double> createEdgeWeightsFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList);

//renders "stitched" together vertices from the 2D panels 
//
//@param[in]    geometry        VertexPositionGeometry  input geometry
//@param[in]    vertexMappings  std::map<int, int>      map that stores global vertex mappings 
//
//@return       none                                    but has a side effect where it calls polyscope to render stitches
void renderStitchedVertices(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs);

//parse json file and render boundary conditions
globalBoundaryConditions parseJson(VertexPositionGeometry& geometry, nlohmann::json& data);

//create a new "glued surface mesh"
//create a mapping from original index to "glued" index (vertex and edges)
//also populate edge lengths of glued mesh
SurfaceMesh * createGluedSurfaceMesh(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs, std::map<int,int>& originalMeshVertexIndexToGluedMeshIndex, 
                                    std::map<int, std::vector<Halfedge>>& gluedOneRingMap);