#pragma once

//C++ includes
#include <cfloat>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
#include <map>

//Eigen includes
#include <Eigen/Sparse>


//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this file contains even more auxilliary functions used for experimentation, accessing basic mesh data etc
//1. getBoundaryVertices()
//2. getBoundaryEdges()
//3. getVertexPositionsandFaceLists()
//4. shortestEdgePathOnBoundary()

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
//@parma[in]     geom               VertexPositionGeometry                          input geometry 
//@param[in]     startVert          Vertex                                          start vertex on the boundary 
//@param[in]     endVert            Vertex                                          end vertex on the boundary 
//
//@return        pair               std::pair<vector<Vertex>, vector<Edge>>         returns a list of vertices and edges on the boundary between the two vertices
std::pair<std::vector<Vertex>, std::vector<Edge>> getVerticesAndEdgesInShortestEdgePathOnBoundary(VertexPositionGeometry& geom, Vertex startVert, Vertex endVert);

//build the global cotan Laplacian for all the panels while accounting for the mapped stitches across panels
void buildGlobalCotanLaplacian(std::map<std::string, std::unique_ptr<VertexPositionGeometry>>& panelMappings, std::map<std::pair<std::string, int>, 
                                std::pair<std::string, int>>& vertexMappings, Eigen::SparseMatrix<double>& L);


//build a vertex mapping map from an input txt file 
std::map<std::pair<std::string, int>, std::pair<std::string, int>> buildVertexMappingMapFromFile(const std::string& filename);

//structure that stores boundary conditions for every patch 
struct PatchBoundaryConditions{
    //boundary vertices in the course direction
    std::vector<Vertex> courseStartBoundaryVertices;
    std::vector<Vertex> courseEndBoundaryVertices;
    //boundary vertices in the wale direction
    std::vector<Vertex> waleStartBoundaryVertices;
    std::vector<Vertex> waleEndBoundaryVertices;

    //boundary edges in the course direction
    std::vector<Edge> courseStartBoundaryEdges;
    std::vector<Edge> courseEndBoundaryEdges;
    //boundary edges in the wale direction
    std::vector<Edge> waleStartBoundaryEdges;
    std::vector<Edge> waleEndBoundaryEdges;
};

