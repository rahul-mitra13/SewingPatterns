//C++ includes
#include <cfloat>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this file contains even more auxilliary functions used for experimentation, accessing basic mesh data etc
//1. getBoundaryVertices()
//2. getBoundaryEdges()
//3. getVertexPositionsandFaceLists()

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

