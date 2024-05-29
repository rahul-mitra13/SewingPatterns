//C++ includes
#include <cfloat>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this file contains even more auxilliary functions used for experimentation, accessing basic mesh data etc
//1. getBoundaryVertices()

//
//get a pair of vertices that are at the extremes of a panel
//
//@param[in]    geometry                    VertexPositionGeometry                  input geometry
//@param[in]    axis                        int                                     axis whose extremes we're interested in (0 -> x, 1 -> y, 2 -> z)
//
//@return       pair of two vertex lists    pair<vector<Vertex>, vector<Vertex>>    pair of vertex lists that represent the vertices at the extremes
std::pair<std::vector<Vertex>, std::vector<Vertex>> getBoundaryVertices(VertexPositionGeometry& geometry, int axis);
