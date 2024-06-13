//C++ includes 
#include <optional>

//geometry-central includes 
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

//file includes
#include "helpers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this file contains helpers we use for the GUI, picking vertices, boundary elements etc
// 1. getAndRenderUserSpecifiedBoundaryVertices()
// 2. getAndRenderUserSpecifiedBoundaryInfo()


//get and render user specified boundary vertices. The path between these vertices defines the boundary conditions 
//
//@param[in]    geometry                VertexPositionGeometry      input geometry
//@param[in]    psMesh                  polyscope::SurfaceMesh      polyscope surface mesh
//
//@return       pair<Vertex, Vertex>    std::pair<Vertex, Vertex>   pair of vertices on the boundary whose path in between will specify the boundary vertices
std::pair<Vertex, Vertex> getAndRenderUserSpecifiedBoundaryVertices(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, int timeVal);

//get and render user specifed boundary paths (open boundaries) after performing a dijkstra search between user-specified boundary vertices 
//
//@param[in]    geometry                                    VertexPositionGeometry                                  input geometry
//@param[in]    psMesh                                      polyscope::SurfaceMesh                                  polyscope surface mesh 
//
//@return       pair<list of vertices, list of edges>       std::pair<std::vector<Vertex>, std::vector<Edge>>       pair(vertices on the boundary, edge on the boundary)
std::pair<std::vector<Vertex>, std::vector<Edge>> getAndRenderUserSpecifiedBoundaryInfo(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMeshes, int timeVal);