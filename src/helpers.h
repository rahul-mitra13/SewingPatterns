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
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/utilities/disjoint_sets.h"

//json include
#include "nlohmann/json.hpp"

//polyscope includes
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

// Include Gurobi
#include "gurobi_c++.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Some useful macros :-)
#define P(x) {if (debug) std::cout << x << std::endl;} 
#define H(x) P(#x << ": " << (x))
const bool debug = 1;

// To print pairs easily
template<class T1, class T2> std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> v) {
  os << "(" << v.first << ", " << v.second << ")";
  return os;
}

// To print vectors easily
template<class T> std::ostream &operator<<(std::ostream &os, std::vector<T> v) {
  os << "["; if (v.size() > 0) os << v[0];
  for(int i = 1; i < v.size(); i++) os << ", " << v[i];
  os << "]";
  return os;
}

// Our own assert() equivalent
#undef ensure
#define ensure(x)                                         \
  if (!(x)) {                                             \
    std::cout << "Assertion failed: " << #x << std::endl; \
    exit(1);                                               \
  }

// Global configuration options
struct Options {

  // Add Polyscope quantities for every iteration of the algorithm
  bool showAllIterations = true;

};

//boundary conditions on the GLUED MESH
struct globalBoundaryConditions{
    std::vector<int> courseStartBoundaryVertices;//vertices where t = 0
    std::vector<int> courseEndBoundaryVertices;//vertices where t = 1

    std::vector<int> courseBdyEdges;//boundary edges in the course direction where sigma = 0
    std::vector<std::vector<double>> waleBdyPathConstraints;//vector of vectors (each of size nEdges()) where each entry stores the weight of the edge in the integration with \sigma (path around a boundary)

    std::vector<std::vector<double>> bdyBdyPathConstraints;//vector of vectors (each of size nEdges()) where each entry stores the weight of the edge in the integration with \sigma (boundary-boundary path)
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
// @param[in]   geometry                IntrinsicGeometryInterface                      input geometry  
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
std::vector<Halfedge> shortestEdgePathOnBoundary(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert);


//get the shortest dijkstra edge path on the boundary
//
//@parma[in]     geom               VertexPositionGeometry                                          input geometry 
//@param[in]     startVert          Vertex                                                          start vertex on the boundary 
//@param[in]     endVert            Vertex                                                          end vertex on the boundary 
//
//@return        tuple              std::tuple<vector<Vertex>, vector<Edge>, vector<double>>        returns a tuple of vertices and edges on the boundary between the two vertices and edge weights on this he path
std::tuple<std::vector<Vertex>, std::vector<Edge>, std::vector<double>> getVerticesAndEdgesInShortestEdgePathOnBoundary(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert);

//same as above but uses the shortest edge path function instead of shortest edge path on boundary function 
std::tuple<std::vector<Vertex>, std::vector<Edge>, std::vector<double>> getVerticesAndEdgesInShortestEdgePath(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert);


//build a vertex mapping map from an input txt file (for a "local" mapping across different models)
//
//@param[in]    filename    const std::string&                                                  name of the file that stores the mappings 
//@return       map         std::map<std::pair<std::string, int>, std::pair<std::string, int>>  this is a map (panel1name, panel1Vertex -> panel2name, panel2Vertex)
std::map<std::pair<std::string, int>, std::pair<std::string, int>> buildLocalVertexMappingMapFromFile(const std::string& filename);

//builds a vector of "stitched together" vertices
//
//@param[in]    filename    const std::string&                  this is the name of the file that stores the vertex mappings in the global sense 
//
//@return       map         std::vector<std::pair<int, int>>    returns a vector of pairs where each pair is a vertex mapping 
std::vector<std::pair<int, int>> buildPairOfStitchedVerticesFromFile(const std::string& filename);

//builds a vector of "stitched together" edges from a vector of "stitched" together vertices 
//
//@param[in]    geometry                VertexPositionGeometry                  input geometry 
//@param[in]    vertexMappingsPairs     std::vector<std::pair<int, int>>        a vector of vertex mappings pairs
//
//@return       vector                  std::vector<std::pair<int, int>>        vector of pairs where each pair denotes edges stitched together
std::vector<std::pair<int, int>> buildPairOfStitchedEdges(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs);

//takes in a vector of vertex indices and returns a list of edge indices making up that vertex list
//
//@param[in]    geometry            VertexPositionGeometry      input geometry 
//@param[in]    vertexList          std::vector<int>            index of vertices we want to create the edge list from 
//
//@return       edgeList            std::vector<int>            vector of edge indices that are the edges in the vertex list 
std::vector<int> creatEdgeListFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList);

//takes in a vector of vertex indices in order and returns edge weights corresponoding to those vertices
//
//@param[in]    geometry        VertexPositionGeometry      input geometry 
//@param[in]    vertexList      std::vector<int>            list of vertices in the edge path
//
//@return       weights         std::vector<double>         a vector of length nEdges that stores the edge weights of the path that the ordered vertices specify
std::vector<double> createEdgeWeightsFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList);

//renders "stitched" together vertices from the 2D panels 
//
//@param[in]    geometry        VertexPositionGeometry  input geometry
//@param[in]    vertexMappings  std::map<int, int>      map that stores global vertex mappings 
//
//@return       none                                    but has a side effect where it calls polyscope to render stitches
void renderStitchedVertices(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs);

//get the indices of the boundary edges in the glued mesh where we should set \sigma_{wale} = 0
//
//@param[in]    globalGeometry      VertexPositionGeometry          global input geometry 
//@param[in]    gluedGeometry       IntrinsicGeometryInterface      glued geometry 
//@param[in]    faceGradients       FaceData<Vector3>               course face gradients 
//@param[in]    edgeMap             std::map<int, int>              map from edges in the global mesh to edges in the glued mesh 
//@param[in]    val                 double                          threshold value
//@param[in]    globalPSMesh        polyscope::SurfaceMesh          polyscope mesh for visualizing/debugging
//
//@return       constrainedEdges    std::vector<int>                indices of edges (in the glued mesh) where sigma_{wale} = 0
std::vector<int> getWaleBdyEdgesInGluedMesh(VertexPositionGeometry& globalGeometry, IntrinsicGeometryInterface& gluedGeometry, FaceData<Vector3>& faceGradients, std::map<int, int>& edgeMap, double val,
                                            polyscope::SurfaceMesh& globalPSMesh);

//parse the global boundary conditions in the glued together mesh 
globalBoundaryConditions parseJson(IntrinsicGeometryInterface& gluedGeometry, nlohmann::json& data, std::map<int,int>& vertexMap, std::map<int, int>& edgeMap);

//create a new "glued" edge length geometry
//create a mapping from original index to "glued" index (vertex and edges)
//also populate edge lengths of glued mesh
// EdgeLengthGeometry * createGluedEdgeLengthGeometry(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs, std::map<int,int>& vertexMap, 
//                                     std::map<int, int>& edgeMap, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

std::unique_ptr<EdgeLengthGeometry> createGluedEdgeLengthGeometry(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs, std::map<int,int>& vertexMap, 
                                    std::map<int, int>& edgeMap, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);


//some utility functions for converting values defined on the glued mesh to values defined on the global mesh 
//used for viz and other routines possibly?
//probably going to be a little slow but semantically I think this is the best thing to do

//convert a function defined on the vertices of the glued mesh to a function defined on the vertices of the global mesh 
VertexData<double> convertGluedToGlobalVertexFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& func, std::map<int, int>& vertexMap);

//convert a line field defined on the global mesh to a line field defined on the glued mesh 
VertexData<Vector2> convertGlobalToGluedVertexLineField(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<Vector2>& globalLineField, std::map<int, int>& vertexMap);

//convert function defined on the edges of the glued mesh to a function defined on the edges of the global mesh 
EdgeData<double> convertGluedToGlobalEdgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& func, std::map<int, int>& edgeMap);

//convert function defined on the edges of the global mesh to a function defined on the edges of the glued mesh 
EdgeData<double> convertGlobalToGluedEdgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& func, std::map<int, int>& edgeMap);

//convert function defined on the halfedges of the glued mesh to a function defined on the halfedges of the global mesh 
HalfedgeData<double> convertGluedToGlobalHalfedgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& func, std::map<int, int>& edgeMap);

//draw the mesh curve network 
void drawMeshCurveNetwork(VertexPositionGeometry& globalGeometry, polyscope::SurfaceMesh& psMesh);

//@clean
//after finding a face singularity, find a singular edge in that face 
//edge that is most aligned with the gradient (or rotated gradient in the wale case) of the time function 
//@param[in]        globalGeometry      VertexPositionGeometry      the global geometry 
//@param[in]        gluedGeometry       EdgeLengthGeometry          glued edge length geometry 
//@param[in]        singFaceIndex       int                         index of the singular face
//@param[in]        globalFaceGradient  Vector3                     the gradient vector at the singular face 
//
//@return           maxEdge             int                         index of the singular edge in the global mesh setting
int findSingularEdgeFromSingularFace(VertexPositionGeometry& globalGeometry, int singFaceIndex, Vector3 globalFaceGradient, double threshold, 
                                    VertexData<double>& globalTimeFunction,
                                    double isoVal);


//@clean 
//hash a floating point number
//still causes collision when numbers are extremely close (hopefully we won't run into this too much)
//@param[in]    f           double      floating point number to hash 
//
//@return       hashed[f]    int         hashed floating point number
int hashFloatQuantized(double f);

//@clean 
//find the number of connected components in a set of faces 
//a single isoline of the time function passes through all these faces
std::vector<std::vector<int>> findConnectedComponents(EdgeLengthGeometry& gluedGeometry, std::vector<int>& faces);

// For some reason I cannot find the right permutation to make addCornerScalarQuantity work correctly.
// I think there is something wrong with CornerData (e.g., note how cornerData.size() != mesh.nCorners()),
// but I don't have the time to look into it.
// In the meantime, here's a wrapper that correctly prepares corner data;
// just call it and feed the output to addCornerScalarQuantity.
std::vector<double> prepareCornerData(CornerData<double> cornerData);

std::vector<double> prepareHalfedgeData(HalfedgeData<double> halfedgeData);

// A Union Find data structure
// We use it for the vertex mappings in the glued mesh.
class UF {
  int *id, cnt, *sz;
  public:
  // Create an empty union find data structure with N isolated sets.
  UF(int N) {
      cnt = N; id = new int[N]; sz = new int[N];
      for (int i = 0; i<N; i++)  id[i] = i, sz[i] = 1; }
  ~UF() { delete[] id; delete[] sz; }

  // Return the id of component corresponding to object p.
  int find(int p) {
      int root = p;
      while (root != id[root])    root = id[root];
      while (p != root) { int newp = id[p]; id[p] = root; p = newp; }
      return root;
  }
  // Replace sets containing x and y with their union.
  void merge(int x, int y) {
      int i = find(x); int j = find(y); if (i == j) return;
      // make smaller root point to larger one
      if (sz[i] < sz[j]) { id[i] = j, sz[j] += sz[i]; }
      else { id[j] = i, sz[i] += sz[j]; }
      cnt--;
  }
  // Are objects x and y in the same set?
  bool connected(int x, int y) { return find(x) == find(y); }
  // Return the number of disjoint sets.
  int count() { return cnt; }
};

// Given a set of continuous numbers, round them to the closest integers while maintaining a given sum
std::vector<int> roundWithSum(std::vector<double> y, int sum);

// Show a set of edges as a curve network
polyscope::CurveNetwork* showEdges(std::string name, std::vector<Edge> edges, VertexPositionGeometry& geometry);