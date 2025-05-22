
//eigen includes
#include <Eigen/Dense>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/surface/stripe_patterns.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"


//C++ includes
#include <cfloat>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

//include gurobi
#include "gurobi_c++.h"



using namespace geometrycentral;
using namespace geometrycentral::surface;

//find maximum dot product of a halfedge with the clockwise rotated face gradients
double maximumDotProduct(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients);

//rotate the face gradients clockwise 
FaceData<Vector3> clockWiseRotatedGradients(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients);

//construct glued halfedge weights 
HalfedgeData<double> constructGluedHalfedgeWeights(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& rotatedFaceGradients,
                                            double maxDotProd);

//construct an edge path between two singular edges 
//e1, e2 are in the global setting 
std::tuple<std::vector<double>, std::vector<double>> constructEdgePath(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Edge e1, Edge e2, 
                                        std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, FaceData<Vector3>& globalFaceGradients, HalfedgeData<double>& gluedHeWeights,
                                        HalfedgeData<double>& sigmaCourseGlued, bool connectSaddles=true);


//update glued halfedge weights 
//set the halfedges that have been take by some edge path i.e., gluedPath to infinity 
void updateGluedHalfedgeWeights(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, std::vector<double>& gluedPath,
                                HalfedgeData<double>& gluedHeWeights);


//given a set of singularities perform an optimal matching between them 
std::vector<std::pair<Vertex, Vertex>> performOptimalMatching(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& heWeights, std::vector<std::pair<Vertex, int>>& singularities);

//compute the path cost between 2 vertices
double computePathCost(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& heWeights, Vertex& startVert, Vertex& endVert);


//--------------New path tracing mechanism-------------------//
//find an isopath between two vertices
std::vector<double> findIsoPath(SurfaceMesh& mesh, VertexData<double>& timeFunction, Vertex vStart, Vertex vEnd);

//find triangle strip containing an isoline
// FaceData<int> traceIsolineToEdge(SurfaceMesh& mesh,
//                                              const VertexData<double>& timeFunction,
//                                              Edge startEdge,
//                                              Edge endEdge,
//                                              double isoVal);


// Traces the isoline at a given isovalue starting from startEdge to endEdge
std::vector<Edge> traceIsoline(
    const SurfaceMesh& mesh,
    const VertexData<double>& timeFunction,
    double isoVal,
    Edge startEdge,
    Edge endEdge
);

// Traces the isoline at a given isovalue starting from startHe to endHe
//same as above just uses halfedges instead of edges
std::vector<Halfedge> traceIsoline(
    const SurfaceMesh& mesh,
    const VertexData<double>& timeFunction,
    double isoVal,
    Halfedge startHe,
    Halfedge endHe
);

//trace the faces that an isoline passes through
std::vector<Face> traceIsolineFaces(
    const VertexPositionGeometry& globalGeometry,
    const VertexData<double>& timeFunction,
    const FaceData<Vector3>& rotatedFaceGradients,
    double isoVal,
    Halfedge startHe,
    Halfedge endHe
);