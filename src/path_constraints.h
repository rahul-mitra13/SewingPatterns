
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
                                        HalfedgeData<double>& sigmaCourseGlued);


//update glued halfedge weights 
//set the halfedges that have been take by some edge path i.e., gluedPath to infinity 
void updateGluedHalfedgeWeights(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, std::vector<double>& gluedPath,
                                HalfedgeData<double>& gluedHeWeights);