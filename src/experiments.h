//file includes
#include "knitting_utils.h"
#include "helpers.h"
#include "stripe_patterns_helpers.h"
#include "iterative_assignment.h"
#include <tuple>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//-----------------------------strategy1-------------------------------------------//
//compute the iterative 1-form and vertex singularity positions trying the strategy 1 
//strategy 1
//1. Sample level sets at intervals of P, say C_i 
//2. For each C_i, consider all the level sets it passes through, say T_i 
//3. Calc average ||\frac{\nabla h}{||\nabla h}||^2 - \grad \sigma_j} (probably area-weighted) over T_i
//4. For C_i with largest average ||\frac{\nabla h}{||\nabla h||, search max/min singularities and insert 
//5. Stop when ||\delta \sigma - \frac{\nabla h}{||\nabla h}||^2 doesn't decrease any more 
std::tuple<HalfedgeData<double>, VertexData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap, polyscope::SurfaceMesh& psMesh,
                                                                    std::map<int, int>& globalToGluedVertexMap, double period);

//compute the per-face gradient of a 1-form
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& sigmaTilde);

//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//find max/min curl vertex for a given isoline 
std::pair<int, int> findVertexSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap, 
                                            double isoVal, bool useAllVertices);