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
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap, 
                                                                    std::vector<std::pair<int, int>>& edgeMappingsPairs, polyscope::SurfaceMesh& psMesh, 
                                                                    globalBoundaryConditions& boundaryConditions, double period);

//compute the per-face gradient of a 1-form
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde);

//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//compute curl per edge in the global setting 
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                FaceData<Vector3>& globalFaceGradients, std::vector<std::pair<int, int>>& edgeMappingsPairs);

//find max/min curl vertex for a given isoline 
std::pair<int, int> findVertexSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllVertices);

//find max/min curl edge for a given isoline (on global setting)
std::pair<int, int> findEdgeSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllEdges);

//solve the optimization problem for strategy 1
std::tuple<HalfedgeData<double>, double> computeStrategy1_oneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, std::map<int, int>& globalToGluedVertexMap);

//find the isoval with max average deviation from \frac{\nabla h}{||h||}
double findIsoValWithMaxDeviation(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& gradSigmaTilde, FaceData<Vector3>& globalFaceGradients,
                                VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, std::map<int, int>& globalToGluedVertexMap);


//-------------------------Experiment 2-----------------------------//
//optimizing min \sum_{{ij} \in E}||\omegaTilde_{ij} - \sigmaTild_{ij}||^2 + \lambda \sum{e \in E}||\sigma_{ij} + \sigma_{ji}||^2
//visualizing where ||\sigma_{ij} + \sigma_{ji}||^2 is highest 
void vizEdgeDifference(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                        FaceData<Vector3>& globalFaceGradients, polyscope::SurfaceMesh& psMesh, 
                        globalBoundaryConditions& boundaryConditions, double period, double lambda,
                        std::map<int,int>& globalToGluedEdgeMap);