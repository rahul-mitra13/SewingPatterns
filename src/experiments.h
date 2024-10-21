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
//2. For each C_i, consider all the edges it passes through, say T_i 
//3. Calc average edge curl over every T_i
//4. For C_i with largest average edge curl, search max/min edge curl and insert singularities
//5. Stop when ||\delta \sigma - \frac{\nabla h}{||\nabla h}||^2 doesn't decrease any more 
std::tuple<HalfedgeData<double>, EdgeData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap, 
                                                                    polyscope::SurfaceMesh& psMesh, globalBoundaryConditions& boundaryConditions, double period);

//trying to solve for a 1-form that only optimizes for equally spaced stripes 
//Matteo's idea: The objective term is || ||\delta \sigma||^2 - 1||^2 
//Have a bdy-bdy path integral for non-collapse and directionality specification
std::tuple<HalfedgeData<double>, double> computeEquallySpacedOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                                std::map<int, int>& globalToGluedVertexMap);

//compute the per-face gradient of a 1-form
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde);

//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//compute curl per edge in the global setting 
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                FaceData<Vector3>& globalFaceGradients, std::map<int, int>& globalToGluedEdgeMap);

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
double findIsoValWithMaxVertexCurlDeviation(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& gradSigmaTilde, FaceData<Vector3>& globalFaceGradients,
                                VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, std::map<int, int>& globalToGluedVertexMap);

//find the isoval with max averga curl in the edge setting
double findIsoValWithMaxAvgEdgeCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& curl, 
                                            VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, 
                                            std::map<int, int>& globalToGluedVertexMap);


//-------------------------Experiment 2-----------------------------//
//optimizing min \sum_{{ij} \in E}||\omegaTilde_{ij} - \sigmaTild_{ij}||^2 + \lambda \sum{e \in E}||\sigma_{ij} + \sigma_{ji}||^2
//visualizing where ||\sigma_{ij} + \sigma_{ji}||^2 is highest 
//for various values of \lambda
void vizEdgeDifference(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                        FaceData<Vector3>& globalFaceGradients, polyscope::SurfaceMesh& psMesh, 
                        globalBoundaryConditions& boundaryConditions, double period, double lambda,
                        std::map<int,int>& globalToGluedEdgeMap);


//------------------------Experiment 3----------------------------//
//Trying to find a harmonic 1-form in the halfedge optimization setting 
//use the thus far inserted singular edges (and edge-based curl of zero elsewhere), 
//enforce face-based curl of zero, and an integral of 1 along some path from either boundary to optimize 
//for a harmonic (half-edge) one-form subject to the singularity placements; do this at each 
//iteration and then normalize to find edge-based curl
std::tuple<HalfedgeData<double>, EdgeData<double>> harmonic1FormImpl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap, std::vector<std::pair<int, int>>& edgeMappingsPairs,
                                                                    polyscope::SurfaceMesh& psMesh, globalBoundaryConditions& boundaryConditions, double period);

//solve the optimization problem for the harmonic 1-form
std::tuple<HalfedgeData<double>, double> computeHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, std::map<int, int>& globalToGluedVertexMap, 
                                                            std::vector<std::pair<int, int>>& edgeMappingsPairs, polyscope::SurfaceMesh& psMesh);

//compute face curl by averaging edge curl over the edges in a face 
FaceData<double> computeAverageEdgeCurlonFaces(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                FaceData<Vector3>& globalFaceGradients, std::map<int, int>& globalToGluedEdgeMap);

//find max/min curl face for a given isoline of the TIME FUNCTION 
//would probably want to change tracing the level sets of the time function with level sets of the harmonic 1-form
std::pair<int, int> findFaceSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllFaces);

//after finding a face singularity, find a singular edge in that face 
//edge that is most aligned with the gradient of the time function
int findSingularEdgeFromSingularFace(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, int singFaceIndex, Vector3 globalFaceGradient);

//find the isoval with max average curl in the face setting
double findIsoValWithMaxAvgFaceCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<double>& curl, 
                                            VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, 
                                            std::map<int, int>& globalToGluedVertexMap);