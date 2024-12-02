#ifndef KNITTING_UTILS_H
#define KNITTING_UTILS_H

//eigen includes
#include <Eigen/Dense>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/surface/stripe_patterns.h"
#include "geometrycentral/surface/barycentric_vector.h"

//libigl includes
//I guess we'll just put all the libigl includes here for now 
#include <igl/cotmatrix.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/grad.h>
#include <igl/average_onto_vertices.h>
#include <igl/hessian_energy.h>
#include <igl/curved_hessian_energy.h>
#include <igl/isolines.h>

//helper files
#include "helpers.h"
#include "stripe_patterns_helpers.h"
#include "path_constraints.h"

//includes to solve the optimization problem
#include "gurobi_c++.h"
#include "model.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//compute time function using a vector of pairs of vertex mappings instead of the map because we miss stitches then 
//
//@param[in]    geometry        VertexPositionGeometry                  input geometry
//@param[in]    vertexMappings  std::vector<std::pair<int,int>>         map that stores global vertex mappings 
//
//@return       timeFunction    VertexData<double>                      the global time function computed over the entire mesh
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::vector<std::pair<int,int>>& vertexMappingsPairs,  globalBoundaryConditions& boundaryConditions, std::map<int, int>& indexMap);

//compute time function directly in the glued mesh setting 
VertexData<double> computeTimeFunction(EdgeLengthGeometry& gluedGeometry, globalBoundaryConditions& bdyConditions);


//compute the per face gradient of some function defined as a scalar over vertices
//
//@param[in]    geometry                VertexPositionGeometry          input geometry
//@param[in]    vertexScalarFunction    VertexData<double>              scalar function defined as a scalar over vertices
//
//@return       faceGradients           FaceData<Vector3>               #F by 3 vector per face
FaceData<Vector3> computeTimeFunctionFaceGrad(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction);

//compute the gradient of the a function defined as a scalar over vertices in the glued mesh setting 
//Don't think the formula I'm using in here is right
FaceData<Vector3> computeTimeFunctionFaceGrad(EdgeLengthGeometry& geometry, VertexData<double>& vertexScalarFunction);

//compute a vector (in ambient space) per vertex that is aligned with the gradient of a scalar field
//
//@param[in]    geometry                VertexPositionGeometry  input geometry
//@param[in]    vertexScalarFunction    VertexData<double>      scalar function defined as a scalar over vertices
//@param[in]    angle                   double                  how much to rotate the vector about the vertex normal axis in radians
//
//@return       vertexValuedField       VertexData<Vector3>     #V by 3 vector per face     
VertexData<Vector3> computeVertexValuedField(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction, double angle);


//compute a direction field per vertex from a vector field (in ambient space) per vertex
//
//@param[in]    geometry            VertexPositionGeometry geometry     input geometry
//@param[in]    vertexValuedField   VertexData<Vector3>                 #V by 3 vector per vertex
//
//@return       directionField      VertexData<Vector2>                 direction field at each vertex
VertexData<Vector2> vertexDirectionField(VertexPositionGeometry& geometry, VertexData<Vector3>& vertexValuedField);

//compute a 1-form that will be used to generate the stripes over the patches
//
//@param[in]    geometry            VertexPositionGeometry              input geometry 
//@param[in]    model               Model                               The gurobi model we're interested in solving
//@param[in]    edgeMappingPairs    std::vector<std::pair<int, int>>    vector of pairs where each pair representes "stitched together" edge indices   
//
//@return       oneForm     EdgeData<double>        The course 1-form used to specify the stripes
EdgeData<double> computeOneForm(VertexPositionGeometry& geometry, Model& gbModel, polyscope::SurfaceMesh& globalPSMesh);

//compute a 1-form that will be used to generate the stripes over the pathces 
//overload from the function above
//it does the optimization in the glued mesh setting 
EdgeData<double> computeOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& vertexMap, std::map<int, int>& edgeMap,
                                        polyscope::SurfaceMesh& psMesh);

//compute a face-based field through the optimization 
//here the singularities are placed on the vertices
HalfedgeData<double> computeVertexSingularityField(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, polyscope::SurfaceMesh& psMesh, std::map<int, int>& vertexMap,
                                                  std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//compute \omega i.e., the 1-form we're trying to match over each edge 
//
//@param[in]    geometry            VertexPositionGeometry  input geometry
//@param[in]    direction           int                     integer representing the direction for omega that we're attempting to generate (0 -> course, 1 -> wale)
//@param[in]    faceGradients       FaceData<Vector3>       the gradients of the time function we're using to solve for the matching 1-form
//
//@return       matchingOneForm     EdgeData<double>        The 1-form we're trying to match 
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3>& faceGradients);

//compute matching 1-form while taking into account "stitched together" edges
//same parameters and return type as the function above 
//additional parameter 
//@param[in]    edgeMappingsPairs       std::vector<std::pair<int, int>>   vector of pairs where each pair stores edges that are stitched together
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry,  int direction, FaceData<Vector3>& faceGradients, std::vector<std::pair<int, int>>& edgeMappingsPairs);

//implement the harmonic 1-form 
//the below two function use the energy min cot_e ||\sigma_e||^2
// std::tuple<CornerData<double>, EdgeData<double>> implCourseHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
//                                                                     VertexData<double>& globalTimeFunction, FaceData<Vector3>& globalTimeFunctionGradientsNormalized,
//                                                                     std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, polyscope::SurfaceMesh& psMesh,
//                                                                     globalBoundaryConditions& boundaryConditions, double period,
//                                                                     Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double, Eigen::RowMajor>& G);

//compute course harmonic 1-form
// std::tuple<HalfedgeData<double>, double> computeHarmonicCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
//                                                         std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh);

//-------------------------------------------------------------------------------------------//

//@clean: tag represents code that should be written cleanly and change minimally going forward :D 

//@clean
//compute the per-face gradient of a 1-form (in global setting)
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde);

//@clean 
//same function as above but just returns a different object
std::vector<std::array<double, 3>> vectorOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde);

//@clean 
//compute a 1-form that will be used to generate the wale stripes 
//@param[in]    globalGeometry      VertexPositionGeometry  geometry of the object in the global setting
//@param[in]    gluedGeometry       EdgeLengthGeometry      edge length geometry in the glued setting 
//@param[in]    model               Model                   model we will be solving using gurobi (specifies the wale constraints)
//@param[in]    vertexMap           std::map<int, int>      a map from global vertex indices to glued vertex indices
//
//@return       sigma               HalfedgeData<double>    1-form value per halfedge representing the wale stripes
HalfedgeData<double> computeWaleOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, 
                                        Eigen::SparseMatrix<double, Eigen::RowMajor>& G, std::map<int, int>& vertexMap);

//@clean 
//return information for wale stripes in the glued mesh setting
//@param[in]    globalGeometry                          VertexPositionGeometry                      global geometry embedded in R^3
//@param[in]    gluedGeometry                           EdgeLengthGeometry                          edge length geometry in the glued setting
//@param[in]    edgeMappingsPaits                       std::vector<std::pair<int, int>>            vector of pairs where each pair stores edges that are stitched together
//@param[in]    edgeMap                                 std::map<int,int>                           a map from edges in the global mesh to edges in the glued mesh
//@param[in]    vertexMap                               std::map<int, int>                          a map from global vertex indices to glued vertex indices
//@param[in]    timeFunctionGlobal                      VertexData<int>                             time function in the global setting
//@param[in]    courseOneFormGrad                       FaceData<Vector3>                           normalized gradient of the final course one form
//@param[in]    period                                  double                                      period for 1-form optimization
//@param[in]    knoppelFrequency                        double                                      period to be used when generating Knoppel stripes
//@param[in]    globalBdyConditions                     globalBoundaryConditions                    boundary conditions specified in the glued mesh setting
//
//@return       striping info       tuple<CornerData<double, EdgeData<double>>     striping info in the glued mesh setting 
std::tuple<CornerData<double>, EdgeData<double>> computeWaleStripeInfo(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    std::vector<std::pair<int, int>>& edgeMappingsPairs, std::map<int, int>& edgeMap, 
                                                                    std::map<int, int>& vertexMap, VertexData<double>& timeFunctionGlobal, FaceData<Vector3>& courseOneFormGrad, 
                                                                    Eigen::SparseMatrix<double, Eigen::RowMajor>& G, double period, double knoppelFrequency, globalBoundaryConditions& globalBdyConditions,
                                                                    EdgeData<double>& courseSingularEdgesGlobal, polyscope::SurfaceMesh& psMesh);


//@clean 
//compute the edge curl in the global setting
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry,
                                FaceData<Vector3>& globalFaceGradients);

//compute edge curl in the global setting 
//compute curl per edge in the global setting 
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients);

//@clean
//compute face curl by averaging edge curl over the edges in a face 
FaceData<double> computeAverageEdgeCurlonFaces(VertexPositionGeometry& globalGeometry, EdgeData<double>& edgeCurl);


//@clean
//find max/min curl face for a given isoline of the TIME FUNCTION 
//would probably want to change tracing the level sets of the time function with level sets of the harmonic 1-form
std::vector<std::pair<int, int>> findFaceSingularityPairs(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, FaceData<double>& curl,
                                            VertexData<double>& globalTimeFunction, polyscope::SurfaceMesh& psMesh,
                                            std::map<int, int>& hashedUsedIsoVals, FaceData<double>& faceSingularities, FaceData<int>& forbiddenFaces,
                                            double isoVal, int numPairs, bool useAllFaces);


//@debugging
//find max/min curl edge for a given isoline of the time function 
//find edge singularity pair 
std::vector<std::pair<int, int>> findEdgeSingularityPairs(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, EdgeData<double>& curl,
                                            VertexData<double>& globalTimeFunction, polyscope::SurfaceMesh& psMesh,
                                            std::map<int, int>& hashedUsedIsoVals, std::map<int, int>& usedEdges, FaceData<double>& faceSingularities, 
                                            double isoVal, int numPairs, bool useAllEdges);


//@debugging
//find edge singularity pairs 
//this method finds the max curl edges and finds an edge on the same isoline of the time function 
//of the opposite sign
std::vector<std::pair<int, int>> findEdgeSingularityPairsUsingMaxCurls(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, EdgeData<double>& curl,
                                            VertexData<double>& globalTimeFunction, std::map<int, int>& hashedUsedIsoVals, int numSingularityPairs);

//@debugging
//find edge singularity pairs 
//this method find the max average edge curl by sampling the isolines of the time function 
std::vector<std::tuple<std::pair<int, int>, double>> findEdgeSingularityPairsUsingTimeFunctionIsoVals(VertexPositionGeometry& globalGeometry, 
                                                        EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
                                                        EdgeData<double>& curl, VertexData<double>& globalTimeFunction, 
                                                        double stepSize, std::map<int, int>& hashedUsedIsoVals, int numSingularityPairs);


//@clean
//get the vertices (of a curve network), edges (of a curve network) and faces that a particular isovalue of the time function passes through 
//generate isolines for the time function given a specific isoVal
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, std::vector<int>> getIsoLine(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& timeFunction, double isoVal);


//@clean
//find the isoval with max average curl in the face setting
double findIsoValWithMaxFaceCurl(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& globalTimeFunction, 
                                FaceData<double>& curl, std::map<int, int>& hashedUsedIsoVals, FaceData<int>& forbiddenFaces, double stepSize);

void revealCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& vertexMap,  Eigen::SparseMatrix<double, Eigen::RowMajor>& G);

//@debugging 
//find isoval with maximum average edge curl 
double findIsoValWithMaxAvgEdgeCurl(VertexPositionGeometry& globalGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& globalTimeFunction, 
                                    EdgeData<double>& curl, std::map<int, int>& hashedUsedIsoVals, double stepSize);

//@clean 
//the below two functions use the energy min ||\del sigma - \nabla \sigma_{i - 1} / ||\nabla \sigma_{i - 1}|| ||^2
//@out  courseOneFormGrad - gradient of the final course 1-form
std::tuple<CornerData<double>, EdgeData<double>> implCourseHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    VertexData<double>& globalTimeFunction, FaceData<Vector3>& globalTimeFunctionGradientsNormalized,
                                                                    std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, polyscope::SurfaceMesh& psMesh,
                                                                    globalBoundaryConditions& boundaryConditions, double period,
                                                                    Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double, Eigen::RowMajor>& G,
                                                                    FaceData<Vector3>& courseOneFormGrad, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//@clean
//compute course 1-form
std::tuple<HalfedgeData<double>, double> computeCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh);

//@clean
//update forbidden faces 
//in particular, if some isoline passes through a set of faces ensure we don't select another pair of faces 
//on the same isoline
void updateForbiddenFaces(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& timeFunction, double isoVal, FaceData<int>& forbiddenFaces);

//@debugging 
//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

//@debugging
//average vertex curl onto edges
EdgeData<double> computeVertexAveragedEdgeCurl(VertexPositionGeometry& globalGeometry, VertexData<double>& vertexCurl);

//@debugging
//compute distance from unit norm of a per-face vector field
double computeDistanceFromUnitNorm(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& gradients);

//@debugging
std::tuple<HalfedgeData<double>, double> computeVirtualSigma(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh);


//@debugging 
//compute course 1-form
std::tuple<HalfedgeData<double>, double> computeHarmonicCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh);

#endif
