//file includes
#include "knitting_utils.h"
#include "helpers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//run our greedy algorithm to optimally place singularities
//input: the global geometry
//input: the time function 
//input: \omega the first 1-form in our iterations (in theory could be either course or wale)
//output is the singularity positions on our mesh 
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, VertexData<double>& timeFunction, EdgeData<double>& omega,
                                            double period, std::vector<std::pair<int, int>>& edgeMappingsPairs, globalBoundaryConditions& globalBdyConditions);
//render d1 * omega 
void showd1Omega(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Eigen::VectorXd& omega, int numPairs);

//find a pair of singularities with similar d1 * omega ("curl") values 
std::pair<int, int> findSingularityPair(VertexPositionGeometry& geometry, VertexData<double>& timeFunction,  Eigen::VectorXd& omega);

//given a function defined on mesh vertices and a particular isovalue
//returns vertices and edges and edges in a curve network that passes through that isovalue 
//additionally, returns a set of faces that the curve network passes through 
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, std::vector<int>> getTimeFunctionIsoLine(VertexPositionGeometry& geometry, VertexData<double>& timeFunction, double isoVal);

//computes an iterative 1-form
//the goal here is to make this 1-form behave like a harmonic function
Eigen::VectorXd computeIterativeOneForm(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Model& gbModel, std::vector<std::pair<int, int>>& edgeMappingsPairs, globalBoundaryConditions& globalBoundaryConditions);


//compute the gradient (per face) of a 1-form 
FaceData<Vector3> computeOneFormGrad(VertexPositionGeometry& geometry, EdgeData<double> oneForm);