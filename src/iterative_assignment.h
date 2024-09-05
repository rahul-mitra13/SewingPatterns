//file includes
#include "knitting_utils.h"
#include "helpers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//-------------------First strategy--------------------------------//
//Trying to create a harmonic 1-form i.e., min ||\sigma||^2 subject to some constraints 
//This strategy is extremely challenging in the extrensic setting since for every cylindrical component, we have 2x bdy-bdy path constraints  
//Need to implement this strategy in the intrinsic setting which will make things a lot easier 

//run our greedy algorithm to optimally place singularities
//input: the global geometry
//input: the time function 
//input: \omega the first 1-form in our iterations (in theory could be either course or wale)
//output is the singularity positions on our mesh 
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, VertexData<double>& timeFunction, EdgeData<double>& omega,
                                            double period, std::map<int, int>& vertexMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, globalBoundaryConditions& globalBdyConditions);
//render d1 * omega 
void showd1Omega(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Eigen::VectorXd& omega, int numPairs);

//find a pair of singularities with similar d1 * omega ("curl") values. Similar values but opposite sign
std::pair<int, int> findSingularityPair(VertexPositionGeometry& geometry, VertexData<double>& timeFunction,  Eigen::VectorXd& omega);

//given a function defined on mesh vertices and a particular isovalue
//returns vertices and edges and edges in a curve network that passes through that isovalue 
//additionally, returns a set of faces that the curve network passes through 
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, std::vector<int>> getTimeFunctionIsoLine(VertexPositionGeometry& geometry, VertexData<double>& timeFunction, double isoVal);

//computes an iterative 1-form
//the goal here is to make this 1-form behave like a harmonic function
Eigen::VectorXd computeIterativeOneForm(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Model& gbModel, std::vector<std::pair<int, int>>& edgeMappingsPairs);

//compute the gradient (per face) of a 1-form 
FaceData<Vector3> computeOneFormGrad(VertexPositionGeometry& geometry, EdgeData<double>& oneForm);

//----------------End of first strategy---------------------------//

//------------------Second strategy------------------------------//
//Here we remove the integrability constraint i.e., d1\sigma = 0 and try and greedily match max/min d1\sigma's of similar magnitude but opposite sign
//In this setting, we're working with the co-differential 
//then we re-optimize by fixing those singular faces 
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, Model& model, 
                                                VertexData<double>& gluedTimeFunction, std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, std::vector<bool>& orientations);

std::pair<int, int> findSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction, const Eigen::VectorXd& sigmaTilde);



//------------------Third strategy------------------------------//
//just sample isolines at equal spacing and introduce singularities at max/min d1\omega or d1\sigma at every specific isonline
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, VertexData<double>& gluedTimeFunction, 
                                            const Eigen::VectorXd& curl, std::map<int, int>& vertexMap);


//------------------------------SOME IDEAS TO TEST---------------------------//
//Visualizing the curl of a face-based vector field 
//Eq. 9 from De Goes SIGGRAPH course, "Vector Field Design"
VertexData<int> computeCurlOnVertex(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, Model& model, 
                                                VertexData<double>& gluedTimeFunction, std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, std::vector<bool>& orientations);