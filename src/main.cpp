//this is to appropriately parse the GarmentCode specification
#include <Python.h>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/edge_length_geometry.h"



#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"
#include "nlohmann/json.hpp"
#include "imgui.h"

//file includes
#include "knitting_utils.h"
#include "stripe_patterns_helpers.h"
#include "iterative_assignment.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//vertex mappings from txt file (between panels)
std::vector<std::pair<int, int>> vertexMappingsPairs;
//edge mappings from txt file
std::vector<std::pair<int, int>> edgeMappingsPairs;
//build an index map from vertices in the original mesh to vertices in the glued mesh 
std::map<int, int> vertexMap;
//build an index map from edges in the orignal mesh to edges in the glued
std::map<int, int> edgeMap;
//one-ring map for vertices in the glued mesh for performance 
std::map<int, std::vector<Halfedge>> gluedOneRingMap;
//build a glued edge length geoemtry to make our life a little easier for some procedures
EdgeLengthGeometry *gluedELG;

std::unique_ptr<ManifoldSurfaceMesh> globalMesh;
std::unique_ptr<VertexPositionGeometry> globalGeometry;
polyscope::SurfaceMesh *globalPSMesh;
//global boundary conditions
globalBoundaryConditions globalBdyConditions;


//1-form optimization period
float period = 1;
//threshold for constraining wale boundary edges 
//I don't really like this and need to figure out a better way of doing this
float threshold = 0.6;

//set the permutation of edges and orientations for 1-form viz
std::vector<size_t> perm;
std::vector<bool> orientations;

//here we will do as much processing as possible directly on the glued together mesh 
void showStripePatterns(){
  
  //require the DEC operators
  gluedELG->requireDECOperators();
  //time function on the glued mesh 
  VertexData<double> timeFunctionGlued = computeTimeFunction(*gluedELG, globalBdyConditions);
  //time function on the global mesh 
  VertexData<double> timeFunctionGlobal = convertGluedToGlobalVertexFunction(*globalGeometry, *gluedELG, timeFunctionGlued, vertexMap);
  globalPSMesh -> addVertexScalarQuantity("time function", timeFunctionGlobal);
  //gradient on the glued/global mesh
  //note that faces have a 1-to-1 mapping from global to glued setting
  FaceData<Vector3> timeFunctionGradientGlobal = computeTimeFunctionFaceGrad(*globalGeometry, timeFunctionGlobal);
  globalPSMesh -> addFaceVectorQuantity("gradient (unnormalized)", timeFunctionGradientGlobal);
  //find boundary edges in the wale direction 
  //doing this here cause the gradient gets rotated later
  std::vector<int> waleBdyEdges = getWaleBdyEdgesInGluedMesh(*globalGeometry, *gluedELG, timeFunctionGradientGlobal, edgeMap, threshold, *globalPSMesh);
  
  //set up the course optimization 
  Model modelCourse;
  modelCourse.setIntegrabilityConstraint(true); 
  modelCourse.useFaceDifferenceViz = true;
  //set the face gradients 
  std::vector<std::array<double, 3>> modelFaceGradientsCourse;
  for (Face f : globalMesh->faces()){
    //normalize the gradients first
    timeFunctionGradientGlobal[f] = timeFunctionGradientGlobal[f].normalize();
    std::array<double, 3> gradient = {timeFunctionGradientGlobal[f][0], timeFunctionGradientGlobal[f][1], timeFunctionGradientGlobal[f][2]}; 
    modelFaceGradientsCourse.push_back(gradient);
  }
  modelCourse.setFaceGradients(modelFaceGradientsCourse);
  modelCourse.setPeriod(period);
  modelCourse.setBdyEdges(globalBdyConditions.courseBdyEdges);
  //visualize the normalized gradients 
  globalPSMesh -> addFaceVectorQuantity("gradient (normalized)", timeFunctionGradientGlobal);

  EdgeData<double> omegaCourseGlobal = computeMatchingOneForm(*globalGeometry, *globalPSMesh, 0, timeFunctionGradientGlobal, edgeMappingsPairs);
  EdgeData<double> omegaCourseGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedELG, omegaCourseGlobal, edgeMap);
  Eigen::Map<Eigen::VectorXd> omegaEig(omegaCourseGlued.raw().data(), (gluedELG->mesh).nEdges());
  std::vector<double> modelMatchingTermsCourse(omegaEig.data(), omegaEig.data() + omegaEig.rows());
  modelCourse.setMatchingTerms(modelMatchingTermsCourse);
  Eigen::VectorXd d1Omega = gluedELG->d1 * omegaEig;
  globalPSMesh->addFaceScalarQuantity("d1(omega0)", d1Omega);
  //visualize the 1-form \omega using Whitney interpolation 
  globalPSMesh->addOneFormTangentVectorQuantity("omega0(Whitney)", omegaCourseGlobal, orientations);
  //sigma in the glued mesh setting 
  EdgeData<double> sigmaCourseGlued = computeOneForm(*globalGeometry, *gluedELG, modelCourse, vertexMap, edgeMap, *globalPSMesh);
  EdgeData<double> sigmaCourseGlobal = convertGluedToGlobalEdgeFunction(*globalGeometry, *gluedELG, sigmaCourseGlued, edgeMap);
  
  //Eigen::Map<Eigen::VectorXd> sigmaEig(sigmaCourseGlued.raw().data(), (gluedELG->mesh).nEdges());
  //Eigen::VectorXd d1Sigma = gluedELG->d1 * omegaEig;
  //globalPSMesh->addFaceScalarQuantity("d1 sigma", d1Sigma);

  //global data
  CornerData<double> stripeValuesSigmaCourse;
  //global data
  FaceData<int> stripeIndicesSigmaCourse;
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(*globalGeometry, *gluedELG, sigmaCourseGlued, period);
  std::tie(positionsCourse, edgesCourse) = generateIsoLines(*globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
  auto courseStripes = polyscope::registerCurveNetwork("course stripe patterns no sings (sigma)", positionsCourse, edgesCourse);
  courseStripes -> setRadius(0.001);
  
  //set up the wale optimization model
  /**  
  Model modelWale; 
  modelWale.setIntegrabilityConstraint(true);
  modelWale.setPeriod(period);
  std::vector<std::array<double, 3>> modelFaceGradientsWale;
  globalGeometry->requireFaceNormals();
  for (Face f : globalMesh->faces()){
    //normalize the gradients first
    timeFunctionGradientGlobal[f] = timeFunctionGradientGlobal[f].normalize();
    timeFunctionGradientGlobal[f] = timeFunctionGradientGlobal[f].rotateAround(globalGeometry->faceNormals[f], PI/2.);
    std::array<double, 3> gradient = {timeFunctionGradientGlobal[f][0], timeFunctionGradientGlobal[f][1], timeFunctionGradientGlobal[f][2]}; 
    modelFaceGradientsWale.push_back(gradient);
  }
  modelWale.setFaceGradients(modelFaceGradientsWale);
  modelWale.setWaleBdyPathConstraints(globalBdyConditions.waleBdyPathConstraints);
  modelWale.setBdyEdges(waleBdyEdges);
  //sigma in the glued mesh setting 
  EdgeData<double> sigmaWaleGlued = computeOneForm(*globalGeometry, *gluedELG, modelWale, vertexMap, *globalPSMesh);
  EdgeData<double> sigmaWaleGlobal = convertGluedToGlobalEdgeFunction(*globalGeometry, *gluedELG, sigmaWaleGlued, edgeMap);
  //global data
  CornerData<double> stripeValuesSigmaWale;
  //global data
  FaceData<int> stripeIndicesSigmaWale;
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(stripeValuesSigmaWale, stripeIndicesSigmaWale) = computeStripeValuesFromOneForm(*globalGeometry, *gluedELG, sigmaWaleGlued, period);
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, stripeValuesSigmaWale, stripeIndicesSigmaWale, period);
  auto waleStripes = polyscope::registerCurveNetwork("wale stripe patterns", positionsWale, edgesWale);
  waleStripes -> setRadius(0.001);
  */

  //greedily placed singularities 
  FaceData<int> greedySingularities = getGreedySingularityPositions(*globalGeometry, *gluedELG, *globalPSMesh, modelCourse, timeFunctionGlued, vertexMap, edgeMap, orientations);
  //globalPSMesh -> addFaceScalarQuantity("greedy singularities", greedySingularities)

  //-----------------------------DEBUGGING STUFF-----------------------------//

  //visualize ||\sigma - \omega_c||^2 dubject to d1 constraints everywhere
  //numerous hard-coded values here 
  Model modelCourseDebug = modelCourse;
  modelCourseDebug.useEdgeAveraging = true; 
  std::vector<int> faceIndices(globalMesh -> nFaces(), 0);
  
  //singularity indices on the bent cylinder 
  faceIndices[7] = 1; 
  faceIndices[10] = -1;

  globalPSMesh -> addFaceScalarQuantity("singular faces", faceIndices);
  modelCourseDebug.setFaceIndices(faceIndices);
  EdgeData<double> omegaDebugGlued = computeOneForm(*globalGeometry, *gluedELG, modelCourseDebug, vertexMap, edgeMap, *globalPSMesh);
  EdgeData<double> omegaDebugGlobal = convertGluedToGlobalEdgeFunction(*globalGeometry, *gluedELG, omegaDebugGlued, edgeMap);
  //global data
  CornerData<double> stripeValuesOmegaDebug;
  //global data
  FaceData<int> stripeIndicesOmegaDebug;
  std::vector<Vector3> positionsDebug;
  std::vector<std::array<int, 2>> edgesDebug;
  std::tie(stripeValuesOmegaDebug, stripeIndicesOmegaDebug) = computeStripeValuesFromOneForm(*globalGeometry, *gluedELG, omegaDebugGlued, period);
  std::tie(positionsDebug, edgesDebug) = generateIsoLines(*globalGeometry, stripeValuesOmegaDebug, stripeIndicesOmegaDebug, period);
  auto debugStripes = polyscope::registerCurveNetwork("debug stripe patterns sings (omega)", positionsDebug, edgesDebug);
  debugStripes -> setRadius(0.001);
  globalPSMesh -> addOneFormTangentVectorQuantity("omega with sings (Whitney)", omegaDebugGlobal, orientations);
    
  //visualize ||\delta \sigma - \frac{\nabla h}{||\nabla h||}||^2
  modelCourseDebug.useEdgeAveraging = false; 
  modelCourseDebug.useFaceDifferenceViz = true;
  //modelCourseDebug.usePsuedoHarmViz = true;
  EdgeData<double> sigmaDebugGlued = computeOneForm(*globalGeometry, *gluedELG, modelCourseDebug, vertexMap, edgeMap,*globalPSMesh);
  EdgeData<double> sigmaDebugGlobal = convertGluedToGlobalEdgeFunction(*globalGeometry, *gluedELG, sigmaDebugGlued, edgeMap);
  //global data
  CornerData<double> stripeValuesSigmaDebug;
  //global data
  FaceData<int> stripeIndicesSigmaDebug;
  std::tie(stripeValuesSigmaDebug, stripeIndicesSigmaDebug) = computeStripeValuesFromOneForm(*globalGeometry, *gluedELG, sigmaDebugGlued, period);
  std::tie(positionsDebug, edgesDebug) = generateIsoLines(*globalGeometry, stripeValuesSigmaDebug, stripeIndicesSigmaDebug, period);
  debugStripes = polyscope::registerCurveNetwork("debug stripe patterns sings (sigma)", positionsDebug, edgesDebug);
  debugStripes -> setRadius(0.001);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma with sings (Whitney)", sigmaDebugGlobal, orientations);

  //don't have ||W(\sigma) - \frac{\nabla h}{||\nabla h||}||^2 - would have to implement Whitney interpolant in this setting 


}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::InputFloat("1-form period", &period);
  //there needs to be a better way to constrain wale edges
  ImGui::InputFloat("Threshold", &threshold);

  if (ImGui::Button("Show Stripe Patterns")){
    showStripePatterns();
  }
}

int main(int argc, char **argv) {
  polyscope::init();
  std::ifstream jsonFile(argv[1]);
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  //run sanity checks
  std::tie(globalMesh, globalGeometry) = readManifoldSurfaceMesh(data["model_path"]);
  if (!(globalMesh -> isManifold())){
    std::cout << "Error: Mesh is not manifold" << std::endl;
    throw std::exception();
  }
  if (!(globalMesh -> isOriented())){
    std::cout << "Error: Meshing is not oriented" << std::endl;
    throw std::exception();
  }
  if (!(globalMesh -> isTriangular())){
    std::cout << "Error: Mesh is not triangular" << std::endl;
    throw std::exception();
  }
  globalPSMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(data["model_path"]), globalGeometry->inputVertexPositions, globalMesh -> getFaceVertexList());
  vertexMappingsPairs = buildPairOfStitchedVerticesFromFile(data["vertex_mappings"]);
  edgeMappingsPairs = buildPairOfStitchedEdges(*globalGeometry, vertexMappingsPairs);
  //set up the orientations for the 1-form viz while we're here
  for (Edge e : globalMesh->edges()){
    if (e.halfedge().tailVertex().getIndex() < e.halfedge().tipVertex().getIndex()){
      perm.push_back(e.getIndex());
      orientations.push_back(true);
    }
    else{
      perm.push_back(e.getIndex());
      orientations.push_back(false);
    }
  }
  //now update the orientations to handle "stitched together edges" which actually represent a single edge 
  //just flip the orientations to make the viz sensible
  for (std::pair<int, int> pair : edgeMappingsPairs){
    orientations[pair.second] = !orientations[pair.second];
  }
  globalPSMesh -> setEdgePermutation(perm);
  //create the glued edge length geometry
  gluedELG = createGluedEdgeLengthGeometry(*globalGeometry, vertexMappingsPairs, vertexMap, edgeMap, gluedOneRingMap);
  //process boundary conditions in the glued mesh setting 
  globalBdyConditions = parseJson(*gluedELG, data, vertexMap, edgeMap);
  //render the stitched vertices
  renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;
  polyscope::show();

  return EXIT_SUCCESS; 
}