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
//build a map from edges in the orignal mesh to edges in the glued
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
float period = 10;
//threshold for constraining wale boundary edges 
//I don't really like this and need to figure out a better way of doing this
float threshold = 0.6;

//set the permutation of edges and orientations for 1-form viz
std::vector<size_t> perm;
std::vector<bool> orientations;

//here we will do as much processing as possible directly on the glued together mesh 
void showStripePatterns(){
  
  //time function on the glued mesh 
  VertexData<double> timeFunctionGlued = computeTimeFunction(*gluedELG, globalBdyConditions);
  //time function on the global mesh 
  VertexData<double> timeFunctionGlobal = convertGluedToGlobalVertexFunction(*globalGeometry, *gluedELG, timeFunctionGlued, vertexMap);
  globalPSMesh -> addVertexScalarQuantity("time function", timeFunctionGlobal);
  //gradient on the glued/global mesh
  //note that faces have a 1-to-1 mapping from global to glued setting
  FaceData<Vector3> timeFunctionGradientGlobal = computeTimeFunctionFaceGrad(*globalGeometry, timeFunctionGlobal);
  globalPSMesh -> addFaceVectorQuantity("gradient", timeFunctionGradientGlobal);
  //find boundary edges in the wale direction 
  //doing this here cause the gradient gets rotated later
  std::vector<int> waleBdyEdges = getWaleBdyEdgesInGluedMesh(*globalGeometry, *gluedELG, timeFunctionGradientGlobal, edgeMap, threshold, *globalPSMesh);
  //set up the course optimization 
  Model modelCourse; 
  //compute matching one form on the global mesh in the course direction
  EdgeData<double> omegaCourseGlobal = computeMatchingOneForm(*globalGeometry, *globalPSMesh, 0, timeFunctionGradientGlobal, edgeMappingsPairs); 
  //matching one form on the glued mesh in the course direction
  EdgeData<double> omegaCourseGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedELG, omegaCourseGlobal, edgeMap);
  modelCourse.setPeriod(period);
  modelCourse.setBdyEdges(globalBdyConditions.courseBdyEdges);
  std::vector<double> modelMatchingTermsCourse; 
  for (Edge e : (gluedELG -> mesh).edges()){
    modelMatchingTermsCourse.push_back(omegaCourseGlued[e]);
  }
  modelCourse.setMatchingTerms(modelMatchingTermsCourse);
  //sigma in the glued mesh setting 
  EdgeData<double> sigmaCourseGlued = computeOneForm(*gluedELG, modelCourse);
  //global data
  CornerData<double> stripeValuesSigmaCourse;
  //global data
  FaceData<int> stripeIndicesSigmaCourse;
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(*globalGeometry, *gluedELG, sigmaCourseGlued, period);
  std::tie(positionsCourse, edgesCourse) = generateIsoLines(*globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
  auto courseStripes = polyscope::registerCurveNetwork("course stripe patterns", positionsCourse, edgesCourse);
  courseStripes -> setRadius(0.001);

  //set up the wale optimization model 
  Model modelWale; 
  modelWale.setPeriod(period);
  //compute matching one form on the global mesh in the wale direction 
  EdgeData<double> omegaWaleGlobal = computeMatchingOneForm(*globalGeometry, *globalPSMesh, 1, timeFunctionGradientGlobal, edgeMappingsPairs);
  //globalPSMesh -> addOneFormTangentVectorQuantity("omega wale", omegaWaleGlobal, orientations);
  EdgeData<double> omegaWaleGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedELG, omegaWaleGlobal, edgeMap);
  std::vector<double> modelMatchingTermsWale; 
  for (Edge e : (gluedELG -> mesh).edges()){
    modelMatchingTermsWale.push_back(omegaWaleGlued[e]);
  }
  //globalPSMesh -> addOneFormTangentVectorQuantity("omega wale ", omegaWaleGlobal, orientations);
  modelWale.setMatchingTerms(modelMatchingTermsWale);
  modelWale.setWaleBdyPathConstraints(globalBdyConditions.waleBdyPathConstraints);
  modelWale.setBdyEdges(waleBdyEdges);
  //sigma in the glued mesh setting 
  EdgeData<double> sigmaWaleGlued = computeOneForm(*gluedELG, modelWale);
  EdgeData<double> sigmaWaleGlobal = convertGluedToGlobalEdgeFunction(*globalGeometry, *gluedELG, sigmaWaleGlued, edgeMap);
  //globalPSMesh -> addOneFormTangentVectorQuantity("sigma wale", sigmaWaleGlobal, orientations);
  //visualize the differences
  EdgeData<double> difference(*globalMesh);
  for (Edge e : globalMesh -> edges()){
    difference[e] = pow(omegaWaleGlobal[e] - sigmaWaleGlobal[e], 2.0);
  }
  globalPSMesh -> addEdgeScalarQuantity("difference", difference);
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

  //greedily placed singularities 
  // FaceData<int> greedySingularities = getGreedySingularityPositions(*globalGeometry, *globalPSMesh, timeFunctionGlobal, omegaCourseGlobal, period, vertexMap, edgeMappingsPairs, globalBdyConditions);
  // globalPSMesh -> addFaceScalarQuantity("Greedily placed singularities", greedySingularities);
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