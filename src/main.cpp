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
std::map<int, int> indexMap;
//build a map from edges in the orignal mesh to edges in the glued mesh 
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

//set the permutation of edges and orientations for 1-form viz
std::vector<size_t> perm;
std::vector<bool> orientations;

//render stripe patterns over the surface
void showStripePatterns(){ 
  VertexData<double> timeFunction = computeTimeFunction(*globalGeometry, vertexMappingsPairs, globalBdyConditions, indexMap);
  globalPSMesh->addVertexScalarQuantity("time function", timeFunction);
  FaceData<Vector3> timeFunctionGradient = computeTimeFunctionFaceGrad(*globalGeometry, timeFunction);
  globalPSMesh -> addFaceVectorQuantity("gradient", timeFunctionGradient);
  EdgeData<double> omegaCourse = computeMatchingOneForm(*globalGeometry, *globalPSMesh, 0, timeFunctionGradient, edgeMappingsPairs);
  globalPSMesh -> addOneFormTangentVectorQuantity("omega course", omegaCourse, orientations);
  
  //set up the optimization model for the course direction
  Model modelCourse;
  modelCourse.setBdyEdges(globalBdyConditions.courseBdyEdges);
  modelCourse.setPeriod(period);
  std::vector<double> modelMatchingTermsCourse; 
  for (Edge e : globalMesh -> edges()){
    modelMatchingTermsCourse.push_back(omegaCourse[e]);
  }
  modelCourse.setMatchingTerms(modelMatchingTermsCourse);
  modelCourse.setEdgeMappingsPairs(edgeMappingsPairs);
  //sigma course over the glued edge length geometry
  EdgeData<double> sigmaCourseELG(*globalMesh);
  sigmaCourseELG = computeOneForm(*globalGeometry, *gluedELG, modelCourse, edgeMap, *globalPSMesh);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma course", sigmaCourseELG, orientations);
  CornerData<double> stripeValuesSigmaCourse;
  FaceData<int> stripeIndicesSigmaCourse;
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(*globalGeometry, sigmaCourseELG, period, gluedELG->mesh, indexMap, 
                                                                  vertexMappingsPairs, edgeMappingsPairs, gluedOneRingMap);
  std::tie(positionsCourse, edgesCourse) = generateIsoLines(*globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
  auto courseStripes = polyscope::registerCurveNetwork("course stripe patterns", positionsCourse, edgesCourse);
  courseStripes -> setRadius(0.004);
  

  //set up the optimization model for the wale direction 
  Model modelWale;
  modelWale.setPeriod(period);
  modelWale.setBdyEdges(globalBdyConditions.waleBdyEdges);
  std::vector<double> modelMatchingTermsWale;
  EdgeData<double> omegaWale = computeMatchingOneForm(*globalGeometry, *globalPSMesh, 1, timeFunctionGradient, edgeMappingsPairs);
  for (Edge e : globalMesh -> edges()){
    modelMatchingTermsWale.push_back(omegaWale[e]);
  }
  globalPSMesh -> addOneFormTangentVectorQuantity("omega wale", omegaWale, orientations);
  // globalPSMesh -> addEdgeScalarQuantity("omega wale", omegaWale);
  modelWale.setMatchingTerms(modelMatchingTermsWale);
  modelWale.setEdgeMappingsPairs(edgeMappingsPairs);
  modelWale.setWaleBdyPathConstraints(globalBdyConditions.waleBdyPathConstraints);
  EdgeData<double> sigmaWale = computeOneForm(*globalGeometry, modelWale, *globalPSMesh);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma wale", sigmaWale, orientations);
  //view sigma as an edge scalar
  // globalPSMesh -> addEdgeScalarQuantity("sigma wale", sigmaWale);
  CornerData<double> stripeValuesSigmaWale;
  FaceData<int> stripeIndicesSigmaWale;
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(stripeValuesSigmaWale, stripeIndicesSigmaWale) = computeStripeValuesFromOneForm(*globalGeometry, sigmaWale, period, gluedELG->mesh, indexMap, 
                                                              vertexMappingsPairs, edgeMappingsPairs, gluedOneRingMap);
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, stripeValuesSigmaWale, stripeIndicesSigmaWale, period);
  auto waleStripes = polyscope::registerCurveNetwork("wale stripe patterns", positionsWale, edgesWale);
  waleStripes -> setRadius(0.004);
 
  //for greedily placing singularities 
  FaceData<int> singPositions(*globalMesh);
  singPositions = getGreedySingularityPositions(*globalGeometry, *globalPSMesh, timeFunction, omegaCourse, period, edgeMappingsPairs, globalBdyConditions);
  globalPSMesh -> addFaceScalarQuantity("greedy singularities", singPositions);
  
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::InputFloat("1-form period", &period);
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
  globalBdyConditions = parseJson(*globalGeometry, data);
  gluedELG = createGluedEdgeLengthGeometry(*globalGeometry, vertexMappingsPairs, indexMap, edgeMap, gluedOneRingMap);
  //render the stitched vertices
  renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;
  polyscope::show();

  return EXIT_SUCCESS; 
}