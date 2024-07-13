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

using namespace geometrycentral;
using namespace geometrycentral::surface;

//vertex mappings from txt file (JUST CALL PYTHON HERE UGHHH)
std::vector<std::pair<int, int>> vertexMappingsPairs;
std::vector<std::pair<int, int>> edgeMappingsPairs;
//build an index map from vertices in the original mesh to vertices in the glued mesh 
std::map<int, int> indexMap;
//one-ring map for vertices in the glued mesh for performance 
std::map<int, std::vector<Halfedge>> gluedOneRingMap;
//build a glued mesh to make our life a little easier for some procedure 
SurfaceMesh * gluedMesh;

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
  //jacket
  // std::vector<int> zeroVertices = {86, 87, 88, 89, 90, 91, 92, 93, 94,
  //                                 151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
  //                                 1, 2, 3, 4, 5, 6, 7, 8};
  // std::vector<int> oneVertices = {113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
  //                                 178, 179, 180, 181, 182, 183, 184, 185, 186, 187};

  VertexData<double> timeFunction = computeTimeFunction(*globalGeometry, vertexMappingsPairs, globalBdyConditions);
  globalPSMesh->addVertexScalarQuantity("time function", timeFunction);
  FaceData<Vector3> timeFunctionGradient = computeTimeFunctionFaceGrad(*globalGeometry, timeFunction);
  globalPSMesh -> addFaceVectorQuantity("gradient", timeFunctionGradient);
  EdgeData<double> omegaCourse = computeMatchingOneForm(*globalGeometry, 0, timeFunctionGradient, edgeMappingsPairs);
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
  EdgeData<double> sigmaCourse = computeOneForm(*globalGeometry, modelCourse);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma course", sigmaCourse, orientations);
  CornerData<double> stripeValuesSigmaCourse;
  FaceData<int> stripeIndicesSigmaCourse;
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(*globalGeometry, sigmaCourse, period, *gluedMesh, indexMap, 
                                                                  vertexMappingsPairs, edgeMappingsPairs, gluedOneRingMap);
  std::tie(positionsCourse, edgesCourse) = generateIsoLines(*globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
  polyscope::registerCurveNetwork("course stripe patterns", positionsCourse, edgesCourse);
  

  //set up the optimization model for the wale direction 
  Model modelWale;
  modelWale.setPeriod(period);
  std::vector<double> modelMatchingTermsWale;
  EdgeData<double> omegaWale = computeMatchingOneForm(*globalGeometry, 1, timeFunctionGradient, edgeMappingsPairs);
  for (Edge e : globalMesh -> edges()){
    modelMatchingTermsWale.push_back(omegaWale[e]);
  }
  globalPSMesh -> addOneFormTangentVectorQuantity("omega wale", omegaWale, orientations);
  modelWale.setMatchingTerms(modelMatchingTermsWale);
  modelWale.setEdgeMappingsPairs(edgeMappingsPairs);
  modelWale.setWaleBdyPathConstraints(globalBdyConditions.waleBdyPathConstraints);
  EdgeData<double> sigmaWale = computeOneForm(*globalGeometry, modelWale);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma wale", sigmaWale, orientations);
  CornerData<double> stripeValuesSigmaWale;
  FaceData<int> stripeIndicesSigmaWale;
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(stripeValuesSigmaWale, stripeIndicesSigmaWale) = computeStripeValuesFromOneForm(*globalGeometry, sigmaWale, period, *gluedMesh, indexMap, 
                                                              vertexMappingsPairs, edgeMappingsPairs, gluedOneRingMap);
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, stripeValuesSigmaWale, stripeIndicesSigmaWale, period);
  polyscope::registerCurveNetwork("wale stripe patterns", positionsWale, edgesWale);

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
  gluedMesh = createGluedSurfaceMesh(*globalGeometry, vertexMappingsPairs, indexMap, gluedOneRingMap);
  //render the stitched vertices
  renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;
  polyscope::show();

  return EXIT_SUCCESS; 
}