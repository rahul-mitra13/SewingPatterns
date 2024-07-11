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
#include "GUI_helpers.h"
#include "stripe_patterns_helpers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

int numPatches;//number of patches

//make a vector to handle all the patches 
//meshes + geometries + polyscope surface meshes
std::vector<std::unique_ptr<ManifoldSurfaceMesh>> meshes;
std::vector<std::unique_ptr<VertexPositionGeometry>> geometries;
std::vector<polyscope::SurfaceMesh *> psMeshes;

//vertex mappings from txt file (JUST CALL PYTHON HERE UGHHH)
std::vector<std::pair<int, int>> vertexMappingsPairs;
std::vector<std::pair<int, int>> edgeMappingsPairs;
//also build a map that stores panel name to geometry (can infer the mesh from the geometry)
//std::map<std::string, std::unique_ptr<VertexPositionGeometry>> panelMappings;

std::unique_ptr<ManifoldSurfaceMesh> globalMesh;
std::unique_ptr<VertexPositionGeometry> globalGeometry;
polyscope::SurfaceMesh *globalPSMesh;
SurfaceMesh *gluedMesh;//represents the glued mesh of all the panel 
std::map<int, int> originalMeshVertexIndexToGluedMeshIndex;
std::map<int, int> originalMeshEdgeIndexToGluedMeshIndex;
std::vector<double> gluedMeshEdgeLengths;

//boundary conditions for each patch 
std::vector<PatchBoundaryConditions> bdyConditions;

//path to json file
std::string jsonFilePath;
bool parseUsingJSon = false;


// Striping frequency
float constantFreq = 10.0;
//1-form optimization period
float period = 0.1;


//testing vertex selection UI (only applies to multiple panels)
void showBoundaryVertexSelectionUI(){
  //do this for every patch
  for (int i = 0; i < numPatches; i++){
    if (parseUsingJSon){
      bdyConditions[i] = getAndRenderBoundaryInfoFromJson(*geometries[i], *psMeshes[i], jsonFilePath);
    }
    else{
      //select course boundary conditions
      std::tie(bdyConditions[i].courseStartBoundaryVertices, bdyConditions[i].courseStartBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 0);
      std::tie(bdyConditions[i].courseEndBoundaryVertices, bdyConditions[i].courseEndBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 1);
      //select wale boundary conditions
      //std::tie(bdyConditions[i].waleStartBoundaryVertices, bdyConditions[i].waleStartBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 0);
      //std::tie(bdyConditions[i].waleEndBoundaryVertices, bdyConditions[i].waleEndBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 1);
    }
  }
}


//render stripe patterns over the surface
void showStripePatterns(){ 

  //pants
  std::vector<int> zeroVertices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                  99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
                                  148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
                                  49, 50, 51, 52, 53, 54, 55, 56, 57, 58};
  std::vector<int> oneVertices = {27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
                                  117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
                                  175, 176, 177, 178, 179, 180, 181, 182, 183, 184,
                                  67, 68, 69, 70, 71, 72, 73, 74, 75, 76};

  //jumpsuit 
  // std::vector<int> zeroVertices = {232, 233, 234, 235, 236, 237, 238, 239, 240, 241,
  //                                 346, 347, 348, 349, 350, 351, 352, 353, 354, 355,
  //                                 404, 405, 406, 407, 408, 409, 410, 411, 412, 413,
  //                                 288, 289, 290, 291, 292, 293, 294, 295, 296, 297};
  // std::vector<int> oneVertices = {172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
  //                                 145, 146, 147, 148, 149, 150, 151, 152, 153, 154};

  //jacket
  // std::vector<int> zeroVertices = {86, 87, 88, 89, 90, 91, 92, 93, 94,
  //                                 151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
  //                                 1, 2, 3, 4, 5, 6, 7, 8};
  // std::vector<int> oneVertices = {113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
  //                                 178, 179, 180, 181, 182, 183, 184, 185, 186, 187};

  globalBoundaryConditions boundaryConditions;
  
  for (int index : zeroVertices){
    boundaryConditions.courseStartBoundaryVertices.push_back(globalMesh->vertex(index));
  }
  for (int index : oneVertices){
    boundaryConditions.courseEndBoundaryVertices.push_back(globalMesh->vertex(index));
  }
  
  VertexData<double> timeFunction = computeTimeFunction(*globalGeometry, vertexMappingsPairs, boundaryConditions);
  globalPSMesh->addVertexScalarQuantity("time function", timeFunction);

  FaceData<Vector3> timeFunctionGradient = computeTimeFunctionFaceGrad(*globalGeometry, timeFunction);
  globalPSMesh -> addFaceVectorQuantity("gradient", timeFunctionGradient);

  EdgeData<double> omega_course = computeMatchingOneForm(*globalGeometry, 0, timeFunctionGradient, edgeMappingsPairs);

  std::vector<size_t> perm;
  //set up 1-form viz 
  std::vector<bool> orientations;
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
  globalPSMesh -> addOneFormTangentVectorQuantity("omega course", omega_course, orientations);

  std::vector<int> courseStartEdges = creatEdgeListFromVertexList(*globalGeometry, zeroVertices);
  std::vector<int> courseEndEdges = creatEdgeListFromVertexList(*globalGeometry, oneVertices);

  std::vector<int> bdyEdges;
  std::merge(courseStartEdges.begin(), courseStartEdges.end(),
           courseEndEdges.begin(), courseEndEdges.end(),
           std::back_inserter(bdyEdges));
  
  //set up the optimization model for the course direction
  Model model;
  float period = 10;
  model.setBdyEdges(bdyEdges);
  model.setPeriod(period);
  std::vector<double> modelMatchingTerms; 
  for (Edge e : globalMesh -> edges()){
    modelMatchingTerms.push_back(omega_course[e]);
  }
  model.setMatchingTerms(modelMatchingTerms);
  model.setEdgeMappingsPairs(edgeMappingsPairs);
  EdgeData<double> sigma_course = computeOneForm(*globalGeometry, model, *globalPSMesh);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma course", sigma_course, orientations);

  CornerData<double> stripeValuesSigma;
  FaceData<int> stripeIndicesSigma;
  std::vector<Vector3> positions;
  std::vector<std::array<int, 2>> edges;
  std::tie(stripeValuesSigma, stripeIndicesSigma) = computeStripeValuesFromOneForm(*globalGeometry, sigma_course, period, vertexMappingsPairs, edgeMappingsPairs, *globalPSMesh);
  std::tie(positions, edges) = generateIsoLines(*globalGeometry, stripeValuesSigma, stripeIndicesSigma, period);
  polyscope::registerCurveNetwork("course stripe patterns", positions, edges);
  

  //set up the optimization model for the wale direction 
  Model modelWale;
  std::vector<int> waleBdyEdges = {281, 252, 277, 279, 278, 265, 283, 298, 247,
                                  201, 271, 204, 208, 217, 222, 231, 246, 236,
                                  238, 228, 250, 268, 255, 264, 241, 260, 261, 
                                  68, 67, 74, 75, 71, 52, 59, 63, 90,
                                  14, 9, 2, 6, 24, 86, 46, 34, 29, 
                                  163, 161, 155, 138, 164, 113, 144, 150, 141};
  modelWale.setPeriod(period);
  std::vector<double> modelMatchingTermsWale;
  EdgeData<double> omega_wale = computeMatchingOneForm(*globalGeometry, 1, timeFunctionGradient, edgeMappingsPairs);
  for (Edge e : globalMesh -> edges()){
    modelMatchingTermsWale.push_back(omega_wale[e]);
  }
  globalPSMesh -> addOneFormTangentVectorQuantity("omega wale", omega_wale, orientations);
  modelWale.setMatchingTerms(modelMatchingTermsWale);
  modelWale.setEdgeMappingsPairs(edgeMappingsPairs);
  EdgeData<double> sigma_wale = computeOneForm(*globalGeometry, modelWale, *globalPSMesh);
  globalPSMesh -> addOneFormTangentVectorQuantity("sigma wale", sigma_wale, orientations);
  CornerData<double> stripeValuesSigmaWale;
  FaceData<int> stripeIndicesSigmaWale;
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(stripeValuesSigmaWale, stripeIndicesSigmaWale) = computeStripeValuesFromOneForm(*globalGeometry, sigma_wale, period, vertexMappingsPairs, edgeMappingsPairs, *globalPSMesh);
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
  //std::string path = argv[1];
  jsonFilePath = argv[1];
  if (parseUsingJSon){
    getAndRenderBoundaryInfoFromJson(jsonFilePath, meshes, geometries, bdyConditions);
    psMeshes.resize(meshes.size());
    for (int i = 0; i < geometries.size(); i++){
      psMeshes[i] = polyscope::registerSurfaceMesh(
        "mesh " + std::to_string(i),
        geometries[i]->inputVertexPositions, meshes[i]->getFaceVertexList(),
        polyscopePermutations(*meshes[i]));
    }
  }
  else{
    //run sanity checks (mostly doing this because stripes are failing)
    std::tie(globalMesh, globalGeometry) = readManifoldSurfaceMesh(argv[1]);
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

    globalPSMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(argv[1]), globalGeometry->inputVertexPositions, globalMesh -> getFaceVertexList());
    vertexMappingsPairs = buildPairOfStitchedVerticesFromFile(argv[2]);
    edgeMappingsPairs = buildPairOfStitchedEdges(*globalGeometry, vertexMappingsPairs);
    //render the stitched vertices
    renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  }

  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;
  polyscope::show();

  return EXIT_SUCCESS;
}