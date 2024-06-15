#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"



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

#include <igl/readOBJ.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

int numPatches;//number of patches

//make a vector to handle all the patches 
//meshes + geometries + polyscope surface meshes
std::vector<std::unique_ptr<ManifoldSurfaceMesh>> meshes;
std::vector<std::unique_ptr<VertexPositionGeometry>> geometries;
std::vector<polyscope::SurfaceMesh *> psMeshes;
//boundary conditions for each patch 
std::vector<PatchBoundaryConditions> bdyConditions;
//path to json file
std::string jsonFilePath;
bool parseUsingJSon = true;


// Striping frequency
float constantFreq = 10.0;
//1-form optimization period
float period = 0.1;


//testing vertex selection UI 
void showBoundaryVertexSelectionUI(){
  //do this for every patch
  for (int i = 0; i < numPatches; i++){
    if (parseUsingJSon){
      bdyConditions[i] = getAndRenderBoundaryInfoFromJson(*geometries[i], *psMeshes[i], jsonFilePath);
      std::cout << "Size of course start vertices " << bdyConditions[i].courseStartBoundaryVertices.size() << std::endl;
      std::cout << "Size of course end boundary vertices " << bdyConditions[i].courseEndBoundaryVertices.size() << std::endl;
    }
    else{
      //select course boundary conditions
      std::tie(bdyConditions[i].courseStartBoundaryVertices, bdyConditions[i].courseStartBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 0);
      std::tie(bdyConditions[i].courseEndBoundaryVertices, bdyConditions[i].courseEndBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 1);
      //select wale boundary conditions
      std::tie(bdyConditions[i].waleStartBoundaryVertices, bdyConditions[i].waleStartBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 0);
      std::tie(bdyConditions[i].waleEndBoundaryVertices, bdyConditions[i].waleEndBoundaryEdges) = getAndRenderUserSpecifiedBoundaryInfo(*geometries[i], *psMeshes[i], 1);
    }
  }
}

//render stripe patterns over the surface
void showStripePatterns() {
  
  for (int i = 0; i < numPatches; i++){
    std::vector<Vertex> zeroVertices, oneVertices;
    zeroVertices = bdyConditions[i].courseStartBoundaryVertices;
    oneVertices = bdyConditions[i].courseEndBoundaryVertices;
    VertexData<double> timeFunction = solveLaplace(*geometries[i], zeroVertices, oneVertices);
    psMeshes[i] -> addVertexScalarQuantity("time function_" + std::to_string(i), timeFunction);
    FaceData<Vector3> faceGrads = computeTimeFunctionFaceGrad(*geometries[i], timeFunction);
    psMeshes[i] -> addFaceVectorQuantity("face gradients_" + std::to_string(i), faceGrads);

    //Compute knoppel course stripe patterns
    /**
    VertexData<Vector3> courseVertexValuedField = computeVertexValuedField(*geometries[i], timeFunction, 0);
    VertexData<Vector2> lineFieldCourse = vertexDirectionField(*geometries[i], courseVertexValuedField);
    VertexData<double> frequencies(*meshes[i], constantFreq);
    CornerData<double> periodicFunc;
    FaceData<int> zeroIndices;
    FaceData<int> branchIndices;
    std::tie(periodicFunc, zeroIndices, branchIndices) =
        computeStripePattern(*geometries[i], frequencies, lineFieldCourse);
    psMeshes[i] -> addFaceScalarQuantity("direction field singularities_" + std::to_string(i), branchIndices);
    psMeshes[i] -> addFaceScalarQuantity("stripe singular faces_" + std::to_string(i), zeroIndices);
    // Extract isolines
    std::vector<Vector3> isolineVerts;
    std::vector<std::array<size_t, 2>> isolineEdges;
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometries[i], periodicFunc, zeroIndices, branchIndices, lineFieldCourse, false);
    //polyscope::registerCurveNetwork("course stripe patterns_" + std::to_string(i), isolineVerts, isolineEdges);

    //Compute wale stripe patterns using knoppel's method
    VertexData<Vector3> waleVertexValuedField = computeVertexValuedField(*geometries[i], timeFunction, PI/2);
    psMeshes[i] -> addVertexVectorQuantity("wale vertex valued field_" + std::to_string(i), waleVertexValuedField);
    VertexData<Vector2> lineFieldWale = vertexDirectionField(*geometries[i], waleVertexValuedField);
    std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometries[i], frequencies, lineFieldWale);
    // Extract isolines
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometries[i], periodicFunc, zeroIndices, branchIndices, lineFieldWale, false);
    //polyscope::registerCurveNetwork("wale stripe patterns_" + std::to_string(i), isolineVerts, isolineEdges);
    */

    //generate course stripe patterns using 1-form approach 
    // std::vector<Edge> startEdges, endEdges;
    // startEdges = bdyConditions[i].courseStartBoundaryEdges;
    // endEdges = bdyConditions[i].courseEndBoundaryEdges;
    // EdgeData<double>  matchingOneForm = computeMatchingOneForm(*geometries[i], 0, faceGrads);
    // //set up the optimization model 
    // Model model;
    // //put all the boundary edges into a single vector 
    // std::vector<int> modelBdyEdges;
    // //put the edge terms we're trying to match in a vector
    // std::vector<double> modelMatchingTerms;

    // for (Edge e : startEdges) modelBdyEdges.push_back(e.getIndex());
    // for (Edge e : endEdges) modelBdyEdges.push_back(e.getIndex());
    // for (Edge e : meshes[i]->edges()) modelMatchingTerms.push_back(matchingOneForm[e]);

    // model.setBdyEdges(modelBdyEdges);
    // model.setMatchingTerms(modelMatchingTerms);
    // model.setPeriod(period);

    // EdgeData<double> courseOneForm = computeOneForm(*geometries[i], model);

    // CornerData<double> stripeValuesSigma;
    // FaceData<int> stripeIndicesSigma;  
    // std::vector<Vector3> positions;
    // std::vector<std::array<int, 2>> edges;

    // //compute course stripe patterns from 1-form
    // std::tie(stripeValuesSigma, stripeIndicesSigma) = computeStripeValuesFromOneForm(*geometries[i], courseOneForm, period);
    // std::tie(positions, edges) = generateIsoLines(*geometries[i], stripeValuesSigma, stripeIndicesSigma, period);
    // polyscope::registerCurveNetwork("course stripe patterns_" + std::to_string(i), positions, edges);

    //generate wale stripe patterns using 1-form approach 
    // startEdges = bdyConditions[i].waleStartBoundaryEdges;
    // endEdges = bdyConditions[i].waleEndBoundaryEdges;
    // matchingOneForm = computeMatchingOneForm(*geometries[i], 1, faceGrads);
    // //set up the optimization model 
    // for (Edge e : startEdges) modelBdyEdges.push_back(e.getIndex());
    // for (Edge e : endEdges) modelBdyEdges.push_back(e.getIndex());
    // for (Edge e : meshes[i]->edges()) modelMatchingTerms.push_back(matchingOneForm[e]);

    // model.setBdyEdges(modelBdyEdges);
    // model.setMatchingTerms(modelMatchingTerms);
    // model.setPeriod(period);

    // EdgeData<double> waleOneForm = computeOneForm(*geometries[i], model);
    // //compute stripe patterns from 1-form
    // std::tie(stripeValuesSigma, stripeIndicesSigma) = computeStripeValuesFromOneForm(*geometries[i], waleOneForm, period);
    // std::tie(positions, edges) = generateIsoLines(*geometries[i], stripeValuesSigma, stripeIndicesSigma, period);
    // polyscope::registerCurveNetwork("wale stripe patterns_" + std::to_string(i), positions, edges);
  }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::InputFloat("1-form period", &period);
  
  if (ImGui::Button("Show Stripe Patterns")){
    showStripePatterns();
  }

  if (ImGui::Button("Select Vertex")){
    showBoundaryVertexSelectionUI();
  }
}

int main(int argc, char **argv) {

  numPatches = argc - 1;

  //resize vectors to account for all the input patches 
  meshes.resize(numPatches);
  geometries.resize(numPatches);
  psMeshes.resize(numPatches);
  bdyConditions.resize(numPatches);

  std::ifstream jsonFile(argv[1]);
  jsonFilePath = argv[1];
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  std::string path = data["path"];
  std::vector<std::vector<int>> courseStartBoundary = data["boundaries"]["course"]["start"];
  std::vector<std::vector<int>> courseEndBoundary = data["boundaries"]["course"]["end"];
  std::vector<std::vector<int>> waleStartBoundary = data["boundaries"]["wale"]["start"];
  std::vector<std::vector<int>> waleEndBoundary = data["boundaries"]["wale"]["end"];

  polyscope::init();
  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  //populate the meshes/geometries vectors
  for (int i = 1; i < argc; i++){
    std::tie(meshes[i - 1], geometries[i - 1]) = readManifoldSurfaceMesh(path);
    psMeshes[i - 1] = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(argv[i]),
      geometries[i - 1]->inputVertexPositions, meshes[i - 1]->getFaceVertexList(),
      polyscopePermutations(*meshes[i - 1]));
  }
  
  // Set the callback function
  polyscope::state::userCallback = callBacks;

  polyscope::show();

  return EXIT_SUCCESS;
}