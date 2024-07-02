//this is to appropriately parse the GarmentCode specification
#include <Python.h>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/vector_heat_method.h"



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
std::map<int, int> vertexMappings;
//also build a map that stores panel name to geometry (can infer the mesh from the geometry)
//std::map<std::string, std::unique_ptr<VertexPositionGeometry>> panelMappings;

std::unique_ptr<ManifoldSurfaceMesh> globalMesh;
std::unique_ptr<VertexPositionGeometry> globalGeometry;
polyscope::SurfaceMesh * globalPSMesh;

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
  VertexData<Vector2> lineFieldCourse = computeSmoothestVertexDirectionField(*globalGeometry, 1);
  globalGeometry->requireVertexTangentBasis();
    VertexData<Vector3> vBasisX(*globalMesh);
    VertexData<Vector3> vBasisY(*globalMesh);
    for(Vertex v : globalMesh->vertices()) {
      vBasisX[v] = globalGeometry->vertexTangentBasis[v][0];
      vBasisY[v] = globalGeometry->vertexTangentBasis[v][1];
    }
    globalPSMesh ->
    addVertexTangentVectorQuantity("vertex tangent vectors", lineFieldCourse, vBasisX, vBasisY);
    VertexData<double> frequencies(*globalMesh, constantFreq);
    CornerData<double> periodicFunc;
    FaceData<int> zeroIndices;
    FaceData<int> branchIndices;
    std::tie(periodicFunc, zeroIndices, branchIndices) =
        computeStripePattern(*globalGeometry, frequencies, lineFieldCourse);
    //Extract isolines
    std::vector<Vector3> isolineVerts;
    std::vector<std::array<size_t, 2>> isolineEdges;
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *globalGeometry, periodicFunc, zeroIndices, branchIndices, lineFieldCourse, false);
    polyscope::registerCurveNetwork("course stripe patterns", isolineVerts, isolineEdges);
  
/**
  for (int i = 0; i < meshes.size(); i++){
    std::vector<Vertex> zeroVertices, oneVertices;
    zeroVertices = bdyConditions[i].courseStartBoundaryVertices;
    oneVertices = bdyConditions[i].courseEndBoundaryVertices;
    
    VertexData<double> timeFunction = solveLaplace(*geometries[i], zeroVertices, oneVertices);
    psMeshes[i] -> addVertexScalarQuantity("time function_" + std::to_string(i), timeFunction);
    FaceData<Vector3> faceGrads = computeTimeFunctionFaceGrad(*geometries[i], timeFunction);
    psMeshes[i] -> addFaceVectorQuantity("face gradients_" + std::to_string(i), faceGrads);
    std::vector<Edge> startEdges, endEdges;
    //set up the optimization model 
    Model model;
    //the 1-form we're trying to match
    EdgeData<double>  matchingOneForm;
    //put all the boundary edges into a single vector 
    std::vector<int> modelBdyEdges;
    //put the edge terms we're trying to match in a vector
    std::vector<double> modelMatchingTerms;
    CornerData<double> stripeValuesSigma;
    FaceData<int> stripeIndicesSigma;  
    std::vector<Vector3> positions;
    std::vector<std::array<int, 2>> edges;

    //Compute knoppel course stripe patterns
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
    polyscope::registerCurveNetwork("wale stripe patterns_" + std::to_string(i), isolineVerts, isolineEdges);
    
    
    //generate course stripe patterns using 1-form approach 
    startEdges = bdyConditions[i].courseStartBoundaryEdges;
    endEdges = bdyConditions[i].courseEndBoundaryEdges;
    matchingOneForm = computeMatchingOneForm(*geometries[i], 0, faceGrads);
    for (Edge e : startEdges) modelBdyEdges.push_back(e.getIndex());
    for (Edge e : endEdges) modelBdyEdges.push_back(e.getIndex());
    for (Edge e : meshes[i]->edges()) modelMatchingTerms.push_back(matchingOneForm[e]);
    model.setBdyEdges(modelBdyEdges);
    model.setMatchingTerms(modelMatchingTerms);
    model.setPeriod(period);
    EdgeData<double> courseOneForm = computeOneForm(*geometries[i], model);
    //compute course stripe patterns from 1-form
    std::tie(stripeValuesSigma, stripeIndicesSigma) = computeStripeValuesFromOneForm(*geometries[i], courseOneForm, period);
    std::tie(positions, edges) = generateIsoLines(*geometries[i], stripeValuesSigma, stripeIndicesSigma, period);
    polyscope::registerCurveNetwork("course stripe patterns_" + std::to_string(i), positions, edges);
    
    //clear out course info before computing wale stripes
    modelBdyEdges.clear();
    modelMatchingTerms.clear();

    //generate wale stripe patterns using 1-form approach 
    startEdges = bdyConditions[i].waleStartBoundaryEdges;
    endEdges = bdyConditions[i].waleEndBoundaryEdges;
    matchingOneForm = computeMatchingOneForm(*geometries[i], 1, faceGrads);
    //set up the optimization model 
    for (Edge e : startEdges) modelBdyEdges.push_back(e.getIndex());
    for (Edge e : endEdges) modelBdyEdges.push_back(e.getIndex());
    for (Edge e : meshes[i]->edges()) modelMatchingTerms.push_back(matchingOneForm[e]);

    model.setBdyEdges(modelBdyEdges);
    model.setMatchingTerms(modelMatchingTerms);
    model.setPeriod(period);

    EdgeData<double> waleOneForm = computeOneForm(*geometries[i], model);
    //compute stripe patterns from 1-form
    std::tie(stripeValuesSigma, stripeIndicesSigma) = computeStripeValuesFromOneForm(*geometries[i], waleOneForm, period);
    std::tie(positions, edges) = generateIsoLines(*geometries[i], stripeValuesSigma, stripeIndicesSigma, period);
    polyscope::registerCurveNetwork("wale stripe patterns_" + std::to_string(i), positions, edges);    
  }
*/
  
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
    std::cout << "Is mesh manifold? " << globalMesh -> isManifold() << std::endl;
    std::cout << "Is mesh oriented? " << globalMesh -> isOriented() << std::endl;
    std::cout << "Is mesh triangular? " << globalMesh -> isTriangular() << std::endl;
    std::cout << "Number of connected components of mesh " << globalMesh -> nConnectedComponents() << std::endl;
    globalPSMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(argv[1]), globalGeometry->inputVertexPositions, globalMesh -> getFaceVertexList());
    vertexMappings = buildGlobalVertexMappingFromFile(argv[2]);
    Eigen::SparseMatrix<double> L;
    buildGlobalCotanLaplacian(*globalGeometry, vertexMappings, L);
  }

  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;
  polyscope::show();

  return EXIT_SUCCESS;
}