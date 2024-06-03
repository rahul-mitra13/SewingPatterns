#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/stripe_patterns.h"

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
#include "helpers.h"//for testing now, gonna remove eventually

using namespace geometrycentral;
using namespace geometrycentral::surface;

int numPatches;//number of patches

//make a vector to handle all the patches 
//meshes + geometries 
std::vector<std::unique_ptr<ManifoldSurfaceMesh>> meshes;
std::vector<std::unique_ptr<VertexPositionGeometry>> geometries;
std::vector<polyscope::SurfaceMesh *> psMeshes;


// Striping frequency
float constantFreq = 40.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void showStripePatterns() {


  //do this for every patch
  for (int i = 0; i < numPatches; i++){  
    std::vector<Vertex> zeroVertices = getBoundaryVertices(*geometries[i], 1).first;
    std::vector<Vertex> oneVertices = getBoundaryVertices(*geometries[i], 1).second;

    VertexData<double> timeFunction = solveLaplace(*geometries[i], zeroVertices, oneVertices);

    psMeshes[i] -> addVertexScalarQuantity("time function_" + std::to_string(i), timeFunction);
    FaceData<Vector3> faceGrads = computeTimeFunctionFaceGrad(*geometries[i], timeFunction);
    psMeshes[i] -> addFaceVectorQuantity("face gradients_" + std::to_string(i), faceGrads);
  

    //Compute course stripe patterns
    VertexData<Vector3> courseVertexValuedField = computeVertexValuedField(*geometries[i], timeFunction, 0);
    psMeshes[i] -> addVertexVectorQuantity("course vertex valued field_" + std::to_string(i), courseVertexValuedField);
    VertexData<Vector2> lineFieldCourse = vertexDirectionField(*geometries[i], courseVertexValuedField);
    VertexData<double> frequencies(*meshes[i], constantFreq);
    CornerData<double> periodicFunc;
    FaceData<int> zeroIndices;
    FaceData<int> branchIndices;
    std::tie(periodicFunc, zeroIndices, branchIndices) =
        computeStripePattern(*geometries[i], frequencies, lineFieldCourse);

    // Extract isolines
    std::vector<Vector3> isolineVerts;
    std::vector<std::array<size_t, 2>> isolineEdges;
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometries[i], periodicFunc, zeroIndices, branchIndices, lineFieldCourse, false);

    polyscope::registerCurveNetwork("course stripe patterns_" + std::to_string(i), isolineVerts, isolineEdges);


    // // Compute wale stripe patterns
    VertexData<Vector3> waleVertexValuedField = computeVertexValuedField(*geometries[i], timeFunction, PI/2);
    psMeshes[i] -> addVertexVectorQuantity("wale vertex valued field_" + std::to_string(i), waleVertexValuedField);
    VertexData<Vector2> lineFieldWale = vertexDirectionField(*geometries[i], waleVertexValuedField);
    std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometries[i], frequencies, lineFieldWale);

    // Extract isolines
    std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
      *geometries[i], periodicFunc, zeroIndices, branchIndices, lineFieldWale, false);

    polyscope::registerCurveNetwork("wale stripe patterns_" + std::to_string(i), isolineVerts, isolineEdges);
  }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::SliderFloat("Striping Frequency", &constantFreq, 0., 100.);
  
  if (ImGui::Button("Show Stripe Patterns")) {
    showStripePatterns();
  }
}

int main(int argc, char **argv) {

  numPatches = argc - 1;

  std::cout << "Testing windows build" << std::endl;

  //resize vectors to account for all the input patches 
  meshes.resize(numPatches);
  geometries.resize(numPatches);
  psMeshes.resize(numPatches);

  polyscope::init();

  //populate the meshes/geometries vectors
  for (int i = 1; i < argc; i++){
    std::tie(meshes[i - 1], geometries[i - 1]) = readManifoldSurfaceMesh(argv[i]);
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