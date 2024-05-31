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

//Let's say for now we're dealing with only 2 patches 

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh1;
std::unique_ptr<VertexPositionGeometry> geometry1;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh1;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh2;
std::unique_ptr<VertexPositionGeometry> geometry2;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh2;

// Some algorithm parameters
float param1 = 42.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void showStripePatterns() {


  std::vector<Vertex> zeroVertices = getBoundaryVertices(*geometry1, 1).first;
  std::vector<Vertex> oneVertices = getBoundaryVertices(*geometry1, 1).second;

  VertexData<double> timeFunction = solveLaplace(*geometry1, zeroVertices, oneVertices);

  psMesh1 -> addVertexScalarQuantity("time function", timeFunction);
  FaceData<Vector3> faceGrads = computeTimeFunctionFaceGrad(*geometry1, timeFunction);
  psMesh1 -> addFaceVectorQuantity("face gradients", faceGrads);
  

  //Compute course stripe patterns
  VertexData<Vector3> courseVertexValuedField = computeVertexValuedField(*geometry1, timeFunction, 0);
  psMesh1 -> addVertexVectorQuantity("course vertex valued field", courseVertexValuedField);
  VertexData<Vector2> lineFieldCourse = vertexDirectionField(*geometry1, courseVertexValuedField);
  double constantFreq = 40.;
  VertexData<double> frequencies(*mesh1, constantFreq);
  CornerData<double> periodicFunc;
  FaceData<int> zeroIndices;
  FaceData<int> branchIndices;
  std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometry1, frequencies, lineFieldCourse);

  // Extract isolines
  std::vector<Vector3> isolineVerts;
  std::vector<std::array<size_t, 2>> isolineEdges;
  std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
    *geometry1, periodicFunc, zeroIndices, branchIndices, lineFieldCourse, false);

  polyscope::registerCurveNetwork("course stripe patterns", isolineVerts, isolineEdges);


  // Compute wale stripe patterns
  VertexData<Vector3> waleVertexValuedField = computeVertexValuedField(*geometry1, timeFunction, PI/2);
  psMesh1 -> addVertexVectorQuantity("wale vertex valued field", waleVertexValuedField);
  VertexData<Vector2> lineFieldWale = vertexDirectionField(*geometry1, waleVertexValuedField);
  std::tie(periodicFunc, zeroIndices, branchIndices) =
      computeStripePattern(*geometry1, frequencies, lineFieldWale);

  // Extract isolines
  std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
    *geometry1, periodicFunc, zeroIndices, branchIndices, lineFieldWale, false);

  polyscope::registerCurveNetwork("wale stripe patterns", isolineVerts, isolineEdges);

}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::SliderFloat("param", &param1, 0., 100.);
  
  if (ImGui::Button("Show Stripe Patterns")) {
    showStripePatterns();
  }
}

int main(int argc, char **argv) {

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = callBacks;

  // Load the first mesh
  std::tie(mesh1, geometry1) = readManifoldSurfaceMesh(argv[1]);

  // Register the mesh with polyscope
  psMesh1 = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(argv[1]),
      geometry1->inputVertexPositions, mesh1->getFaceVertexList(),
      polyscopePermutations(*mesh1));

  // Load the second mesh
  std::tie(mesh2, geometry2) = readManifoldSurfaceMesh(argv[2]);

  // Register the mesh with polyscope
  psMesh2 = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(argv[2]),
      geometry2->inputVertexPositions, mesh2->getFaceVertexList(),
      polyscopePermutations(*mesh2));
 
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}