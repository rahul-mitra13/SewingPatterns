//geometry-central includes 
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/stripe_patterns.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

//external libs
#include "args/args.hxx"
#include "nlohmann/json.hpp"
#include "imgui.h"

//file includes
#include "knitting_utils.h"
#include "stripe_patterns_helpers.h"
#include "iterative_assignment.h"
#include "experiments.h"
#include "KnitGraph.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//vertex mappings from txt file (between panels)
std::vector<std::pair<int, int>> vertexMappingsPairs;
//edge mappings from txt file (between panels)
std::vector<std::pair<int, int>> edgeMappingsPairs;
//build an index map from vertices in the original mesh to vertices in the glued mesh 
std::map<int, int> vertexMap;
//build an index map from edges in the orignal mesh to edges in the glued
std::map<int, int> edgeMap;
//one-ring map for vertices in the glued mesh for performance 
std::map<int, std::vector<Halfedge>> gluedOneRingMap;
//build a glued edge length geoemtry to make our life a little easier for some procedures
//EdgeLengthGeometry *gluedELG;
std::unique_ptr<EdgeLengthGeometry> gluedELG;

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
//knoppel frequency in the stripe patterns
float knoppelFrequency = 0.0;

//set the permutation of edges and orientations for 1-form viz
std::vector<size_t> perm;
std::vector<bool> orientations;

//eigen matrix vertex positions
Eigen::MatrixXd V;
//eigen matrix face lists
Eigen::MatrixXi F;
// Compute the global gradient operator: #F*3 by #V
Eigen::SparseMatrix<double> grad;
//gradient operator in the row major format
Eigen::SparseMatrix<double, Eigen::RowMajor> G;

//here we will do as much processing as possible directly on the glued together mesh 
void showStripePatterns(){
  
  //time function on the glued mesh 
  VertexData<double> timeFunctionGlued = computeTimeFunction(*gluedELG, globalBdyConditions);
  //time function on the global mesh 
  VertexData<double> timeFunctionGlobal = convertGluedToGlobalVertexFunction(*globalGeometry, *gluedELG, timeFunctionGlued, vertexMap);
  globalPSMesh -> addVertexScalarQuantity("time function", timeFunctionGlobal);

  // Compute course stripes using Knöppel's method
  VertexData<Vector3> vertexVectorField = computeVertexValuedField(*globalGeometry, timeFunctionGlobal, 0.0);
  VertexData<Vector2> lineField = vertexDirectionField(*globalGeometry, vertexVectorField);
  VertexData<double> freq(globalGeometry->mesh, 1./(period));
  CornerData<double> stripeValues(globalGeometry->mesh);
  FaceData<int> stripeSingularities(globalGeometry->mesh);
  FaceData<int> fieldSingularities(globalGeometry->mesh);
  std::tie(stripeValues, stripeSingularities, fieldSingularities) = computeStripePattern(*globalGeometry, freq, lineField); // this is a GC call
  // Do some visualization
  globalPSMesh->addVertexVectorQuantity("vertexVectorField", vertexVectorField);
  globalPSMesh->addFaceScalarQuantity("knoppel face singularities", stripeSingularities);
  globalPSMesh->addFaceScalarQuantity("knoppel field singularities", fieldSingularities);
  globalPSMesh->addCornerScalarQuantity("stripeValues", prepareCornerData(stripeValues));
  std::vector<Vector3> knoppelPos; 
  std::vector<std::array<size_t, 2>> knoppelEdges; 
  std::tie(knoppelPos, knoppelEdges) = extractPolylinesFromStripePattern(*globalGeometry, stripeValues, stripeSingularities,
                                          fieldSingularities, lineField, false);
  auto knoppelStripes = polyscope::registerCurveNetwork("knoppel wale stripes stripes", knoppelPos, knoppelEdges);
  knoppelStripes->setRadius(0.001);
  knoppelStripes->setEnabled(false);

  //gradient on the glued/global mesh
  //note that faces have a 1-to-1 mapping from global to glued setting
  FaceData<Vector3> timeFunctionGradientGlobal = computeTimeFunctionFaceGrad(*globalGeometry, timeFunctionGlobal);
  //normalize the gradient of the timeFunction 
  FaceData<Vector3> timeFunctionGradientGlobalNormalized(*globalMesh);
  for (Face f : globalMesh -> faces()){
    timeFunctionGradientGlobalNormalized[f] = timeFunctionGradientGlobal[f].normalize();
  }
  globalPSMesh -> addFaceVectorQuantity("normalized time function gradient", timeFunctionGradientGlobalNormalized);


  // // Call Matteo's revealCurl
  // Model model;
  // model.setPeriod(period);
  // model.setBdyEdges(globalBdyConditions.courseBdyEdges);
  // revealCurl(*globalGeometry, *gluedELG, model, vertexMap, G);
  // P("revealCurl done");
  // polyscope::show();

  //-------iteratively find course stripes and course singularities----------//
  CornerData<double> courseStripeValues(globalGeometry -> mesh);
  EdgeData<double> courseSingularEdgesGlobal(globalGeometry -> mesh);
  G = grad;
  FaceData<Vector3> courseOneFormGrad(globalGeometry -> mesh);
  std::tie(courseStripeValues, courseSingularEdgesGlobal) = implCourseHarmonic1Form(*globalGeometry, *gluedELG, timeFunctionGlobal,
                                                                    timeFunctionGradientGlobalNormalized, vertexMap, edgeMap, *globalPSMesh,
                                                                    globalBdyConditions, period, V, F, G, courseOneFormGrad, gluedOneRingMap);

  std::tie(courseStripeValues, courseSingularEdgesGlobal);

  globalPSMesh -> addEdgeScalarQuantity("course singular edges", courseSingularEdgesGlobal);


  //WALE STRIPES
  //find Knöppel singularities in the WALE DIRECTION 
  //just run Knoppel's algorithm on these models 
  //and then run our 1-form optimization with the singularities
  CornerData<double> waleStripeValues;
  EdgeData<double> waleSingularEdgesGlobal;
  FaceData<int> waleSingularFaces(globalGeometry -> mesh, 0);//there are no face singularities
  std::tie(waleStripeValues, waleSingularEdgesGlobal) = computeWaleStripeInfo(*globalGeometry, *gluedELG, 
                                                                    edgeMappingsPairs, edgeMap, vertexMap, timeFunctionGlobal, 
                                                                    courseOneFormGrad, G, period, knoppelFrequency, globalBdyConditions, 
                                                                    courseSingularEdgesGlobal, *globalPSMesh);
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, waleStripeValues, waleSingularFaces, period);
  auto waleStripes = polyscope::registerCurveNetwork("wale stripes", positionsWale, edgesWale);
  globalPSMesh -> addEdgeScalarQuantity("wale singularities", waleSingularEdgesGlobal);
  waleStripes -> setRadius(0.001);
  waleStripes -> setEnabled(false);


  //show vertex curl and edge curl of normalized time function gradient
  // VertexData<double> vertexCurl = computeVertexCurl(*globalGeometry, *gluedELG, 
  //                                   timeFunctionGradientGlobalNormalized, gluedOneRingMap);
  // EdgeData<double> edgeCurl = computeVertexAveragedEdgeCurl(*globalGeometry, vertexCurl);
  // globalPSMesh->addVertexScalarQuantity("normalized time function gradient vertex curl", vertexCurl);
  // globalPSMesh->addEdgeScalarQuantity("normalized time function gradient edge curl", edgeCurl);

  //generate the knit graph
  KnitGraph graph = KnitGraph(*globalGeometry, *gluedELG, *globalPSMesh, period, 
                      courseStripeValues, courseSingularEdgesGlobal, waleStripeValues, waleSingularEdgesGlobal,
                      edgeMap);
  graph.buildGraph();

}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::InputFloat("1-form period", &period);
  //there needs to be a better way to constrain wale edges
  ImGui::InputFloat("Threshold", &threshold);
  //frequency for knoppel stripes
  ImGui::InputFloat("Knoppel frequency", &knoppelFrequency);

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
  fixDelaunay(*globalMesh, *globalGeometry); // we make the mesh approximately Delaunay
  std::tie(V, F) = getVertexPositionsandFaceLists(*globalGeometry);
  igl::grad(V,F,grad);

  // Parse stripe period, if available
  if (argc > 2) {
    sscanf(argv[2], "%f", &period);
    knoppelFrequency = 1.0 / period; // to match course and wale periods
  }

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

  // Internally, Polyscope numbers the edges by looping over faces.
  // Since our numbering is different than that after fixDelaunay, we need to specify the new numbering by providing a permutation.
  // Note that this is only useful for EdgeData visualization; apart from that everything works fine.
  std::set<Edge> visited;
  for (Face f : globalMesh->faces()) {
    for (Edge e : f.adjacentEdges()) {
      if (visited.count(e) == 0) {
        perm.push_back(e.getIndex());
        visited.insert(e);
      }
    }
  }
  vertexMappingsPairs = buildPairOfStitchedVerticesFromFile(data["vertex_mappings"]);
  edgeMappingsPairs = buildPairOfStitchedEdges(*globalGeometry, vertexMappingsPairs);
  //set up the orientations for the 1-form viz while we're here
  for (Edge e : globalMesh->edges()){
    if (e.halfedge().tailVertex().getIndex() < e.halfedge().tipVertex().getIndex()){
      orientations.push_back(true);
    }
    else{
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

  // Halfedge permutation (global -> glued)
  // Note: this works as long as the global and glued faces have the same order.
  // If not, you should iterate over global mesh halfedges and find the corresponding half-edge in the glued mesh.
  // I was just too lazy to implement that.
  std::vector<size_t> halfedgePerm;
  for (Face f : gluedELG->mesh.faces())
    for (Halfedge he : f.adjacentHalfedges())
      halfedgePerm.push_back(he.getIndex());  
  globalPSMesh->setHalfedgePermutation(halfedgePerm);

  //process boundary conditions in the glued mesh setting 
  globalBdyConditions = parseJson(*gluedELG, data, vertexMap, edgeMap);

  //view the global boundary conditions
  // VertexData<double> oneVertices(*globalMesh, 0);
  // VertexData<double> zeroVertices(*globalMesh, 0);
  // for (int v : globalBdyConditions.courseStartBoundaryVertices){
  //   zeroVertices[v] = 1.0;
  // }
  // for (int v : globalBdyConditions.courseEndBoundaryVertices){
  //   oneVertices[v] = 1.0;
  // }
  // globalPSMesh->addVertexScalarQuantity("ones", oneVertices);
  // globalPSMesh->addVertexScalarQuantity("zeros", zeroVertices);

  //render the stitched vertices
  renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  // Disable the ground plane
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;

  if (argc > 2) // period was provided
    showStripePatterns();

  polyscope::show();

  return EXIT_SUCCESS; 
}