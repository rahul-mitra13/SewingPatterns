//geometry-central includes 
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/stripe_patterns.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"

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
//#include "experiments.h"
#include "KnitGraph.h"

#include "homology_generators.h"

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

//pre-processed meshes and geometries
std::unique_ptr<ManifoldSurfaceMesh> globalMeshPre;
std::unique_ptr<VertexPositionGeometry> globalGeometryPre;

std::unique_ptr<ManifoldSurfaceMesh> globalMesh;
std::unique_ptr<VertexPositionGeometry> globalGeometry;
polyscope::SurfaceMesh *globalPSMesh;
//global boundary conditions
globalBoundaryConditions globalBdyConditions;

//1-form optimization course Period
float coursePeriod = 1;
float walePeriod;//the ratio of width to height should be 1:1.6 (Kui)
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

//the knit graph over the model 
KnitGraph graph;

//here we will do as much processing as possible directly on the glued together mesh 
void showStripePatterns(){

  H(coursePeriod);
  
  //set the wale period 
  walePeriod = ((1.0/1.6) * coursePeriod); // rendering setting
  // walePeriod = 1.0 * coursePeriod; // knitting setting
  //time function on the glued mesh 
  VertexData<double> timeFunctionGlued = computeTimeFunction(*gluedELG, globalBdyConditions);
  //time function on the global mesh 
  VertexData<double> timeFunctionGlobal = convertGluedToGlobalVertexFunction(*globalGeometry, *gluedELG, timeFunctionGlued, vertexMap);
  globalPSMesh -> addVertexScalarQuantity("time function", timeFunctionGlobal);

  // Compute a smooth 1-direction field
  FaceData<Vector2> smoothDirectionField = computeSmoothestBoundaryAlignedFaceDirectionField(*globalGeometry, 1);
  FaceData<Vector3> basisX(globalGeometry->mesh);
  FaceData<Vector3> basisY(globalGeometry->mesh);
  globalGeometry->requireFaceTangentBasis();
  for (Face f : globalGeometry->mesh.faces()) {
    basisX[f] = globalGeometry->faceTangentBasis[f][0];
    basisY[f] = globalGeometry->faceTangentBasis[f][1];
  }
  globalPSMesh->addFaceTangentVectorQuantity("smooth direction field", smoothDirectionField, basisX, basisY);

  // Compute course stripes using Knöppel's method
  VertexData<Vector3> vertexVectorField = computeVertexValuedField(*globalGeometry, timeFunctionGlobal, 0.0);
  VertexData<Vector2> usedRoot(globalGeometry->mesh);
  VertexData<Vector2> lineField = vertexDirectionField(*globalGeometry, vertexVectorField, usedRoot);
  VertexData<double> freq(globalGeometry->mesh, 1./(coursePeriod));
  CornerData<double> stripeValues(globalGeometry->mesh);
  FaceData<int> stripeSingularities(globalGeometry->mesh);
  FaceData<int> fieldSingularities(globalGeometry->mesh);
  std::tie(stripeValues, stripeSingularities, fieldSingularities) = computeStripePattern(*globalGeometry, freq, lineField); // this is a GC call
  // Do some visualization
  globalPSMesh->addVertexVectorQuantity("vertexVectorField", vertexVectorField);
  globalPSMesh->addFaceScalarQuantity("knoppel course face singularities", stripeSingularities);
  globalPSMesh->addFaceScalarQuantity("knoppel course field singularities", fieldSingularities);
  globalPSMesh->addCornerScalarQuantity("knoppel course stripe values", prepareCornerData(stripeValues));
  std::vector<Vector3> knoppelPos; 
  std::vector<std::array<size_t, 2>> knoppelEdges; 
  std::tie(knoppelPos, knoppelEdges) = extractPolylinesFromStripePattern(*globalGeometry, stripeValues, stripeSingularities,
                                          fieldSingularities, lineField, false);
  auto knoppelStripes = polyscope::registerCurveNetwork("knoppel course stripes", knoppelPos, knoppelEdges);
  knoppelStripes->setRadius(0.001);
  knoppelStripes->setEnabled(false);

  //loops from the saddle vertex for meshes with genus 
  std::vector<std::vector<double>> allSaddleLoops;
  //all the homology generators 
  std::vector<std::vector<double>> homologyGenerators;
  if (globalGeometry->mesh.nConnectedComponents() == 1){//do the homology generator stuff for just global 3D meshes for now 
      std::vector<Vertex> saddleVertices = getSaddleVertices(*globalGeometry, timeFunctionGlobal);
      allSaddleLoops = findAllSaddleLoops(*globalGeometry, saddleVertices, timeFunctionGlobal);
      homologyGenerators = buildHomologyGeneratorsVector(*globalGeometry, *globalMesh);
  }


  //gradient on the glued/global mesh
  //note that faces have a 1-to-1 mapping from global to glued setting
  FaceData<Vector3> timeFunctionGradientGlobal = computeTimeFunctionFaceGrad(*globalGeometry, timeFunctionGlobal);
  //normalize the gradient of the timeFunction 
  FaceData<Vector3> timeFunctionGradientGlobalNormalized(*globalMesh);
  for (Face f : globalMesh -> faces()){
    timeFunctionGradientGlobalNormalized[f] = timeFunctionGradientGlobal[f].normalize();
  }
  globalPSMesh -> addFaceVectorQuantity("normalized time function gradient", timeFunctionGradientGlobalNormalized);

  // Compute intrinsic time function gradient
  gluedELG->requireEdgeLengths();
  FaceData<Vector2> timeFunctionGradientIntrinsic = computeTimeFunctionFaceGradIntrinsic(*gluedELG, timeFunctionGlued);
  for (Face f : gluedELG->mesh.faces()) timeFunctionGradientIntrinsic[f] = timeFunctionGradientIntrinsic[f].normalize();
  VertexData<double> timeFunctionAligned = computeTimeFunctionAligned(*gluedELG, globalBdyConditions, timeFunctionGradientIntrinsic);
  VertexData<double> timeFunctionAlignedGlobal = convertGluedToGlobalVertexFunction(*globalGeometry, *gluedELG, timeFunctionAligned, vertexMap);
  globalPSMesh->addVertexScalarQuantity("time function aligned", timeFunctionAlignedGlobal)->setEnabled(false);

  G = grad;

  //-------iteratively find COURSE stripes and course singularities----------//

  CornerData<double> courseStripeValues(globalGeometry -> mesh);
  EdgeData<double> courseSingularEdgesGlobal(globalGeometry -> mesh);
  FaceData<Vector3> courseOneFormGrad(globalGeometry -> mesh);

  std::string courseDataFile;
  courseDataFile = "../data/fertility/courseData.ply";

  if (courseDataFile.size() > 0) { // if a file is specified, load it
    RichSurfaceMeshData courseData(globalGeometry->mesh, courseDataFile);
    courseStripeValues = courseData.getCornerProperty<double>("courseStripeValues");
    courseSingularEdgesGlobal = courseData.getEdgeProperty<double>("courseSingularEdges");
  } else {
    std::tie(courseStripeValues, courseSingularEdgesGlobal) = implCourseHarmonic1Form(*globalGeometry, *gluedELG, timeFunctionAlignedGlobal,
                                                                      timeFunctionGradientGlobalNormalized, vertexMap, edgeMap, *globalPSMesh,
                                                                      globalBdyConditions, coursePeriod, V, F, G, courseOneFormGrad, gluedOneRingMap, 
                                                                      allSaddleLoops, homologyGenerators);
    RichSurfaceMeshData courseData(globalGeometry->mesh);
    courseData.addCornerProperty("courseStripeValues", courseStripeValues);
    courseData.addEdgeProperty("courseSingularEdges", courseSingularEdgesGlobal);
    courseData.write("courseData.ply");
  }
  
  // Store course singular edges for rendering later (Rahul's code)
  RichSurfaceMeshData richData(globalGeometry->mesh);
  richData.addMeshConnectivity();
  richData.addGeometry(*globalGeometry);
  richData.addEdgeProperty("courseSingularEdge", courseSingularEdgesGlobal);

  globalPSMesh -> addEdgeScalarQuantity("course singular edges", courseSingularEdgesGlobal);

  // Draw course singular edges
  std::vector<Vector3> singularEdgePointsPos;
  std::vector<Vector3> singularEdgePointsNeg;
  std::vector<std::array<int,2>> singularEdgesPos;
  std::vector<std::array<int,2>> singularEdgesNeg;
  for (Edge e : globalMesh->edges()) if (courseSingularEdgesGlobal[e] != 0) {
      if (courseSingularEdgesGlobal[e] > 0) {
        singularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
        singularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
        singularEdgesPos.push_back({(int)singularEdgePointsPos.size()-2, (int)singularEdgePointsPos.size()-1});
      } else {
        singularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
        singularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
        singularEdgesNeg.push_back({(int)singularEdgePointsNeg.size()-2, (int)singularEdgePointsNeg.size()-1});
      }
  }
  polyscope::registerCurveNetwork("course singular edges (+1)", singularEdgePointsPos, singularEdgesPos)->setRadius(0.003)->setColor({1,0,0})->setEnabled(false);
  polyscope::registerCurveNetwork("course singular edges (-1)", singularEdgePointsNeg, singularEdgesNeg)->setRadius(0.003)->setColor({0,0,1})->setEnabled(false);


  FaceData<int> stripeIndicesSigmaCourse(*globalMesh, 0);
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  // std::tie(positionsCourse, edgesCourse) = generateIsoLines(*globalGeometry, courseStripeValues, stripeIndicesSigmaCourse, coursePeriod);

  std::unordered_map<size_t, std::vector<PolyLinePoint>> components; // don't really need this
  std::tie(positionsCourse, edgesCourse) = findStripeConnectedComponents(*globalGeometry, *gluedELG, courseStripeValues, stripeIndicesSigmaCourse, coursePeriod, edgeMap, components);
  auto courseStripes = polyscope::registerCurveNetwork("course stripes ", positionsCourse, edgesCourse)->setRadius(0.001)->setEnabled(false);
  globalPSMesh->addCornerScalarQuantity("courseStripeValues", prepareCornerData(courseStripeValues));

  // Plot stripe values with offset (to debug knit graph)
  CornerData<double> courseStripeValuesWithOffset(courseStripeValues);
  for (Corner co : globalMesh->corners())
      courseStripeValuesWithOffset[co] -= coursePeriod/4;
  std::tie(positionsCourse, edgesCourse) = findStripeConnectedComponents(*globalGeometry, *gluedELG, courseStripeValuesWithOffset, stripeIndicesSigmaCourse, coursePeriod, edgeMap, components);
  //std::cout << "Number of components after " << std::to_string(numSingularities) << " singularity insertions is " << components.size() << std::endl;
  polyscope::registerCurveNetwork("course stripes with offset", positionsCourse, edgesCourse)->setRadius(0.002)->setColor({50.0/255, 205.0/255, 50.0/255})->setEnabled(false);

  // polyscope::show();

  //-------iteratively find WALE stripes and course singularities----------//

  //find Knöppel singularities in the WALE DIRECTION 
  //just run Knoppel's algorithm on these models 
  //and then run our 1-form optimization with the singularities
  CornerData<double> waleStripeValues;
  EdgeData<double> waleSingularEdgesGlobal;
  FaceData<int> waleSingularFaces(globalGeometry -> mesh, 0);//there are no face singularities

  std::string waleDataFile;
  // waleDataFile = "bobWaleData10pairs.ply";
  // waleDataFile = "bobPeriod2.0/waleData.ply";
  waleDataFile = "../data/fertility/waleData.ply";

  if (waleDataFile.size() > 0) { // if a file is specified, load it
    RichSurfaceMeshData waleData(globalGeometry->mesh, waleDataFile);
    waleStripeValues = waleData.getCornerProperty<double>("waleStripeValues");
    waleSingularEdgesGlobal = waleData.getEdgeProperty<double>("waleSingularEdges");
  } else {

    std::tie(waleStripeValues, waleSingularEdgesGlobal) = computeWaleStripeInfo(*globalGeometry, *gluedELG, 
                                                                      edgeMappingsPairs, edgeMap, vertexMap, timeFunctionGlobal, timeFunctionGlued,
                                                                      courseOneFormGrad, G, walePeriod, knoppelFrequency, globalBdyConditions, 
                                                                      courseSingularEdgesGlobal, gluedOneRingMap,*globalPSMesh, allSaddleLoops,
                                                                      homologyGenerators);

    RichSurfaceMeshData waleData(globalGeometry->mesh);
    waleData.addEdgeProperty("waleSingularEdges", waleSingularEdgesGlobal);
    waleData.addCornerProperty("waleStripeValues", waleStripeValues);
    waleData.write("waleData.ply");
  }

  // Store wale singular edges for rendering later (Rahul's stuff)
  richData.addEdgeProperty("waleSingularEdges", waleSingularEdgesGlobal);
  // richData.write("singularEdges.ply");

  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, waleStripeValues, waleSingularFaces, walePeriod);
  auto waleStripes = polyscope::registerCurveNetwork("wale stripes", positionsWale, edgesWale);
  globalPSMesh -> addEdgeScalarQuantity("wale singularities", waleSingularEdgesGlobal);
  waleStripes -> setRadius(0.001);
  waleStripes -> setEnabled(false);

  // Draw wale singular edges as a curve network
  std::vector<Vector3> waleSingularEdgePointsPos;
  std::vector<Vector3> waleSingularEdgePointsNeg;
  std::vector<std::array<int,2>> waleSingularEdgesPos;
  std::vector<std::array<int,2>> waleSingularEdgesNeg;
  for (Edge e : globalMesh->edges()) if (waleSingularEdgesGlobal[e] != 0) {
      if (waleSingularEdgesGlobal[e] > 0) {
        waleSingularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
        waleSingularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
        waleSingularEdgesPos.push_back({(int)waleSingularEdgePointsPos.size()-2, (int)waleSingularEdgePointsPos.size()-1});
      } else {
        waleSingularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
        waleSingularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
        waleSingularEdgesNeg.push_back({(int)waleSingularEdgePointsNeg.size()-2, (int)waleSingularEdgePointsNeg.size()-1});
      }
  }
  polyscope::registerCurveNetwork("wale singular edges (+1)", waleSingularEdgePointsPos, waleSingularEdgesPos)->setRadius(0.003)->setColor({1,0,0})->setEnabled(false);
  polyscope::registerCurveNetwork("wale singular edges (-1)", waleSingularEdgePointsNeg, waleSingularEdgesNeg)->setRadius(0.003)->setColor({0,0,1})->setEnabled(false);


  globalPSMesh->addCornerScalarQuantity("wale stripe values", prepareCornerData(waleStripeValues));

  // Plot wale stripe values with offset (to debug knit graph)
  CornerData<double> waleStripeValuesWithOffset(waleStripeValues);
  for (Corner co : globalGeometry->mesh.corners())
      waleStripeValuesWithOffset[co] -= walePeriod/4;
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, waleStripeValuesWithOffset, waleSingularFaces, walePeriod);
  polyscope::registerCurveNetwork("wale stripes with offset", positionsWale, edgesWale)->setRadius(0.002)->setColor({1,140./255,0})->setEnabled(false);
  
  std::cout << "Done with wale stripes" << std::endl;
  // polyscope::show();

  // generate the knit graph
  graph = KnitGraph(*globalGeometry, *gluedELG, *globalPSMesh, coursePeriod, walePeriod,
                      courseStripeValues, courseSingularEdgesGlobal, waleStripeValues, waleSingularEdgesGlobal,
                      edgeMap);
  graph.buildGraph();

  graph.writeKnitGraphToTxtFile("model.obj");


}

//repair a knit graph vertex that's missing connections
//ugh 
void manualKnitGraphRepair(){

  std::vector<int> newInfo;
  newInfo = repairKnitGraphVertex();

  int id = newInfo[0];
  int row_in = newInfo[1];
  int row_out = newInfo[2];
  int col_in_1 = newInfo[3];
  int col_in_2 = newInfo[4];
  int col_out_1 = newInfo[5];
  int col_out_2 = newInfo[6];

  //update the changes 
  graph.getRealVertices()[id].row_in = row_in;
  graph.getRealVertices()[row_in].row_out = id;

  graph.getRealVertices()[id].row_out = row_out;
  graph.getRealVertices()[row_out].row_in = id;

  graph.getRealVertices()[id].col_in[0] = col_in_1;
  graph.getRealVertices()[col_in_1].col_out[0] = id;

  graph.getRealVertices()[id].col_in[1] = col_in_2;
  graph.getRealVertices()[col_in_2].col_out[1] = id;

  graph.getRealVertices()[id].col_out[0] = col_out_1;
  graph.getRealVertices()[col_out_1].col_in[0] = id;

  graph.getRealVertices()[id].col_out[1] = col_out_2;
  graph.getRealVertices()[col_out_2].col_in[1] = id;

  graph.renderGraph();
  graph.sanityCheck();

  //write the graph to a txt file 
  graph.writeKnitGraphToTxtFile("model.obj");

}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void callBacks() {

  ImGui::InputFloat("Course 1-form period", &coursePeriod);
  //there needs to be a better way to constrain wale edges
  //ImGui::InputFloat("Threshold", &threshold);
  //frequency for knoppel stripes
  ImGui::InputFloat("Knoppel frequency", &knoppelFrequency);

  if (ImGui::Button("Show Stripe Patterns")){
    showStripePatterns();
  }

  if (ImGui::Button("Manaul Repair")){
    manualKnitGraphRepair();
  }

  if (ImGui::Button("Save view as JSON string")){
    std::cout << polyscope::view::getViewAsJson() << std::endl;
  }
}

int main(int argc, char **argv) {
  polyscope::init();
  std::ifstream jsonFile(argv[1]);
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  //run sanity checks
  std::tie(globalMeshPre, globalGeometryPre) = readManifoldSurfaceMesh(data["model_path"]);
  
  //make the mesh Delaunay 
  fixDelaunay(*globalMeshPre, *globalGeometryPre); // we make the mesh approximately Delaunay

  // // Remesh
  // remesh(*globalMeshPre, *globalGeometryPre);

  std::tie(V, F) = getVertexPositionsandFaceLists(*globalGeometryPre);
  globalMesh = std::make_unique<ManifoldSurfaceMesh>(F);
  globalGeometry = std::make_unique<VertexPositionGeometry>(*globalMesh, V);

  for (auto bl : globalMesh -> boundaryLoops()){
    Face f = bl.asFace();
    std::cout << "boundary vertex = " << f.halfedge().tailVertex() << std::endl;
  }

  //EdgeData<double> negativeWeights(*globalMesh, 0.0);
  
  //find the gradient operator (want to do this just once)
  std::tie(V, F) = getVertexPositionsandFaceLists(*globalGeometry);
  igl::grad(V,F,grad);
  //find the max edge length so that default coursePeriod = 2 * max_e 
  double maxLength = -DBL_MAX;
  for (Edge e : globalMesh -> edges()){
    double length = norm(globalGeometry->vertexPositions[e.halfedge().tipVertex()] - globalGeometry->vertexPositions[e.halfedge().tailVertex()]);
    if (length > maxLength)
      maxLength = length;
  }
  H(maxLength);
  coursePeriod = 1.8 * maxLength; // let's push it :-)
  // coursePeriod = 2.0 * maxLength; //set the default period to be twice the max edge length
  // coursePeriod = 1.0 * maxLength; // skirt with seam

  // Parse stripe period, if available
  if (argc > 2) {
    sscanf(argv[2], "%f", &coursePeriod);
    knoppelFrequency = 1.0 / coursePeriod; // to match course and wale periods
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

  H(coursePeriod);

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

  //visualizing the cotan weights
  // globalGeometry->requireEdgeCotanWeights();
  // for (Edge e : globalMesh -> edges()){
  //   if (globalGeometry->edgeCotanWeights[e] < 0) negativeWeights[e] = 1.0;
  // }
  // globalPSMesh->addEdgeScalarQuantity("negative cotan weights", negativeWeights);

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

  //render the stitched vertices
  renderStitchedVertices(*globalGeometry, vertexMappingsPairs);
  // Disable the ground plane
  //polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // Set the callback function
  polyscope::state::userCallback = callBacks;

  // for rendering figures
  //   has to go before the show() call
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  // std::string viewerString = R"({"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Perspective","viewMat":[0.888747692108154,3.49245965480804e-10,0.458396464586258,0.169118270277977,0.229787424206734,0.865281701087952,-0.445515632629395,0.0717917159199715,-0.396641701459885,0.501285016536713,0.769016802310944,-2.19324994087219,0.0,0.0,0.0,1.0],"windowHeight":916,"windowWidth":1512})"; // bob
  // std::string viewerString = R"({"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Perspective","viewMat":[0.939642131328583,8.14907252788544e-10,-0.342151880264282,1.5130295753479,-0.145037829875946,0.905707657337189,-0.398314714431763,34.9959487915039,0.309889197349548,0.423899203538895,0.851041853427887,-132.280212402344,0.0,0.0,0.0,1.0],"windowHeight":916,"windowWidth":1512})"; // skirt

  // polyscope::view::setViewFromJson(viewerString, false);

  if (argc > 2) // period was provided
    showStripePatterns();


  polyscope::show();

  return EXIT_SUCCESS; 
}