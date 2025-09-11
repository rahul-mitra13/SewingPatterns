#include <omp.h>

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
#include "KG.h"

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
double coursePeriod = 1;
double walePeriod;//the ratio of width to height should be 1:1.6 (Kui)
//threshold for constraining wale boundary edges 
//I don't really like this and need to figure out a better way of doing this
double threshold = 0.6;
//knoppel frequency in the stripe patterns
double knoppelFrequency = 0.0;

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

// Global configuration options
// TODO: setup command line interface with CLI11
Options opts; 


//here we will do as much processing as possible directly on the glued together mesh 
void showStripePatterns(){
  
  //set the wale period 
  //walePeriod = ((1.0/1.6) * coursePeriod);
  walePeriod = coursePeriod;
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

  //if doing the homology generators on the 3D models
  if (globalGeometry->mesh.nConnectedComponents() == 1){
      //std::vector<Vertex> saddleVertices = getSaddleVertices(*globalGeometry, timeFunctionGlobal);
      //allSaddleLoops = findAllSaddleLoops(*globalGeometry, saddleVertices, timeFunctionGlobal);
      homologyGenerators = buildHomologyGeneratorsVector(*globalGeometry, *globalMesh);
  }

  //get the saddle loops on general patch models
  std::vector<Vertex> saddleVertices = getSaddleVertices(*gluedELG, timeFunctionGlued);
  allSaddleLoops = findAllSaddleLoops(*gluedELG, saddleVertices, timeFunctionGlued);

  // std::cout << "size of all saddle loops = " << allSaddleLoops.size() << std::endl;


  //gradient on the glued/global mesh
  //note that faces have a 1-to-1 mapping from global to glued setting
  FaceData<Vector3> timeFunctionGradientGlobal = computeTimeFunctionFaceGrad(*globalGeometry, timeFunctionGlobal);
  //normalize the gradient of the timeFunction 
  FaceData<Vector3> timeFunctionGradientGlobalNormalized(*globalMesh);
  FaceData<Vector3> timeFunctionGradientGlobalRotated(*globalMesh);
  for (Face f : globalMesh -> faces()){
    timeFunctionGradientGlobalNormalized[f] = timeFunctionGradientGlobal[f].normalize();
    timeFunctionGradientGlobalRotated[f] = timeFunctionGradientGlobal[f].rotateAround(globalGeometry->faceNormals[f], PI/2.);
  }
  globalPSMesh -> addFaceVectorQuantity("normalized time function gradient", timeFunctionGradientGlobalNormalized);


  G = grad;

  // // Call Matteo's revealCurl
  // Model model;
  // model.setPeriod(period);
  // model.setBdyEdges(globalBdyConditions.courseBdyEdges);
  // revealCurl(*globalGeometry, *gluedELG, model, vertexMap, G);
  // P("revealCurl done");
  // polyscope::show();

  std::string richDataFile;
  // richDataFile = "info.ply";
  // richDataFile = "bunny_new_2.0.ply";
  // richDataFile = "sock_info_0.0035.ply";
  // richDataFile = "bunny_2.0_passed.ply";

  //-------iteratively find course stripes and course singularities----------//
  CornerData<double> courseStripeValues(globalGeometry -> mesh);
  EdgeData<double> courseSingularEdgesGlobal(globalGeometry -> mesh);
  FaceData<Vector3> courseOneFormGrad(globalGeometry -> mesh);
  CornerData<double> waleStripeValues;
  EdgeData<double> waleSingularEdgesGlobal;
  FaceData<int> waleSingularFaces(globalGeometry -> mesh, 0);//there are no face singularities
  if (richDataFile.size() > 0) { // if a file is specified, load it
    std::cout << "Reading a file..." << std::endl;
    RichSurfaceMeshData richData(globalGeometry->mesh, richDataFile);
    courseStripeValues = richData.getCornerProperty<double>("courseStripeValues");
    courseSingularEdgesGlobal = richData.getEdgeProperty<double>("courseSingularEdges");
    waleStripeValues = richData.getCornerProperty<double>("waleStripeValues");
    waleSingularEdgesGlobal = richData.getEdgeProperty<double>("waleSingularEdges");
  } else {
    std::tie(courseStripeValues, courseSingularEdgesGlobal) = implCourseHarmonic1Form(*globalGeometry, *gluedELG, timeFunctionGlobal,
                                                                    timeFunctionGradientGlobalNormalized, vertexMap, edgeMap, *globalPSMesh,
                                                                    globalBdyConditions, coursePeriod, V, F, G, courseOneFormGrad, gluedOneRingMap, 
                                                                    allSaddleLoops, homologyGenerators, opts);
    
    std::tie(waleStripeValues, waleSingularEdgesGlobal) = computeWaleStripeInfo(*globalGeometry, *gluedELG, 
                                                                      edgeMappingsPairs, edgeMap, vertexMap, timeFunctionGlobal, timeFunctionGlued,
                                                                      courseOneFormGrad, G, walePeriod, knoppelFrequency, globalBdyConditions, 
                                                                      courseSingularEdgesGlobal, gluedOneRingMap,*globalPSMesh, allSaddleLoops,
                                                                      homologyGenerators);

    
    // // globalPSMesh->addEdgeScalarQuantity("course singular edges", courseSingularEdgesGlobal);
    // // globalPSMesh->addEdgeScalarQuantity("wale singular edges", waleSingularEdgesGlobal);

    globalPSMesh -> addEdgeScalarQuantity("course singular edges", courseSingularEdgesGlobal);
    // Store course singular edges for rendering later
    RichSurfaceMeshData richData(globalGeometry->mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*globalGeometry);
    richData.addEdgeProperty("courseSingularEdges", courseSingularEdgesGlobal);
    richData.addCornerProperty("courseStripeValues", courseStripeValues);
    richData.addEdgeProperty("waleSingularEdges", waleSingularEdgesGlobal);
    richData.addCornerProperty("waleStripeValues", waleStripeValues);
    richData.write("info.ply");
  }

  // Draw course singular edges as a curve network. Sort them by time function order
  std::vector<Edge> courseSingularEdgePos, courseSingularEdgeNeg;
  for (Edge e : globalMesh->edges()) if (courseSingularEdgesGlobal[e] != 0) {
    if (courseSingularEdgesGlobal[e] > 0)
      courseSingularEdgePos.push_back(e);
    else
      courseSingularEdgeNeg.push_back(e);
  }
  std::sort(courseSingularEdgePos.begin(), courseSingularEdgePos.end(), [&](Edge& a, Edge& b) { return courseSingularEdgesGlobal[a] < courseSingularEdgesGlobal[b]; });
  std::sort(courseSingularEdgeNeg.begin(), courseSingularEdgeNeg.end(), [&](Edge& a, Edge& b) { return -courseSingularEdgesGlobal[a] < -courseSingularEdgesGlobal[b]; });
  std::vector<Vector3> courseSingularEdgePointsPos;
  std::vector<std::array<int,2>> courseSingularEdgesPos;
  for (Edge e : courseSingularEdgePos) {
    courseSingularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
    courseSingularEdgePointsPos.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
    courseSingularEdgesPos.push_back({(int)courseSingularEdgePointsPos.size()-2, (int)courseSingularEdgePointsPos.size()-1});
  }
  std::vector<Vector3> courseSingularEdgePointsNeg;
  std::vector<std::array<int,2>> courseSingularEdgesNeg;
  for (Edge e : courseSingularEdgeNeg) {
    courseSingularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.firstVertex()]);
    courseSingularEdgePointsNeg.push_back(globalGeometry->vertexPositions[e.secondVertex()]);
    courseSingularEdgesNeg.push_back({(int)courseSingularEdgePointsNeg.size()-2, (int)courseSingularEdgePointsNeg.size()-1});
  }
  polyscope::registerCurveNetwork("course singular edges (+1)", courseSingularEdgePointsPos, courseSingularEdgesPos)->setRadius(0.001)->setColor({0.5,0.5,0})->setEnabled(false);
  polyscope::registerCurveNetwork("course singular edges (-1)", courseSingularEdgePointsNeg, courseSingularEdgesNeg)->setRadius(0.001)->setColor({0,0.5,0.5})->setEnabled(false);


  // Plot stripe values with offset (to debug knit graph)
  FaceData<int> stripeIndicesSigmaCourse(*globalMesh, 0); // face singularities (none)
  std::vector<Vector3> positionsCourse;
  std::vector<std::array<int, 2>> edgesCourse;
  std::vector<StripeConnectedComponent> components; // don't really need this
  CornerData<double> stripeValuesWithOffset = courseStripeValues;
  for (Corner co : globalMesh->corners())
      stripeValuesWithOffset[co] -= coursePeriod/4;
  std::tie(positionsCourse, edgesCourse, components) = findStripeConnectedComponents(*globalGeometry, *gluedELG, stripeValuesWithOffset, stripeIndicesSigmaCourse, coursePeriod, edgeMap);
  polyscope::registerCurveNetwork("course stripes with offset", positionsCourse, edgesCourse)->setRadius(0.0005)->setColor({50.0/255, 205.0/255, 50.0/255})->setEnabled(false);

  // Plot wale stripe values with offset (to debug knit graph)
  CornerData<double> waleStripeValuesWithOffset(waleStripeValues);
  for (Corner co : globalGeometry->mesh.corners())
      waleStripeValuesWithOffset[co] -= walePeriod/4;
  std::vector<Vector3> positionsWale;
  std::vector<std::array<int, 2>> edgesWale;
  std::tie(positionsWale, edgesWale) = generateIsoLines(*globalGeometry, waleStripeValuesWithOffset, waleSingularFaces, walePeriod);
  polyscope::registerCurveNetwork("wale stripes with offset", positionsWale, edgesWale)->setRadius(0.0005)->setColor({1,140./255,0})->setEnabled(false);

  

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
  polyscope::registerCurveNetwork("wale singular edges (+1)", waleSingularEdgePointsPos, waleSingularEdgesPos)->setRadius(0.001)->setColor({1,0,0})->setEnabled(false);
  polyscope::registerCurveNetwork("wale singular edges (-1)", waleSingularEdgePointsNeg, waleSingularEdgesNeg)->setRadius(0.001)->setColor({0,0,1})->setEnabled(false);


  // generate the knit graph
  // graph = KnitGraph(*globalGeometry, *gluedELG, *globalPSMesh, coursePeriod, walePeriod,
  //                     courseStripeValues, courseSingularEdgesGlobal, waleStripeValues, waleSingularEdgesGlobal,
  //                     edgeMap);
  // graph.buildGraph();

  // graph.writeKnitGraphLineElement();
  // graph.writeKnitGraphToTxtFile("model.obj");

  //transfer the stripe texturing coordinates to the glued setting 
  CornerData<double> courseStripeValuesGlued(gluedELG->mesh);
  CornerData<double> waleStripeValuesGlued(gluedELG->mesh);
  for (Corner c : gluedELG->mesh.corners()){
    courseStripeValuesGlued[c] = courseStripeValues[c.getIndex()];
  }
  for (Corner c : gluedELG->mesh.corners()){
    waleStripeValuesGlued[c] = waleStripeValues[c.getIndex()];
  }

  // generate the new knit graph
  KG graph = KG(*globalGeometry, *gluedELG, *globalPSMesh, coursePeriod, walePeriod,
                      courseStripeValuesGlued, courseSingularEdgesGlobal, waleStripeValuesGlued, waleSingularEdgesGlobal,
                      edgeMap);
  graph.buildGraph();

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

  ImGui::InputDouble("Course 1-form period", &coursePeriod);
  //there needs to be a better way to constrain wale edges
  //ImGui::InputFloat("Threshold", &threshold);
  //frequency for knoppel stripes
  ImGui::InputDouble("Knoppel frequency", &knoppelFrequency);

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

  // Check how many CPUs are available
#if defined(_OPENMP)
  int nThreads = omp_get_max_threads();
  P("OpenMP is available with " << nThreads << " threads");
#else
  int nThreads = 1;
  P("OpenMP is not available!");
#endif


  polyscope::init();
  std::ifstream jsonFile(argv[1]);
  nlohmann::json data = nlohmann::json::parse(jsonFile);

  //run sanity checks
  std::tie(globalMeshPre, globalGeometryPre) = readManifoldSurfaceMesh(data["model_path"]);


  //make the mesh Delaunay 
  //fixDelaunay(*globalMeshPre, *globalGeometryPre); // we make the mesh approximately Delaunay

  // // Remesh
  // remesh(*globalMeshPre, *globalGeometryPre);

  std::tie(V, F) = getVertexPositionsandFaceLists(*globalGeometryPre);
  globalMesh = std::make_unique<ManifoldSurfaceMesh>(F);
  globalGeometry = std::make_unique<VertexPositionGeometry>(*globalMesh, V);


  //EdgeData<double> negativeWeights(*globalMesh, 0.0);
  
  //find the gradient operator (want to do this just once)
  std::tie(V, F) = getVertexPositionsandFaceLists(*globalGeometry);
  igl::grad(V,F,grad);
  //find the max edge length so that default coursePeriod = 2 * max_e 
  double maxLength = -DBL_MAX;
  double sum = 0;
  double avgEdgeLength;
  for (Edge e : globalMesh -> edges()){
    double length = norm(globalGeometry->vertexPositions[e.halfedge().tipVertex()] - globalGeometry->vertexPositions[e.halfedge().tailVertex()]);
    sum += length;
    if (length > maxLength)
      maxLength = length;
  }
  coursePeriod = 2.0 * maxLength;//set the default length to twice the course period
  avgEdgeLength = sum / globalMesh -> nEdges();

  // Parse stripe period, if available
  if (argc > 2) {
    sscanf(argv[2], "%lf", &coursePeriod);
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

  std::cout << "Number of boundary loops = " << gluedELG -> mesh.nBoundaryLoops() << std::endl;

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

  if (argc > 2) // period was provided
    showStripePatterns();

  //for rendering figures
  //has to go before the show() call
  // std::string viewerString = R"({"farClipRatio":20.0,
  //                               "fov":45.0,
  //                               "nearClipRatio":0.005,
  //                               "projectionMode":
  //                               "Perspective",
  //                               "viewMat":[0.0205091927200556,-1.87355908565223e-10,-0.999801218509674,
  //                               -0.399929404258728,-0.0461873635649681,0.998910784721375,
  //                               -0.000946979271247983,-0.528441607952118,0.99871826171875,
  //                               0.0461944416165352,0.0204852763563395,-3.25812482833862,0.0,0.0,0.0,1.0],
  //                               "windowHeight":1440,"windowWidth":2560})";
  // polyscope::view::setViewFromJson(viewerString, false);
  
  polyscope::show();

  return EXIT_SUCCESS; 
}