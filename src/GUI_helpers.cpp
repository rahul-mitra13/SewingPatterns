#include "GUI_helpers.h"


//get a pair of vertices where the shortest edge path between them specifies the boundary conditions
std::pair<Vertex, Vertex> getAndRenderUserSpecifiedBoundaryVertices(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, int timeVal){
    std::vector<Vector3> pointCloud(2);
    SurfaceMesh& mesh = geometry.mesh;
    Vertex v1 = mesh.vertex(psMesh.selectVertex()); 
    Vertex v2 = mesh.vertex(psMesh.selectVertex());
    pointCloud[0] = geometry.vertexPositions[v1];
    pointCloud[1] = geometry.vertexPositions[v2];
    auto psCloud = polyscope::registerPointCloud(std::to_string(timeVal) + " boundary vertices", pointCloud);
    psCloud->setPointColor(polyscope::render::RGB_RED);
    return std::make_pair(v1, v2);
}

std::pair<std::vector<Vertex>, std::vector<Edge>> getAndRenderUserSpecifiedBoundaryInfo(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMeshes, int timeVal){

    SurfaceMesh& mesh = geometry.mesh;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Vector3> pointCloud;

    static int numComponents = 0;
    // Register the callback which creates the UI and does the hard work
    auto focusedPopupUI = [&]() {
      { // Create a window with instruction and a close button.
        static bool showWindow = true;
        ImGui::SetNextWindowSize(ImVec2(500, 0), ImGuiCond_Once);
        ImGui::Begin("Select boundary components", &showWindow);

        ImGui::PushItemWidth(300);
        ImGui::TextUnformatted(("Enter number of " + std::to_string(timeVal) + "-boundary components").c_str());
        ImGui::Separator();

        // Choose by number
        ImGui::PushItemWidth(300);
        ImGui::InputInt("number of components", &numComponents);
        ImGui::PopItemWidth();
        ImGui::Separator();
        if (ImGui::Button("Ok")) {
          polyscope::popContext();
        }
        ImGui::End();
      }
    };
    // Pass control to the context we just created
    polyscope::pushContext(focusedPopupUI);
  
    for (int i = 0; i < numComponents; i++){
      Vertex startVertex, endVertex;
      std::tie(startVertex, endVertex) = getAndRenderUserSpecifiedBoundaryVertices(geometry, psMeshes, timeVal);
      std::vector<Halfedge> heList = shortestEdgePathOnBoundary(geometry, startVertex, endVertex);
      //add the first vertex 
      pointCloud.push_back(geometry.vertexPositions[heList[0].tailVertex()]);
      vertices.push_back(heList[0].tailVertex());
      for (Halfedge he : heList){
        pointCloud.push_back(geometry.vertexPositions[he.tipVertex()]); 
        vertices.push_back(he.tipVertex());
        edges.push_back(he.edge());
      }
    }
    auto cloud = polyscope::registerPointCloud("boundary " + std::to_string(timeVal) + " vertices", pointCloud);
    return std::make_pair(vertices, edges);
}

PatchBoundaryConditions getAndRenderBoundaryInfoFromJson(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMeshes, const std::string& jsonFilePath){

  SurfaceMesh& mesh = geometry.mesh;
  PatchBoundaryConditions toReturn;

  std::ifstream jsonFile(jsonFilePath);
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  std::vector<std::vector<int>> courseStartBoundary = data["boundaries"]["course"]["start"];
  std::vector<std::vector<int>> courseEndBoundary = data["boundaries"]["course"]["end"];
  std::vector<std::vector<int>> waleStartBoundary = data["boundaries"]["wale"]["start"];
  std::vector<std::vector<int>> waleEndBoundary = data["boundaries"]["wale"]["end"];
  
  //handle course start case
  std::vector<Vector3> courseStartVerticesPointCloud;
  for (int i = 0; i < courseStartBoundary.size(); i++){
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    Vertex startVertex = mesh.vertex(courseStartBoundary[i][0]);
    Vertex endVertex = mesh.vertex(courseStartBoundary[i][1]);
    std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(geometry, startVertex, endVertex);
    for (Vertex v : vertices){
     courseStartVerticesPointCloud.push_back(geometry.vertexPositions[v]);
    }
    toReturn.courseStartBoundaryVertices.insert(toReturn.courseStartBoundaryVertices.end(), vertices.begin(), vertices.end());
    toReturn.courseStartBoundaryEdges.insert(toReturn.courseStartBoundaryEdges.end(), edges.begin(), edges.end());
  }
  auto courseStartPC = polyscope::registerPointCloud("course start vertices", courseStartVerticesPointCloud);
  courseStartPC->setPointColor(polyscope::render::RGB_RED);

  //handle course end case
  std::vector<Vector3> courseEndVerticesPointCloud;
  for (int i = 0; i < courseEndBoundary.size(); i++){
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (int j = 0; j < 2; j++){
      Vertex startVertex = mesh.vertex(courseEndBoundary[i][0]);
      Vertex endVertex = mesh.vertex(courseEndBoundary[i][1]);
      std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(geometry, startVertex, endVertex);
    }
    for (Vertex v : vertices){
     courseEndVerticesPointCloud.push_back(geometry.vertexPositions[v]);
    }
    toReturn.courseEndBoundaryVertices.insert(toReturn.courseEndBoundaryVertices.end(), vertices.begin(), vertices.end());
    toReturn.courseEndBoundaryEdges.insert(toReturn.courseEndBoundaryEdges.end(), edges.begin(), edges.end());
  }
  auto courseEndPC = polyscope::registerPointCloud("course end vertices", courseEndVerticesPointCloud);
  courseEndPC->setPointColor(polyscope::render::RGB_RED);

  //handle wale start case
  std::vector<Vector3> waleStartVerticesPointCloud;
  for (int i = 0; i < waleStartBoundary.size(); i++){
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (int j = 0; j < 2; j++){
      Vertex startVertex = mesh.vertex(waleStartBoundary[i][0]);
      Vertex endVertex = mesh.vertex(waleStartBoundary[i][1]);
      std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(geometry, startVertex, endVertex);
    }
    for (Vertex v : vertices){
     waleStartVerticesPointCloud.push_back(geometry.vertexPositions[v]);
    }
    toReturn.waleStartBoundaryVertices.insert(toReturn.waleStartBoundaryVertices.end(), vertices.begin(), vertices.end());
    toReturn.waleStartBoundaryEdges.insert(toReturn.waleStartBoundaryEdges.end(), edges.begin(), edges.end());
  }
  auto waleStartPC = polyscope::registerPointCloud("wale start vertices", waleStartVerticesPointCloud);
  waleStartPC->setPointColor(polyscope::render::RGB_BLUE);

  //handle wale end case
  std::vector<Vector3> waleEndVerticesPointCloud;
  for (int i = 0; i < waleEndBoundary.size(); i++){
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (int j = 0; j < 2; j++){
      Vertex startVertex = mesh.vertex(waleEndBoundary[i][0]);
      Vertex endVertex = mesh.vertex(waleEndBoundary[i][1]);
      std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(geometry, startVertex, endVertex);
    }
    for (Vertex v : vertices){
     waleEndVerticesPointCloud.push_back(geometry.vertexPositions[v]);
    }
    toReturn.waleEndBoundaryVertices.insert(toReturn.waleEndBoundaryVertices.end(), vertices.begin(), vertices.end());
    toReturn.waleEndBoundaryEdges.insert(toReturn.waleEndBoundaryEdges.end(), edges.begin(), edges.end());
  }
  auto waleEndPC = polyscope::registerPointCloud("wale end vertices", waleEndVerticesPointCloud);
  waleEndPC->setPointColor(polyscope::render::RGB_BLUE);

  return toReturn;

}

void getAndRenderBoundaryInfoFromJson(const std::string& jsonFilePath, std::vector<std::unique_ptr<ManifoldSurfaceMesh>>& meshes,
                                      std::vector<std::unique_ptr<VertexPositionGeometry>>& geometries,
                                      std::vector<PatchBoundaryConditions>& bdyConditions){

                                        
  std::ifstream jsonFile(jsonFilePath);
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  auto numPanels = data["panels"].size();
  //resize everything 
  meshes.resize(numPanels);
  geometries.resize(numPanels);
  bdyConditions.resize(numPanels);
  for (int i = 0; i < numPanels; i++){
    auto obj_path = data["panels"][std::to_string(i)]["obj"];
    std::tie(meshes[i], geometries[i]) = readManifoldSurfaceMesh(obj_path);
    auto courseStartBoundary = data["panels"][std::to_string(i)]["boundaries"]["course"]["start"];
    auto courseEndBoundary = data["panels"][std::to_string(i)]["boundaries"]["course"]["end"];
    std::cout << "size of course start boundary = " << courseStartBoundary.size() << std::endl;
    std::cout << "size of course end boundary = " << courseEndBoundary.size() << std::endl;
    //handle course constraints
    for (int j = 0; j < courseStartBoundary.size(); j++){
      Vertex v1 = meshes[i]->vertex(courseStartBoundary[j][0]);
      Vertex v2 = meshes[i]->vertex(courseStartBoundary[j][1]);
      std::vector<Vertex> vertices;
      std::vector<Edge> edges;
      std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(*geometries[i], v1, v2);
      bdyConditions[i].courseStartBoundaryVertices.insert(bdyConditions[i].courseStartBoundaryVertices.end(), vertices.begin(), vertices.end());
      bdyConditions[i].courseStartBoundaryEdges.insert(bdyConditions[i].courseStartBoundaryEdges.end(), edges.begin(), edges.end());
    }
    for (int j = 0; j < courseEndBoundary.size(); j++){
      Vertex v1 = meshes[i]->vertex(courseEndBoundary[j][0]);
      Vertex v2 = meshes[i]->vertex(courseEndBoundary[j][1]);
      std::vector<Vertex> vertices;
      std::vector<Edge> edges;
      std::tie(vertices, edges) = getVerticesAndEdgesInShortestEdgePathOnBoundary(*geometries[i], v1, v2);
      bdyConditions[i].courseEndBoundaryVertices.insert(bdyConditions[i].courseEndBoundaryVertices.end(), vertices.begin(), vertices.end());
      bdyConditions[i].courseEndBoundaryEdges.insert(bdyConditions[i].courseEndBoundaryEdges.end(), edges.begin(), edges.end());
    }
  }
}