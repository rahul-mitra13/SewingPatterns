#include "GUI_helpers.h"


//get a pair of vertices where the shortest edge path between them specifies the boundary conditions
std::pair<Vertex, Vertex> getAndRenderUserSpecifiedBoundaryVertices(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, int timeVal){
    
    polyscope::warning("Please select 2 boundary vertices for " + std::to_string(timeVal) + " boundary conditions");

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

