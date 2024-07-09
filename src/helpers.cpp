#include "helpers.h"

//get a pair of vertex lists that are at the extremes of a panel
std::pair<std::vector<Vertex>, std::vector<Vertex>> getBoundaryVertices(VertexPositionGeometry& geometry, int axis){
    
    SurfaceMesh& mesh = geometry.mesh;
    double eps = 1e-4;
    double min = DBL_MAX;
    double max = DBL_MIN;
    std::vector<Vertex> lowestVertices;
    std::vector<Vertex> highestVertices;

    //find the lowest and highest axis value
    for (Vertex v : mesh.vertices()){
        Vector3 pos = geometry.vertexPositions[v];
        if (pos[axis] < min){
            min = pos[axis];
        }
        if (pos[axis] > max){
            max = pos[axis];
        }
    }
    //find all vertices which are close to those values
    for (Vertex v : mesh.vertices()){
        if (std::fabs(geometry.vertexPositions[v][axis] - min) < eps){
            lowestVertices.push_back(v);
        }
        if (std::fabs(geometry.vertexPositions[v][axis] - max) < eps){
            highestVertices.push_back(v);
        }
    }

    return std::make_pair(lowestVertices, highestVertices);
}

//project a vector onto a given plane
Vector2 projectOntoPlane(const Eigen::Vector3d &vec, const Eigen::Vector3d &normal, const Eigen::Vector3d &axis){
    Eigen::Vector3d eY = normal.cross(axis).normalized();
    Eigen::Vector3d eX = eY.cross(normal).normalized();
    Eigen::Vector2d toReturn =  {eX.dot(vec), eY.dot(vec)};
    return Vector2{toReturn(0), toReturn(1)};
}


//get a pair of boundary edges that are at the extremes of a panel
std::pair<std::vector<Edge>, std::vector<Edge>> getBoundaryEdges(VertexPositionGeometry& geometry, int axis){

    SurfaceMesh& mesh = geometry.mesh;

    std::vector<Edge> lowestEdges;
    std::vector<Edge> highestEdges;
    
    std::vector<Vertex> lowestVertices = getBoundaryVertices(geometry, axis).first;
    std::vector<Vertex> highestVertices = getBoundaryVertices(geometry, axis).second;

    for (Vertex v : lowestVertices){
        for (Halfedge he : v.outgoingHalfedges()){
            if (std::find(lowestVertices.begin(), lowestVertices.end(), he.tipVertex()) != lowestVertices.end()){//the outgoing halfedge points to a vertex in the set
                if (std::find(lowestEdges.begin(), lowestEdges.end(), he.edge()) == lowestEdges.end())
                    lowestEdges.push_back(he.edge());
            }
        } 
    }
    //do the same for highest vertices 
    for (Vertex v : highestVertices){
        for (Halfedge he : v.outgoingHalfedges()){
            if (std::find(highestVertices.begin(), highestVertices.end(), he.tipVertex()) != highestVertices.end()){//the outgoing halfedge points to a vertex in the set
                if (std::find(highestEdges.begin(), highestEdges.end(), he.edge()) == highestEdges.end())
                    highestEdges.push_back(he.edge());
            }
        } 
    }

    return std::make_pair(lowestEdges, highestEdges);

}

// Get vertex positions and face lists from an input geometry
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> getVertexPositionsandFaceLists(VertexPositionGeometry& geometry){

    SurfaceMesh& mesh = geometry.mesh;
    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    for (Vertex v : mesh.vertices()){
        Vector3 pos = geometry.vertexPositions[v];
        V(v.getIndex(), 0) = pos.x;
        V(v.getIndex(), 1) = pos.y;
        V(v.getIndex(), 2) = pos.z;
    }

    std::vector<std::vector<size_t>> faceVertexIndices = mesh.getFaceVertexList();

    for (int i = 0; i < faceVertexIndices.size(); i++){
        F(i, 0) = faceVertexIndices[i][0];
        F(i, 1) = faceVertexIndices[i][1];
        F(i, 2) = faceVertexIndices[i][2];
    }

    return std::make_pair(V, F);
}

//get shortest edge path between two vertices by taking only boundary edges of the mesh 
std::vector<Halfedge> shortestEdgePathOnBoundary(VertexPositionGeometry& geom, Vertex startVert, Vertex endVert){

    // Early out for empty case
  if (startVert == endVert) {
    return std::vector<Halfedge>();
  }

  // Gather values
  SurfaceMesh& mesh = geom.mesh;
  geom.requireEdgeLengths();

  // Search state: incoming halfedges to each vertex, once discovered
  std::unordered_map<Vertex, Halfedge> incomingHalfedge;

  // Search state: visible neighbors eligible to expand to
  using WeightedHalfedge = std::tuple<double, Halfedge>;
  std::priority_queue<WeightedHalfedge, std::vector<WeightedHalfedge>, std::greater<WeightedHalfedge>> pq;

  // Helper to add a vertex's
  auto vertexDiscovered = [&](Vertex v) {
    return v == startVert || incomingHalfedge.find(v) != incomingHalfedge.end();
  };
  auto enqueueVertexNeighbors = [&](Vertex v, double dist) {
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!vertexDiscovered(he.twin().vertex())) {
        double len = geom.edgeLengths[he.edge()];
        double targetDist = dist + len;
        pq.emplace(targetDist, he);
      }
    }
  };

  // Add initial halfedges
  enqueueVertexNeighbors(startVert, 0.);

  while (!pq.empty()) {

    // Get the next closest neighbor off the queue
    double currDist = std::get<0>(pq.top());
    Halfedge currIncomingHalfedge = std::get<1>(pq.top());
    pq.pop();

    Vertex currVert = currIncomingHalfedge.twin().vertex();
    if (vertexDiscovered(currVert)) continue;
    if (!(currVert).isBoundary()) continue;//also consider boundary vertices on this path
    // Accept the neighbor
    incomingHalfedge[currVert] = currIncomingHalfedge;

    // Found path! Walk backwards to reconstruct it and return
    if (currVert == endVert) {
      std::vector<Halfedge> path;
      Vertex walkV = currVert;
      while (walkV != startVert) {
        Halfedge prevHe = incomingHalfedge[walkV];
        path.push_back(prevHe);
        walkV = prevHe.vertex();
      }

      std::reverse(std::begin(path), std::end(path));

      geom.unrequireEdgeLengths();
      return path;
    }

    // Enqueue neighbors
    enqueueVertexNeighbors(currVert, currDist);
  }

  // Didn't find path
  geom.unrequireEdgeLengths();
  return std::vector<Halfedge>();
}

std::pair<std::vector<Vertex>, std::vector<Edge>> getVerticesAndEdgesInShortestEdgePathOnBoundary(VertexPositionGeometry& geom, Vertex startVert, Vertex endVert){
    // Gather values
    SurfaceMesh& mesh = geom.mesh;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Halfedge> heList = shortestEdgePathOnBoundary(geom, startVert, endVert);
    vertices.push_back(heList[0].tailVertex());
    for (Halfedge he : heList){
        vertices.push_back(he.tipVertex());
        edges.push_back(he.edge());
    }
    return std::make_pair(vertices, edges);
}

//this builds a mapping
//panel1name -> pair(panel2name, (vertex from panel1 -> vertex from panel 2))
std::map<std::pair<std::string, int>, std::pair<std::string, int>> buildLocalVertexMappingMapFromFile(const std::string& filename){
    
    std::map<std::pair<std::string, int>, std::pair<std::string, int>> vertexMappings;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return vertexMappings;
    }
    std::string line;
    while (std::getline(file, line)) {
        // Remove parentheses
        line.erase(std::remove(line.begin(), line.end(), '('), line.end());
        line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> parsedLine;
        while (std::getline(ss, item, ',')) {
            parsedLine.push_back(item);
        }
        vertexMappings.insert({std::make_pair(parsedLine[0], std::stoi(parsedLine[1])), std::make_pair(parsedLine[2], std::stoi(parsedLine[3]))});
    }
    return vertexMappings;
}

//this builds a mapping of stitched together vertices
//(vertex1 -> vertex2)
std::map<int, int> buildGlobalVertexMappingFromFile(const std::string& filename){
    
    std::map<int, int> vertexMappings;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return vertexMappings;
    }
    std::string line;
    while (std::getline(file, line)) {
        // Remove parentheses
        line.erase(std::remove(line.begin(), line.end(), '('), line.end());
        line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> parsedLine;
        while (std::getline(ss, item, ',')) {
            parsedLine.push_back(item);
        }
        vertexMappings.insert({std::stoi(parsedLine[0]), std::stoi(parsedLine[1])});
    }

    return vertexMappings;

}

//returns a vector of pairs of vertex mappings
//(vertex1, vertex2)
std::vector<std::pair<int, int>> buildPairOfStitchedVerticesFromFile(const std::string& filename){

    std::vector<std::pair<int, int>> vertexMappings;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return vertexMappings;
    }
    std::string line;
    while (std::getline(file, line)) {
        // Remove parentheses
        line.erase(std::remove(line.begin(), line.end(), '('), line.end());
        line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> parsedLine;
        while (std::getline(ss, item, ',')) {
            parsedLine.push_back(item);
        }
        vertexMappings.push_back(std::make_pair(std::stoi(parsedLine[0]), std::stoi(parsedLine[1])));
    }

    return vertexMappings;
}

//render the stitched together vertices in polyscope
void renderStitchedVertices(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs){

    SurfaceMesh& mesh = geometry.mesh;

    std::vector<Vector3> nodes;
    std::vector<std::array<int, 2>> edges;
    int index = 0;
    for (auto entry: vertexMappingsPairs){
        nodes.push_back(geometry.vertexPositions[mesh.vertex(entry.first)]);
        nodes.push_back(geometry.vertexPositions[mesh.vertex(entry.second)]);
        std::array<int, 2> edge = {index, index + 1};
        edges.push_back(edge);
        index += 2;
    }
    auto stitches = polyscope::registerCurveNetwork("stitched vertices", nodes, edges);
    stitches->setRadius(0.001);
}

//builds a vector of "stitched together" edges from a vector of "stitched" together vertices 
std::vector<std::pair<int, int>> buildPairOfStitchedEdges(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs){
    
    SurfaceMesh& mesh = geometry.mesh;
    std::vector<std::pair<int, int>> mappedEdges;
    int ctr = 0;

    //this is so inefficient and is going to suck for larger meshes
    for (std::pair<int, int> p : vertexMappingsPairs){
        Vertex v1 = mesh.vertex(p.first);
        Vertex v2 = mesh.vertex(p.second);
        //iterate over the outgoing halfedges of v1 
        for (Halfedge v1he : v1.outgoingHalfedges()){
            //iterate over the outgoing halfedges of v2
            for (Halfedge v2he : v2.outgoingHalfedges()){
                if ((std::find(vertexMappingsPairs.begin(), vertexMappingsPairs.end(), std::make_pair(v1he.tipVertex().getIndex(), v2he.tipVertex().getIndex()))
                != vertexMappingsPairs.end()) || (std::find(vertexMappingsPairs.begin(), vertexMappingsPairs.end(), std::make_pair(v2he.tipVertex().getIndex(), v1he.tipVertex().getIndex()))
                != vertexMappingsPairs.end()))
                mappedEdges.push_back(std::make_pair(v1he.edge().getIndex(), v2he.edge().getIndex()));
                ctr++;
                break;
            }
        }
    }
    return mappedEdges;
}

std::vector<int> creatEdgeListFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList){

    SurfaceMesh& mesh = geometry.mesh;
    std::vector<int> edgeList;

    for (Edge e : mesh.edges()){
        if (std::find(vertexList.begin(), vertexList.end(), e.halfedge().tailVertex().getIndex()) != vertexList.end() && 
        std::find(vertexList.begin(), vertexList.end(), e.halfedge().tipVertex().getIndex()) != vertexList.end()){
            edgeList.push_back(e.getIndex());
        }
    }

    return edgeList;
}

SurfaceMesh * createGluedSurfaceMesh(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs, std::map<int,int>& originalMeshVertexIndexToGluedMeshIndex,
                                        std::map<int, int>& originalMeshEdgeIndexToGluedMeshIndex, std::vector<double>& gluedMeshEdgeLengths){
    SurfaceMesh& mesh = geometry.mesh;
    //find original index to glued mesh index for vertices
    int numUniqueVertices = 0;
    for (Vertex v : mesh.vertices()){
        size_t iV = v.getIndex();
        //iterate over the mappings 
        for (auto p : vertexMappingsPairs){
            if (p.first == iV || p.second == iV){
                if (originalMeshVertexIndexToGluedMeshIndex.find(p.first) == originalMeshVertexIndexToGluedMeshIndex.end()
                && originalMeshVertexIndexToGluedMeshIndex.find(p.second) == originalMeshVertexIndexToGluedMeshIndex.end()){
                    originalMeshVertexIndexToGluedMeshIndex.insert({p.first, numUniqueVertices});
                    originalMeshVertexIndexToGluedMeshIndex.insert({p.second, numUniqueVertices});
                    numUniqueVertices++;
                }
                if (originalMeshVertexIndexToGluedMeshIndex.find(p.first) != originalMeshVertexIndexToGluedMeshIndex.end()&&
                    originalMeshVertexIndexToGluedMeshIndex.find(p.second) == originalMeshVertexIndexToGluedMeshIndex.end()){
                    originalMeshVertexIndexToGluedMeshIndex.insert({p.second, originalMeshVertexIndexToGluedMeshIndex.at(p.first)});
                }
                if (originalMeshVertexIndexToGluedMeshIndex.find(p.second) != originalMeshVertexIndexToGluedMeshIndex.end() &&
                    originalMeshVertexIndexToGluedMeshIndex.find(p.first) == originalMeshVertexIndexToGluedMeshIndex.end()){
                    originalMeshVertexIndexToGluedMeshIndex.insert({p.first, originalMeshVertexIndexToGluedMeshIndex.at(p.second)});
                }
            }
        }
        if (originalMeshVertexIndexToGluedMeshIndex.find(iV) == originalMeshVertexIndexToGluedMeshIndex.end()){
            originalMeshVertexIndexToGluedMeshIndex.insert({iV, numUniqueVertices});
            numUniqueVertices++;
        }
    }
    std::vector<std::vector<size_t>> polygons;
    geometry.requireEdgeLengths();
    for (Face f : mesh.faces()){
        if (f.isBoundaryLoop()) continue;
        std::vector<size_t> currPolygon;
        int i = originalMeshVertexIndexToGluedMeshIndex[f.halfedge().tailVertex().getIndex()];
        int j = originalMeshVertexIndexToGluedMeshIndex[f.halfedge().next().tailVertex().getIndex()];
        int k = originalMeshVertexIndexToGluedMeshIndex[f.halfedge().next().next().tailVertex().getIndex()];
        // int i = f.halfedge().tailVertex().getIndex();
        // int j = f.halfedge().next().tailVertex().getIndex();
        // int k = f.halfedge().next().next().tailVertex().getIndex();
        currPolygon.push_back(i);
        currPolygon.push_back(j);
        currPolygon.push_back(k);
        polygons.emplace_back(currPolygon);
    }
    //for some reason this is not manifold and has duplicate edges? 
    SurfaceMesh * gluedMesh = new SurfaceMesh(polygons);
    std::cout << "isOriented? " << gluedMesh->isOriented() << std::endl;
    std::cout << "isEdgeManifold? " << gluedMesh->isEdgeManifold() << std::endl;
    std::cout << "isManifold? " << gluedMesh->isManifold() << std::endl;

    std::cout << "Number of boundary loops in the original mesh " << mesh.nBoundaryLoops() << std::endl;
    std::cout << "Number of boundary loops in the glued mesh " << gluedMesh->nBoundaryLoops() << std::endl;
    return gluedMesh;
}