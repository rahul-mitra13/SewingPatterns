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