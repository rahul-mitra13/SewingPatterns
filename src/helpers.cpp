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
std::vector<Halfedge> shortestEdgePathOnBoundary(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert){

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
    if (!(currVert).isBoundary()) continue;//also consider only boundary vertices on this path
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

std::tuple<std::vector<Vertex>, std::vector<Edge>, std::vector<double>> getVerticesAndEdgesInShortestEdgePathOnBoundary(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert){
    // Gather values
    SurfaceMesh& mesh = geom.mesh;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<double> weights(mesh.nEdges(), 0.0);
    std::vector<Halfedge> heList = shortestEdgePathOnBoundary(geom, startVert, endVert);
    vertices.push_back(heList[0].tailVertex());
    for (Halfedge he : heList){
        vertices.push_back(he.tipVertex());
        edges.push_back(he.edge());
        if (he.orientation()) weights[he.edge().getIndex()] = 1.0;
        else weights[he.edge().getIndex()] = -1.0;
    }
    return std::make_tuple(vertices, edges, weights);
}

std::tuple<std::vector<Vertex>, std::vector<Edge>, std::vector<double>> getVerticesAndEdgesInShortestEdgePath(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert){
    // Gather values
    SurfaceMesh& mesh = geom.mesh;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<double> weights(mesh.nEdges(), 0.0);
    std::vector<Halfedge> heList = shortestEdgePath(geom, startVert, endVert);
    vertices.push_back(heList[0].tailVertex());
    for (Halfedge he : heList){
        vertices.push_back(he.tipVertex());
        edges.push_back(he.edge());
        if (he.orientation()) weights[he.edge().getIndex()] = 1.0;
        else weights[he.edge().getIndex()] = -1.0;
    }
    return std::make_tuple(vertices, edges, weights);
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

//get the boundary edges in the wale direction where \sigma = 0
std::vector<int> getWaleBdyEdgesInGluedMesh(VertexPositionGeometry& globalGeometry, IntrinsicGeometryInterface& gluedGeometry, FaceData<Vector3>& faceGradients, std::map<int, int>& edgeMap, double val,
                                            polyscope::SurfaceMesh& globalPSMesh){
    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    std::vector<int> edgeIndicesGlobal;
    std::vector<int> edgeIndicesGlued;
    EdgeData<double> constrainedEdges(globalMesh, 0.0);
    for (Edge e : globalMesh.edges()){
        if (!e.isBoundary()) continue;
        Vector3 v1 = (globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()]);
        v1 = v1.normalize();
        Vector3 v2 = (globalGeometry.vertexPositions[e.halfedge().tailVertex()] - globalGeometry.vertexPositions[e.halfedge().tipVertex()]);
        v2 = v2.normalize();
        Vector3 grad1 = faceGradients[e.halfedge().face()];
        grad1 = grad1.normalize();
        if (dot(v1, grad1) > val || dot(v2, grad1) > val){
            edgeIndicesGlobal.push_back(e.getIndex());
            constrainedEdges[e] = 1.0;
        }
    }
    globalPSMesh.addEdgeScalarQuantity("constrained edges", constrainedEdges);

    for (int e : edgeIndicesGlobal){
        edgeIndicesGlued.push_back(edgeMap[e]);
    }
    return edgeIndicesGlued;
}

//builds a vector of "stitched together" edges from a vector of "stitched" together vertices 
std::vector<std::pair<int, int>> buildPairOfStitchedEdges(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs){
    
    SurfaceMesh& mesh = geometry.mesh;
    std::vector<std::pair<int, int>> mappedEdges;
    
    //this is so inefficient and is going to suck for larger meshes
    for (std::pair<int, int> p : vertexMappingsPairs){
        Vertex v1 = mesh.vertex(p.first);
        Vertex v2 = mesh.vertex(p.second);
        //iterate over the outgoing halfedges of v1 
        for (Halfedge v1he : v1.outgoingHalfedges()){
            //iterate over the outgoing halfedges of v2
            for (Halfedge v2he : v2.outgoingHalfedges()){
                if ((std::find(vertexMappingsPairs.begin(), vertexMappingsPairs.end(), std::make_pair(static_cast<int> (v1he.tipVertex().getIndex()), static_cast<int> (v2he.tipVertex().getIndex())))
                != vertexMappingsPairs.end()) || (std::find(vertexMappingsPairs.begin(), vertexMappingsPairs.end(), std::make_pair(static_cast<int> (v2he.tipVertex().getIndex()), 
                static_cast<int>(v1he.tipVertex().getIndex())))
                != vertexMappingsPairs.end()))
                mappedEdges.push_back(std::make_pair(v1he.edge().getIndex(), v2he.edge().getIndex()));
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

std::vector<double> createEdgeWeightsFromVertexList(VertexPositionGeometry& geometry, std::vector<int>& vertexList){

    SurfaceMesh& mesh = geometry.mesh;
    std::vector<double> edgeWeights(mesh.nEdges(), 0.0);
    for (int i = 0; i < vertexList.size(); i++){
        int v1 = vertexList[i];
        int v2 = vertexList[(i + 1) % vertexList.size()];
        for (Halfedge he : mesh.halfedges()){
            if (he.tailVertex().getIndex() == v1 && he.tipVertex().getIndex() == v2){
                if (he.orientation()) edgeWeights[he.edge().getIndex()] = 1.0;
                else edgeWeights[he.edge().getIndex()] = -1.0;
                break;
            }
        }
    }

    return edgeWeights;
}

std::unique_ptr<EdgeLengthGeometry> createGluedEdgeLengthGeometry(VertexPositionGeometry& geometry, std::vector<std::pair<int, int>>& vertexMappingsPairs, std::map<int,int>& vertexMap, 
                                    std::map<int, int>& edgeMap, std::map<int, std::vector<Halfedge>>& gluedOneRingMap){
    SurfaceMesh& mesh = geometry.mesh;

    // vertexMap will map: original mesh -> glued mesh.
    // Store another map from glued mesh index to original mesh indices
    // (it's the inverse of vertexMap)
    std::map<int, std::vector<int>> gluedMeshToOriginal;

    // Setup vertexMap and its inverse gluedMeshVertexIndexToOriginalMeshIndex
    UF uf(mesh.nVertices()); // union find
    for (auto [v1,v2] : vertexMappingsPairs)
        uf.merge(v1, v2); 
    
    // Map each root to a new vertex in the glued mesh
    int numUniqueVertices = 0;
    for (size_t i = 0; i < mesh.nVertices(); i++) {
        int root = uf.find(i);
        if (vertexMap.count(root) == 0) { // this set has not been mapped yet
            vertexMap[root] = numUniqueVertices;
            gluedMeshToOriginal[numUniqueVertices] = {root};
            numUniqueVertices++;
        } else if (i != root) { // i has a glued counterpart, but is not the root
            vertexMap[i] = vertexMap[root];
            gluedMeshToOriginal[vertexMap[i]].push_back(i);
        }
    }
    
    // Copy triangulation with new vertex indices, and create the glued mesh
    std::vector<std::vector<size_t>> polygons;
    geometry.requireEdgeLengths();
    for (Face f : mesh.faces()){
        if (f.isBoundaryLoop()) continue;
        std::vector<size_t> currPolygon;
        for (Vertex v : f.adjacentVertices())
            currPolygon.push_back(vertexMap[v.getIndex()]);
        polygons.emplace_back(currPolygon);
    }
    ManifoldSurfaceMesh * gluedMesh = new ManifoldSurfaceMesh(polygons);

    //create a map (vertex in original mesh -> outgoing global halfedges in glued context)
    for (Vertex v : mesh.vertices())
        gluedOneRingMap[v.getIndex()] = {};
    for (Halfedge he : mesh.halfedges())
        for (int iv : gluedMeshToOriginal[vertexMap[he.tailVertex().getIndex()]])
            gluedOneRingMap[iv].push_back(he);

    // Setup a map for glued mesh: (v1, v2) -> he
    std::map<std::pair<int,int>, Halfedge> vertexpair2halfedge;
    for (Halfedge he : gluedMesh->halfedges())
        vertexpair2halfedge[{he.tailVertex().getIndex(), he.tipVertex().getIndex()}] = he;

    //build a map from edges in the original mesh to edges in the glued mesh 
    //also set edge lengths in the glued mesh 
    EdgeData<double> edgeLengths(*gluedMesh);
    geometry.requireEdgeLengths();
    for (Halfedge he : mesh.halfedges()) {
        Halfedge heGlued = vertexpair2halfedge[{vertexMap[he.tailVertex().getIndex()], vertexMap[he.tipVertex().getIndex()]}];
        edgeMap[he.edge().getIndex()] = heGlued.edge().getIndex();
        edgeLengths[heGlued.edge()] = geometry.edgeLengths[he.edge()];
    }

    std::cout << "Number of faces in the original mesh " << mesh.nFaces() << std::endl;
    std::cout << "Number of faces in the glued mesh " << gluedMesh -> nFaces() << std::endl;
    std::cout << "Number of vertices in the original mesh " << mesh.nVertices() << std::endl;
    std::cout << "Number of vertices in the glued mesh " << gluedMesh -> nVertices() << std::endl;
    std::cout << "Number of edges in the original mesh " << mesh.nEdges() << std::endl;
    std::cout << "Number of edges in the glued mesh " << gluedMesh -> nEdges() << std::endl;
    std::cout << "Number of halfedges in the original mesh " << mesh.nHalfedges() << std::endl;
    std::cout << "Number of halfedges in the glued mesh " << gluedMesh -> nHalfedges() << std::endl;
    std::cout << "Number of corners in the original mesh " << mesh.nCorners() << std::endl;
    std::cout << "Number of corners in the glued mesh " << gluedMesh -> nCorners() << std::endl;
    std::cout << "Number of boundary loops in the original mesh " << mesh.nBoundaryLoops() << std::endl;
    std::cout << "Number of boundary loops in the glued mesh " << gluedMesh -> nBoundaryLoops() << std::endl;
    std::cout << "Number of connected components in the original mesh " << mesh.nConnectedComponents() << std::endl;
    std::cout << "Number of connected components ih the glued mesh " << gluedMesh -> nConnectedComponents() << std::endl;
    std::cout << "Is original mesh oriented " << mesh.isOriented() << std::endl;
    std::cout << "Is glued mesh oriented " << gluedMesh -> isOriented() << std::endl;

    //return a pointer
    return std::make_unique<EdgeLengthGeometry>(*gluedMesh, edgeLengths);
}


//parse the global boundary conditions in the glued mesh setting
//returns the boundary conditions on the glued together mesh  
globalBoundaryConditions parseJson(IntrinsicGeometryInterface& gluedGeometry, nlohmann::json& data, std::map<int,int>& vertexMap, std::map<int, int>& edgeMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    globalBoundaryConditions toReturn;

    std::vector<int> knittingStartVertices;
    std::vector<int> knittingEndVertices;

    //course conditions 
    bool courseUseBdyLoops = data["boundaries"]["course"]["useBoundaryLoops"];
    //first handle course boundary conditions
    if (courseUseBdyLoops){//just use the boundary loops of the glued mesh to specify the boundary conditions 
        std::vector<int> startVertices = data["boundaries"]["course"]["startVertices"];
        std::vector<int> endVertices = data["boundaries"]["course"]["endVertices"];
                
        //knitting start conditions
        for (int i = 0; i < startVertices.size(); i++){
            knittingStartVertices.push_back(startVertices[i]);
            //grab the boundary face
            BoundaryLoop bLoop = gluedMesh.vertex(vertexMap[startVertices[i]]).halfedge().twin().face().asBoundaryLoop();
            for (Vertex v : bLoop.adjacentVertices()){
                toReturn.courseStartBoundaryVertices.push_back(v.getIndex());
            }
            double pathLength = 0;
            for (Edge e : bLoop.adjacentEdges()){
                toReturn.courseBdyEdges.push_back(e.getIndex());
                pathLength += gluedGeometry.edgeLengths[e];
            }
            std::cout << "Knitting start boundary loop length: " << pathLength << std::endl;
            //build weights for wale boundary conditions 
            std::vector<double> weights(gluedMesh.nEdges(), 0.0);
            for (Halfedge he : bLoop.adjacentHalfedges()){
                if (he.orientation()){ // he and he.edge() point in the same direction
                    weights[he.edge().getIndex()] = 1.0;
                }
                else{ // he and he.edge() point in opposite directions
                    weights[he.edge().getIndex()] = -1.0;
                }
            }
            toReturn.waleBdyPathConstraints.push_back(weights);
        }
        //knitting end conditions
        for (int i = 0; i < endVertices.size(); i++){
            knittingEndVertices.push_back(endVertices[i]);
            //grab the boundary face 
            BoundaryLoop bLoop = gluedMesh.vertex(vertexMap[endVertices[i]]).halfedge().twin().face().asBoundaryLoop();
            for (Vertex v : bLoop.adjacentVertices()){
                toReturn.courseEndBoundaryVertices.push_back(v.getIndex());
            }
            double pathLength = 0;
            for (Edge e : bLoop.adjacentEdges()){
                toReturn.courseBdyEdges.push_back(e.getIndex());
                pathLength += gluedGeometry.edgeLengths[e];
            }
            std::cout << "Knitting end boundary loop length: " << pathLength << std::endl;
            //build weights for wale boundary conditions 
            std::vector<double> weights(gluedMesh.nEdges(), 0.0);
            for (Halfedge he : bLoop.adjacentHalfedges()){
                if (he.orientation()){
                    weights[he.edge().getIndex()] = 1.0;
                }
                else{
                    weights[he.edge().getIndex()] = -1.0;
                }
            }
            toReturn.waleBdyPathConstraints.push_back(weights);
        }
    }
    else{//can't use boundary loops for both knitting start and knitting end  
        std::vector<std::vector<int>> startVertices = data["boundaries"]["course"]["startVertices"];
        std::vector<std::vector<int>> endVertices = data["boundaries"]["course"]["endVertices"];
        //knitting start conditions
        for (int i = 0; i < startVertices.size(); i++){
            if (startVertices[i].size() == 1){//this is just a vertex on a boundary loop, do the same as above 
                knittingStartVertices.push_back(startVertices[i][0]);
                BoundaryLoop bLoop = gluedMesh.vertex(vertexMap[startVertices[i][0]]).halfedge().twin().face().asBoundaryLoop();
                for (Vertex v : bLoop.adjacentVertices()){
                    toReturn.courseStartBoundaryVertices.push_back(v.getIndex());
                }
                for (Edge e : bLoop.adjacentEdges()){
                    toReturn.courseBdyEdges.push_back(e.getIndex());
                }
                //build weights for wale boundary conditions
                std::vector<double> weights(gluedMesh.nEdges(), 0.0);
                for (Halfedge he : bLoop.adjacentHalfedges()){
                    if (he.orientation()){
                        weights[he.edge().getIndex()] = 1.0;
                    }
                    else{
                        weights[he.edge().getIndex()] = -1.0;
                    }
                }
                toReturn.waleBdyPathConstraints.push_back(weights);
            }
            else{//two vertices specify a path on the boundary 
                knittingStartVertices.push_back(startVertices[i][0]);
                Vertex startVertex = gluedMesh.vertex(vertexMap[startVertices[i][0]]);
                Vertex endVertex = gluedMesh.vertex(vertexMap[startVertices[i][1]]);
                std::vector<double> weights;
                std::vector<Vertex> vertices;
                std::vector<Edge> edges;
                std::tie(vertices, edges, weights) = getVerticesAndEdgesInShortestEdgePathOnBoundary(gluedGeometry, startVertex, endVertex);
                for (Vertex v : vertices){
                    toReturn.courseStartBoundaryVertices.push_back(v.getIndex());
                }
                for (Edge e : edges){
                    toReturn.courseBdyEdges.push_back(e.getIndex());
                }
                toReturn.waleBdyPathConstraints.push_back(weights);
            }
        }
        //knitting end conditions 
        for (int i = 0; i < endVertices.size(); i++){
            if (endVertices[i].size() == 1){//this is just a vertex on a boundary loop, do the same as above 
                knittingEndVertices.push_back(endVertices[i][0]);
                BoundaryLoop bLoop = gluedMesh.vertex(vertexMap[endVertices[i][0]]).halfedge().twin().face().asBoundaryLoop();
                for (Vertex v : bLoop.adjacentVertices()){
                    toReturn.courseEndBoundaryVertices.push_back(v.getIndex());
                }
                for (Edge e : bLoop.adjacentEdges()){
                    toReturn.courseBdyEdges.push_back(e.getIndex());
                }
                //build weights for wale boundary conditions
                std::vector<double> weights(gluedMesh.nEdges(), 0.0);
                for (Halfedge he : bLoop.adjacentHalfedges()){
                    if (he.orientation()){
                        weights[he.edge().getIndex()] = 1.0;
                    }
                    else{
                        weights[he.edge().getIndex()] = -1.0;
                    }
                }
                toReturn.waleBdyPathConstraints.push_back(weights);
            }
            else{//two vertices specify a path on the boundary 
                knittingEndVertices.push_back(endVertices[i][0]);
                Vertex startVertex = gluedMesh.vertex(vertexMap[endVertices[i][0]]);
                Vertex endVertex = gluedMesh.vertex(vertexMap[endVertices[i][1]]);
                std::vector<double> weights;
                std::vector<Vertex> vertices;
                std::vector<Edge> edges;
                std::tie(vertices, edges, weights) = getVerticesAndEdgesInShortestEdgePathOnBoundary(gluedGeometry, startVertex, endVertex);
                for (Vertex v : vertices){
                    toReturn.courseEndBoundaryVertices.push_back(v.getIndex());
                }
                for (Edge e : edges){
                    toReturn.courseBdyEdges.push_back(e.getIndex());
                }
                toReturn.waleBdyPathConstraints.push_back(weights);
            }
        }
    }
    //create boundary-boundary intergation paths 
    for (int i = 0; i < knittingStartVertices.size(); i++){
        for (int j = 0; j < knittingEndVertices.size(); j++){
            Vertex v1 = gluedMesh.vertex(vertexMap[knittingStartVertices[i]]);
            Vertex v2 = gluedMesh.vertex(vertexMap[knittingEndVertices[j]]);
            std::vector<Vertex> vertices;
            std::vector<Edge> edges;
            std::vector<double> weights; 
            std::tie(vertices, edges, weights) = getVerticesAndEdgesInShortestEdgePath(gluedGeometry, v1, v2);
            toReturn.bdyBdyPathConstraints.push_back(weights);
        }
    }

    return toReturn;
}

//convert a function defined on the vertices of the glued mesh to a function defined on the vertices of the global mesh 
VertexData<double> convertGluedToGlobalVertexFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& func, std::map<int, int>& vertexMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> globalVertexFunction(globalMesh);

    for (Vertex v : globalMesh.vertices()){
        globalVertexFunction[v] = func[gluedMesh.vertex(vertexMap[v.getIndex()])];
    }

    return globalVertexFunction; 
}

//convert a line field defined on the global mesh to a line field defined on the glued mesh 
VertexData<Vector2> convertGlobalToGluedVertexLineField(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<Vector2>& globalLineField, std::map<int, int>& vertexMap){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<Vector2> gluedLineField(gluedMesh);

    for (Vertex v : globalMesh.vertices()){
        gluedLineField[vertexMap[v.getIndex()]] = globalLineField[v];
    }

    return gluedLineField;
}

//convert function defined on the edges of the glued mesh to a function defined on the edges of the global mesh 
EdgeData<double> convertGluedToGlobalEdgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& func, std::map<int, int>& edgeMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    EdgeData<double> globalEdgeFunction(globalMesh);

    for (Edge e : globalMesh.edges()){
        globalEdgeFunction[e] = func[gluedMesh.edge(edgeMap[e.getIndex()])];
    }

    return globalEdgeFunction; 
}

//convert function defined on the edges of the global mesh to a function defined on the edges of the glued mesh 
EdgeData<double> convertGlobalToGluedEdgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& func, std::map<int, int>& edgeMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    EdgeData<double> gluedEdgeFunction(gluedMesh);
    for (Edge e : globalMesh.edges()){
        gluedEdgeFunction[edgeMap[e.getIndex()]] = func[e];
    }

    return gluedEdgeFunction;
}

//convert function defined on the halfedges of the glued mesh to a function defined on the halfedges of the global mesh 
HalfedgeData<double> convertGluedToGlobalHalfedgeFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& func, std::map<int, int>& edgeMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> globalHalfedgeFunction;
    for (Edge e : globalMesh.edges()){
        if (e.halfedge().isInterior())
            globalHalfedgeFunction[e.halfedge()] = func[gluedMesh.edge(edgeMap[e.getIndex()]).halfedge()];
        if (e.halfedge().twin().isInterior())
            globalHalfedgeFunction[e.halfedge().twin()] = func[gluedMesh.edge(edgeMap[e.getIndex()]).halfedge().twin()];
    }

    return globalHalfedgeFunction;
}

//draw the mesh curve network 
void drawMeshCurveNetwork(VertexPositionGeometry& globalGeometry, polyscope::SurfaceMesh& psMesh){

    SurfaceMesh& globalMesh = globalGeometry.mesh;

    //register the point cloud that is the mesh
    std::vector<Vector3> points(globalMesh.nVertices());  
    std::vector<std::array<int, 2>> edges(globalMesh.nEdges());
    std::vector<Vector3> edgeVectors(globalMesh.nEdges());
    globalGeometry.requireVertexPositions();
    for (Vertex v : globalMesh.vertices()){
        points[v.getIndex()] = globalGeometry.vertexPositions[v];
    }
    for (Edge e : globalMesh.edges()){
        std::array<int, 2> edge = {(int) e.halfedge().tailVertex().getIndex(), (int) e.halfedge().tipVertex().getIndex()};
        edges[e.getIndex()] = edge;
        edgeVectors[e.getIndex()] = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
    }
    auto meshNetwork = polyscope::registerCurveNetwork("Mesh curve network", points, edges);
    meshNetwork -> setRadius(0.001);
    meshNetwork -> addEdgeVectorQuantity("Canonical Edge Directions", edgeVectors);
}

//@clean
//after finding a face singularity, find a singular edge in that face 
//edge that is most aligned with the gradient (or rotated gradient in the wale case) of the time function 
//return the index of the edge in the global setting
int findSingularEdgeFromSingularFace(VertexPositionGeometry& globalGeometry, int singFaceIndex, Vector3 globalFaceGradient, double threshold, 
                                    VertexData<double>& globalTimeFunction,
                                    double isoVal){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 

    double eps = 1e-8;
    double max = -DBL_MAX;
    Edge maxEdge;
    double p1, p2;
    //normalize the global face gradient
    globalFaceGradient = globalFaceGradient.normalize();

    for (Edge e : globalMesh.face(singFaceIndex).adjacentEdges()){
        if (e.isBoundary()) continue;//skip boundary edges
        Halfedge he1 = e.halfedge();
        Halfedge he2 = e.halfedge().twin();
        double p1 = dot((globalGeometry.vertexPositions[he1.tipVertex()] - globalGeometry.vertexPositions[he1.tailVertex()]).normalize(), globalFaceGradient);
        double p2 = dot((globalGeometry.vertexPositions[he2.tipVertex()] - globalGeometry.vertexPositions[he2.tailVertex()]).normalize(), globalFaceGradient);
        if (p1 > max){
            max = p1;
            maxEdge = he1.edge();
        }
        else if (p2 > max){
            max = p2;
            maxEdge = he2.edge();
        }
    }
    //check if the selected edges are highly aligned 
    //don't take edges that already have a singularity on them 
    //or are on faces that have singularities
    if (max > threshold){
        return maxEdge.getIndex();
    }
    else{
        std::cout << "searching in adjacent faces..." << std::endl;
        //search in the 1-ring of the selected face to find an edge that is highly aligned
        //with the time function gradient
        for (Halfedge he : globalMesh.face(singFaceIndex).adjacentHalfedges()){
            if ((isoVal > globalTimeFunction[he.tailVertex()] && isoVal < globalTimeFunction[he.tipVertex()]) || 
                (isoVal > globalTimeFunction[he.tipVertex()] && isoVal < globalTimeFunction[he.tailVertex()])){//only search in halfedges that contain isoval
                Face neighborFace = he.twin().face();
                for (Edge e : neighborFace.adjacentEdges()){
                    if (e.isBoundary()) continue;//skip boundary edges
                    Halfedge he1 = e.halfedge();
                    Halfedge he2 = e.halfedge().twin();
                    double p1 = dot((globalGeometry.vertexPositions[he1.tipVertex()] - globalGeometry.vertexPositions[he1.tailVertex()]).normalize(), globalFaceGradient);
                    double p2 = dot((globalGeometry.vertexPositions[he2.tipVertex()] - globalGeometry.vertexPositions[he2.tailVertex()]).normalize(), globalFaceGradient);
                    if (p1 > max){
                        max = p1;
                        maxEdge = he1.edge();
                    }
                    else if (p2 > max){
                        max = p2;
                        maxEdge = he2.edge();
                    }
                }
            }
        }
        //check if the selected edges are highly aligned 
        //don't take edges that already have a singularity on them 
        //or are on faces that have singularities
        if (max > threshold){
            return maxEdge.getIndex();
        }
    }
    std::cout << "skipping edge selection..." << std::endl;
    //couldn't find a very well-aligned edge, skip
    //probably not the most principled thing to do here
    return -1;
}

//@clean
//hash a floating point number
int hashFloatQuantized(double f) {
    int precision = 1e5;
    int quantized = static_cast<int>(std::round(f * precision));
    return std::hash<int>{}(quantized);
}

//@clean 
//find the connected components in a time function isoline 
std::vector<std::vector<int>> findConnectedComponents(EdgeLengthGeometry& gluedGeometry, std::vector<int>& faces){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    std::vector<std::vector<int>> connectedComponents;
    FaceData<int> componentLabels(gluedMesh, -1);
    int componentIndex = 0;

    for (int val : faces){
        Face f = gluedMesh.face(val);
        if (componentLabels[f] != -1 || std::find(faces.begin(), faces.end(), f.getIndex()) == faces.end()) continue;
        // Perform BFS to find all faces in this connected component
        std::queue<Face> queue;
        queue.push(f);
        componentLabels[f] = componentIndex;

        while (!queue.empty()) {
            Face currentFace = queue.front();
            queue.pop();

            // Explore all neighboring faces across edges
            for (Halfedge he : currentFace.adjacentHalfedges()) {
                Face neighborFace = he.twin().face();
                if (componentLabels[neighborFace] == -1 && std::find(faces.begin(), faces.end(), neighborFace.getIndex()) != faces.end()) {
                    componentLabels[neighborFace] = componentIndex;
                    queue.push(neighborFace);
                }
            }
        }

        // Increment component index for the next component
        componentIndex++;
    }

    //extract the faces 
    std::unordered_map<int, std::vector<int>> components;
    for (int val : faces){
        Face f = gluedMesh.face(val);
        int label = componentLabels[f];
        if (label != -1) {
            components[label].push_back(f.getIndex());
        }
    }
    std::vector<std::vector<int>> faceComponents(componentIndex);
    for (auto component : components){
        faceComponents[component.first] = component.second;
    }
    return faceComponents;
}
    
std::vector<double> prepareCornerData(CornerData<double> cornerData) {
  std::vector<double> preparedData;
  for (Face f : cornerData.getMesh()->faces())
    for (Corner co : f.adjacentCorners())
      preparedData.push_back(cornerData[co]);
  return preparedData;
}

std::vector<int> roundWithSum(std::vector<double> y, int sum) {

    std::vector<int> res(y.size());
    
    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env);

        // Define integer variables
        std::vector<GRBVar> x;
        for (int i = 0; i < y.size(); i++)
            x.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_INTEGER));
        
        // Define objective, ||x-y||^2
        // Should we try sum of absolute values?
        GRBQuadExpr obj = 0;
        for (int i = 0; i < x.size(); i++)
            obj += (x[i]-y[i]) * (x[i]-y[i]);
        model.setObjective(obj, GRB_MINIMIZE);

        // Define sum constraint
        GRBLinExpr s;
        for (int i = 0; i < x.size(); i++)
            s += x[i];
        model.addConstr(s == sum);

        model.optimize();

        for (int i = 0; i < x.size(); i++)
            res[i] = std::round(x[i].get(GRB_DoubleAttr_X));
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }


    return res;
}