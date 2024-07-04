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

//build the global cotan Laplacian for all the panels while accounting for the mapped stitches across panels
//reads the local mappings
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::map<int, int>& vertexMappings, globalBoundaryConditions& bdyConditions){
    
    SurfaceMesh& mesh = geometry.mesh;
    int numVertices = mesh.nVertices(); 
    int numStitchedVertices = vertexMappings.size();
    Eigen::SparseMatrix<double> L(numVertices - numStitchedVertices, numVertices - numStitchedVertices);
    //require the cotan edge weights 
    geometry.requireEdgeCotanWeights();

    //store mappings between index in the original mesh and index in the Laplacian matrix
    std::map<int, int> originalIndexToLaplacianMatrixIndex;
    //store mappings between index in the Laplacian matrix to index in the orignal mesh 
    std::map<int, int> laplacianMatrixIndexToOriginalIndex;
    //number of "unique vertices" i.e., only consider one vertex per stitch
    int numUniqueVertices = 0;
    std::vector<int> seenVertices;
    for (Vertex v : mesh.vertices()){
        size_t iV = v.getIndex();
        if (vertexMappings.find(iV) != vertexMappings.end()){
            seenVertices.push_back(iV);
            seenVertices.push_back(vertexMappings.at(iV));

            //emulating bi-directional mapping here
            originalIndexToLaplacianMatrixIndex.insert({iV, numUniqueVertices});
            originalIndexToLaplacianMatrixIndex.insert({vertexMappings.at(iV), numUniqueVertices});
            laplacianMatrixIndexToOriginalIndex.insert({numUniqueVertices, iV});
            numUniqueVertices++;
        }
        else if (originalIndexToLaplacianMatrixIndex.find(iV) == originalIndexToLaplacianMatrixIndex.end()){//we've never seen this index before
            originalIndexToLaplacianMatrixIndex.insert({iV, numUniqueVertices});
            laplacianMatrixIndexToOriginalIndex.insert({numUniqueVertices, iV});
            numUniqueVertices++;
        }  
    }
    std::vector<Eigen::Triplet<double>> tripletList;

    //keep a set of indices you've already populated in the Laplacian 
    std::set<int> setIndices;//set of vertex indices we've already set in the laplacian

    for (Vertex v : mesh.vertices()){
        if (std::find(setIndices.begin(), setIndices.end(), originalIndexToLaplacianMatrixIndex.at(v.getIndex())) != setIndices.end()) continue;//we've handled this already
        double L_diag = 0.0;//diagonal entries of L
        //iterate over the one-ring of the vertex 
        for (Halfedge he : v.outgoingHalfedges()){
            //off diagonal entries 
            tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(he.tipVertex().getIndex()),
            -geometry.edgeCotanWeights[he.edge()]);
            setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));

            if (vertexMappings.find(v.getIndex()) != vertexMappings.end()){//handle off diagonal entry for stitched vertex
                Vertex mappedVertex = mesh.vertex(vertexMappings.at(v.getIndex()));
                //iterate over the 1-ring of the mapped vertex 
                for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                    tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), 
                                        originalIndexToLaplacianMatrixIndex.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                    setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
                }
            }

            //handle the diagonal entry case
            L_diag += geometry.edgeCotanWeights[he.edge()];
            if (vertexMappings.find(v.getIndex()) != vertexMappings.end()){//handle diagonal entry for stitched vertex
                Vertex mappedVertex = mesh.vertex(vertexMappings.at(v.getIndex()));
                //iterate over the 1-ring of the mapped vertex 
                for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                    L_diag += geometry.edgeCotanWeights[mappedVertexHalfedge.edge()];
                }
            }
        }
        tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(v.getIndex()), L_diag);
        setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
    }

    L.setFromTriplets(tripletList.begin(), tripletList.end());

    //force boundary conditions
    Eigen::VectorXd b = Eigen::VectorXd::Zero(numUniqueVertices);
    
    for (Vertex v : bdyConditions.courseStartBoundaryVertices){
        int updatedIndex = originalIndexToLaplacianMatrixIndex.at(v.getIndex());
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 0.0;
    }

    for (Vertex v : bdyConditions.courseEndBoundaryVertices){
        int updatedIndex = originalIndexToLaplacianMatrixIndex.at(v.getIndex());
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 1.0;
    }
    
    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed" << std::endl;
    }
    Eigen::VectorXd u = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed" << std::endl;
    }
    
    VertexData<double> testFunction(mesh);
    for (Vertex v : mesh.vertices()){
        testFunction[v] = u(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
    }

    return testFunction;
}

//compute time function using a vector of pairs of vertex mappings instead of the map because we miss stitches then 
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::vector<std::pair<int,int>>& vertexMappingsPairs){

    SurfaceMesh& mesh = geometry.mesh;
    VertexData<double> timeFunction(mesh);
    int numVertices = mesh.nVertices(); 
    int numStitchedVertices = vertexMappingsPairs.size();
    Eigen::SparseMatrix<double> L(numVertices - numStitchedVertices, numVertices - numStitchedVertices);
    //require the cotan edge weights 
    geometry.requireEdgeCotanWeights();
    //store mappings between index in the original mesh and index in the Laplacian matrix
    std::map<int, int> originalIndexToLaplacianMatrixIndex;
    //store mappings between index in the Laplacian matrix to index in the orignal mesh 
    std::map<int, int> laplacianMatrixIndexToOriginalIndex;
    //all the pairs that have been seen so far
    std::vector<std::pair<int, int>> seenPairs;
    //number of "unique vertices" i.e., only consider one vertex per stitch
    int numUniqueVertices = 0;
    /**
    for (Vertex v : mesh.vertices()){
        size_t iV = v.getIndex();
        bool isMappedVertex = false;
        std::pair<int, int> pairing;
        //search for mapping
        for (auto p : vertexMappingsPairs){
            if ((p.first == iV) && (std::find(seenPairs.begin(), seenPairs.end(), p) == seenPairs.end())){//this vertex has a mapping and we have not seen this pair before
                isMappedVertex = true;
                pairing = p;
                seenPairs.push_back(p);
                break;
            }
        }
        if (isMappedVertex){
            //emulating a bidirectional map
            originalIndexToLaplacianMatrixIndex.insert({iV, numUniqueVertices});
            originalIndexToLaplacianMatrixIndex.insert({pairing.second, numUniqueVertices});
            laplacianMatrixIndexToOriginalIndex.insert({numUniqueVertices, iV});
            numUniqueVertices++;
        }
        else if (originalIndexToLaplacianMatrixIndex.find(iV) == originalIndexToLaplacianMatrixIndex.end()){//we've never seen this index before
            originalIndexToLaplacianMatrixIndex.insert({iV, numUniqueVertices});
            laplacianMatrixIndexToOriginalIndex.insert({numUniqueVertices, iV});
            numUniqueVertices++;
        }
    }
    std::cout << "Number of vertices " << numVertices << std::endl;
    std::cout << "Number of stitches " << numStitchedVertices << std::endl;
    std::cout << "Number of unique vertices " << numUniqueVertices << std::endl;

    for (auto entry : originalIndexToLaplacianMatrixIndex){
        std::cout << entry.first << " mapped to " << entry.second << std::endl;
    }
    seenPairs.clear();
    std::vector<Eigen::Triplet<double>> tripletList;
    //keep a set of indices you've already populated in the Laplacian 
    std::set<int> setIndices;//set of vertex indices we've already set in the laplacian
    for (Vertex v : mesh.vertices()){
        size_t iV = v.getIndex();
        if (std::find(setIndices.begin(), setIndices.end(), originalIndexToLaplacianMatrixIndex.at(v.getIndex())) != setIndices.end()) continue;//we've handled this already
        double L_diag = 0.0;//diagonal entries of L
        bool isMappedVertex = false;
        std::pair<int, int> pairing;
        //iterate over the one-ring of the vertex 
        for (Halfedge he : v.outgoingHalfedges()){
            std::cout << "for vertex " << v << " " << originalIndexToLaplacianMatrixIndex.at(v.getIndex()) << std::endl;
            //off diagonal entries 
            tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(he.tipVertex().getIndex()),
            -geometry.edgeCotanWeights[he.edge()]);
            setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
            //search for mappings 
            for (auto p : vertexMappingsPairs){
                if ((p.first == iV) && (std::find(seenPairs.begin(), seenPairs.end(), p) == seenPairs.end())){//this vertex has a mapping and we have not seen this pair before
                    isMappedVertex = true;
                    pairing = p;
                    seenPairs.push_back(p);
                    break;
                }
            }
            if (isMappedVertex){
                Vertex mappedVertex = mesh.vertex(pairing.second);
                //iterate over the 1-ring of the mapped vertex 
                for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                    std::cout << "for vertex " << v << " " << originalIndexToLaplacianMatrixIndex.at(v.getIndex()) << std::endl;
                    tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), 
                                        originalIndexToLaplacianMatrixIndex.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                    setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
                }
            }
            //handle the diagonal entry case
            L_diag += geometry.edgeCotanWeights[he.edge()];
            if (isMappedVertex){//handle diagonal entry for stitched vertex
                Vertex mappedVertex = mesh.vertex(pairing.second);
                //iterate over the 1-ring of the mapped vertex 
                for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                    L_diag += geometry.edgeCotanWeights[mappedVertexHalfedge.edge()];
                }
            }
        }
        std::cout << "for vertex " << v << " " << originalIndexToLaplacianMatrixIndex.at(v.getIndex()) << std::endl;
        tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(v.getIndex()), L_diag);
        setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
    }
    std::cout << "size of triplet list = " << tripletList.size() << std::endl;
    for (auto entry : tripletList){
        std::cout << "row " << entry.row() << " col " << entry.col() << "val " << entry.value() << std::endl;
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
    */
    return timeFunction;

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