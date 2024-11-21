#include "stripe_patterns_helpers.h"

std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& geometry, EdgeData<double>& sigma, float period, SurfaceMesh& gluedMesh,
                                                                            std::vector<std::pair<int, int>>& vertexMappingsPairs, std::vector<std::pair<int, int>>& edgeMappingsPairs,
                                                                            std::map<int, std::vector<Halfedge>>& gluedOneRingMap){
  
  SurfaceMesh& mesh = geometry.mesh;

  //store the sigma values per halfedge
  HalfedgeData<double> sigma_halfedges(mesh);
  
  //assign halfedge sigma values taking orientation into account
  for (Halfedge he : mesh.halfedges()){
  	if (he.orientation()){
  		sigma_halfedges[he] = sigma[he.edge()];
  	}
  	else{
  		sigma_halfedges[he] = -1.0 * sigma[he.edge()];
  	}
  }

  //update the 1-form for edges that are "stitched together"
  //update the stitched together edges to account for correct orientation
  for (std::pair<int, int> p : edgeMappingsPairs){
    Edge e2 = mesh.edge(p.second);
    for (Halfedge he : mesh.halfedges()){
      if (he.edge().getIndex() == e2.getIndex()){
        sigma_halfedges[he] = -1.0 * sigma_halfedges[he];
      }
    }
  }

  //integrate up the sigma values now onto vertices
  VertexData<double> sigma_mod(mesh); 
  VertexData<bool> visited(mesh, false);
  //vertex to start integration
  Vertex startVertex; 
  for (Vertex v : mesh.vertices()){
    //start integrating from a boundary vertex
    if (v.isBoundary()){
      for (std::pair<int, int> p : vertexMappingsPairs){//ensure you're not starting at a stitched vertex (probably don't need this)
        if (v.getIndex() != p.first && v.getIndex() != p.second){
          startVertex = v;
          break;
        }
      }
    }
  }
  //start at a boundary vertex if you want stripes to align perfectly with the boundary
  //perform BFS on the global connected mesh
  startVertex = mesh.vertex(0);
  std::queue<Vertex> Q;
  Q.push(startVertex);
  sigma_mod[startVertex] = 0.0;
  visited[startVertex] = true;
  while(!Q.empty()){
    Vertex vi = Q.front(); Q.pop();
    for (Halfedge he : gluedOneRingMap[vi.getIndex()]){
      Vertex vj = he.twin().vertex();
      if (!visited[vj]){
        sigma_mod[vj] = fmod(sigma_mod[vi] + sigma_halfedges[he], period);
        visited[vj] = true;
        Q.push(vj);
      }
    }
  }

  //assign texture coordinates
  CornerData<double> textureCoordinates(mesh);
  FaceData<int> paramIndices(mesh);

  for (Face f : mesh.faces()){

    //skip boundary faces
    if (f.isBoundaryLoop()) continue;
    
    //grab the halfedges
    Halfedge hij = f.halfedge();
    Halfedge hjk = hij.next();
    Halfedge hki = hjk.next();

    //grab the sigmas
    double sigma_ij = sigma_halfedges[hij];
    double sigma_jk = sigma_halfedges[hjk];
    double sigma_ki = sigma_halfedges[hki];

    //compute alpha values at triangle corners
    double alphaI = sigma_mod[f.halfedge().vertex()];
    double alphaJ = alphaI + sigma_ij;
    double alphaK = alphaJ + sigma_jk;
    double alphaL = alphaK + sigma_ki;

    // store the coordinates
    textureCoordinates[hij.corner()] = alphaI;
    textureCoordinates[hjk.corner()] = alphaJ;
    textureCoordinates[hki.corner()] = alphaK;
    paramIndices[f] = std::round((alphaL - alphaI) / (period));

  }

  return std::tie(textureCoordinates, paramIndices);
}

//carry out the integration in the glued mesh setting
//sigma is defined over the EDGES of glued mesh 
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& globalGeometry, 
                                                              EdgeLengthGeometry& gluedGeometry, EdgeData<double>& sigma, float period){

  SurfaceMesh& globalMesh = globalGeometry.mesh; 
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  //store the value per halfedge 
  HalfedgeData<double> sigma_halfedges(gluedMesh);

  //assign halfedge sigma values taking orientation into account 
  for (Halfedge he : gluedMesh.halfedges()){
  	if (he.orientation()){
  		sigma_halfedges[he] = sigma[he.edge()];
  	}
  	else{
  		sigma_halfedges[he] = -1.0 * sigma[he.edge()];
  	}
  }

  //integrate up the sigma values now onto vertices
  VertexData<double> sigma_mod(gluedMesh); 
  VertexData<bool> visited(gluedMesh, false);
  //vertex to start integration
  Vertex startVertex; 
  for (Vertex v : gluedMesh.vertices()){
    //start integrating from a boundary vertex
    if (v.isBoundary()){
      startVertex = v;
      break;
    }
  }

  std::queue<Vertex> Q;
  Q.push(startVertex);
  // Floating point thing for knit graph
  sigma_mod[startVertex] = 0.0 + 1e-16;
  visited[startVertex] = true;
  while(!Q.empty()){
    Vertex vi = Q.front(); Q.pop();
    Halfedge h = vi.halfedge();
    do{
      Vertex vj = h.twin().vertex();
      if (!visited[vj]){
        sigma_mod[vj] = fmod(sigma_mod[vi] + sigma_halfedges[h], period);
        visited[vj] = true;
        Q.push(vj);
      }
      h = h.twin().next();
    }
    while( h != vi.halfedge());
  }

  //assign texture coordinates to the glued mesh 
  CornerData<double> textureCoordinates(gluedMesh);
  FaceData<int> paramIndices(gluedMesh);
  
  //assign texture coordinates to the global mesh 
  CornerData<double> textureCoordinatesGlobal(globalMesh);
  FaceData<int> paramIndicesGlobal(globalMesh);

  for (Face f : gluedMesh.faces()){

    //skip boundary faces
    if (f.isBoundaryLoop()) continue;
    
    //grab the halfedges
    Halfedge hij = f.halfedge();
    Halfedge hjk = hij.next();
    Halfedge hki = hjk.next();

    //grab the sigmas
    double sigma_ij = sigma_halfedges[hij];
    double sigma_jk = sigma_halfedges[hjk];
    double sigma_ki = sigma_halfedges[hki];

    //compute alpha values at triangle corners
    double alphaI = sigma_mod[f.halfedge().vertex()];
    double alphaJ = alphaI + sigma_ij;
    double alphaK = alphaJ + sigma_jk;
    double alphaL = alphaK + sigma_ki;
    // store the coordinates in the glued mesh 
    textureCoordinates[hij.corner()] = alphaI;
    textureCoordinates[hjk.corner()] = alphaJ;
    textureCoordinates[hki.corner()] = alphaK;
    paramIndices[f] = std::round((alphaL - alphaI) / (period));

    //store the coordiantes in the global mesh 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().corner()] = alphaI; 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().next().corner()] = alphaJ; 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().next().next().corner()] = alphaK; 
    paramIndicesGlobal[globalMesh.face(f.getIndex())] = std::round((alphaL - alphaI) / (period));
  }

  return std::tie(textureCoordinatesGlobal, paramIndicesGlobal);
}

//carry out the integration in the glued mesh setting
//sigma is defined over the HALFEDGES of the glued mesh
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde, float period){

  SurfaceMesh& globalMesh = globalGeometry.mesh; 
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;
  
  //integrate up the sigma values now onto vertices
  VertexData<double> sigma_mod(gluedMesh); 
  VertexData<bool> visited(gluedMesh, false);
  //vertex to start integration
  Vertex startVertex; 
  for (Vertex v : gluedMesh.vertices()){
    //start integrating from a boundary vertex
    if (v.isBoundary()){
      startVertex = v;
      break;
    }
  }

  std::queue<Vertex> Q;
  Q.push(startVertex);
  // Floating point thing for knit graph
  sigma_mod[startVertex] = 0.0 + 1e-16;
  visited[startVertex] = true;
  while(!Q.empty()){
    Vertex vi = Q.front(); Q.pop();
    Halfedge h = vi.halfedge();
    do{
      Vertex vj = h.twin().vertex();
      if (!visited[vj]){
        sigma_mod[vj] = fmod(sigma_mod[vi] + sigmaTilde[h], period);
        visited[vj] = true;
        Q.push(vj);
      }
      h = h.twin().next();
    }
    while( h != vi.halfedge());
  }

  //assign texture coordinates to the glued mesh 
  CornerData<double> textureCoordinates(gluedMesh);
  FaceData<int> paramIndices(gluedMesh);
  
  //assign texture coordinates to the global mesh 
  CornerData<double> textureCoordinatesGlobal(globalMesh);
  FaceData<int> paramIndicesGlobal(globalMesh);

  for (Face f : gluedMesh.faces()){

    //skip boundary faces
    if (f.isBoundaryLoop()) continue;
    
    //grab the halfedges
    Halfedge hij = f.halfedge();
    Halfedge hjk = hij.next();
    Halfedge hki = hjk.next();

    //grab the sigmas
    double sigma_ij = sigmaTilde[hij];
    double sigma_jk = sigmaTilde[hjk];
    double sigma_ki = sigmaTilde[hki];

    //compute alpha values at triangle corners
    double alphaI = sigma_mod[f.halfedge().vertex()];
    double alphaJ = alphaI + sigma_ij;
    double alphaK = alphaJ + sigma_jk;
    double alphaL = alphaK + sigma_ki;
    // store the coordinates in the glued mesh 
    textureCoordinates[hij.corner()] = alphaI;
    textureCoordinates[hjk.corner()] = alphaJ;
    textureCoordinates[hki.corner()] = alphaK;
    paramIndices[f] = std::round((alphaL - alphaI) / (period));

    //store the coordiantes in the global mesh 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().corner()] = alphaI; 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().next().corner()] = alphaJ; 
    textureCoordinatesGlobal[globalMesh.face(f.getIndex()).halfedge().next().next().corner()] = alphaK; 
    paramIndicesGlobal[globalMesh.face(f.getIndex())] = std::round((alphaL - alphaI) / (period));
  }

  return std::tie(textureCoordinatesGlobal, paramIndicesGlobal);
}

//returns a vector per face that stores all the isovalues (modulo period) passing through this face
FaceData<std::vector<double>> getFaceIsoValues(IntrinsicGeometryInterface& geometry,
                                              const CornerData<double>& stripeValues,
                                              const FaceData<int>& stripesIndices, double period){
  
  SurfaceMesh& mesh = geometry.mesh;

  FaceData<std::vector<double>> toReturn(mesh);

  //pertubation to avoid floating point error
  double pert = 0.0000001;
  
  for (Face f : mesh.faces()){
    
    std::vector<double> isovalues;
    
    //skip singular faces 
    if (stripesIndices[f] != 0) continue;

    Halfedge hij = f.halfedge();
    Halfedge hjk = hij.next();
    Halfedge hki = hjk.next();

    //grab the stripe values
    double alphaI = stripeValues[hij.corner()];
    double alphaJ = stripeValues[hjk.corner()];
    double alphaK = stripeValues[hki.corner()];

    double min_alpha_val = std::min({alphaI, alphaJ, alphaK});
    double max_alpha_val = std::max({alphaI, alphaJ, alphaK});
    
    min_alpha_val -= pert;
    max_alpha_val += pert;

    double start = std::ceil(min_alpha_val/period) * period;
    double end = std::floor(max_alpha_val/period) * period;
    double currIsoVal = start;

    while (currIsoVal < end + pert){
      isovalues.push_back(currIsoVal);
      currIsoVal += period;
    }

    toReturn[f] = isovalues;
  }
  return toReturn;
}

//returns a vector per halfedge that stores all the isovalues (modulo period) passing through this halfedge
HalfedgeData<std::vector<double>> getHalfEdgeIsoValues(IntrinsicGeometryInterface& geometry,
                                              const CornerData<double>& stripeValues,
                                              const FaceData<int>& stripesIndices, double period){
  
  geometry.requireFaceIndices();
  
  SurfaceMesh& mesh = geometry.mesh;

  HalfedgeData<std::vector<double>> toReturn(mesh);

  //pertubation to avoid floating point error
  double pert = 0.00001;
  
  for (Halfedge he : mesh.halfedges()){
    
    std::vector<double> isovalues;
    double tipValue;
    double tailValue;

    tailValue = stripeValues[he.corner()];
    tipValue = stripeValues[he.next().corner()];

    double min_val = std::min({tailValue, tipValue});
    double max_val = std::max({tailValue, tipValue});

    min_val -= pert;
    max_val += pert;

    double start = std::ceil(min_val/period) * period;
    double end = std::floor(max_val/period) * period;
    double currIsoVal = start;

    while (currIsoVal < end + pert){
      isovalues.push_back(currIsoVal);
      currIsoVal += period;
    }

    toReturn[he] = isovalues;
  }

  return toReturn;
}

//if the halfedge contains the passed iso value
//returns the ratio along the halfedge that the intersection exists
//otherwise return false
bool halfedgeContainsLevelSet(double val1, double val2, double& bary, double currIsoVal){
  
  //pertubation
  double pert = 0.00001;
  
  if (std::abs(val1 - val2) < pert) return false;
  
  //just check if the isovalue lies in between the two values 
  if (currIsoVal > (std::min(val1, val2) - pert) && currIsoVal < (std::max(val1, val2) + pert)){
    bary = (currIsoVal - val2) / (val1 - val2);
    return true;
  }

  return false;

}

//extract isolines when multiple isolines may pass through a face (using face-based matching)
std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>> generateIsoLines(EmbeddedGeometryInterface& geometry,
                                                      const CornerData<double>& stripeValues,
                                                      const FaceData<int>& stripesIndices, double period){

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireFaceIndices();
  geometry.requireVertexPositions();

  //pertubation
  double pert = 0.000001;

  std::vector<PolyLinePoint> isoline_points;

  HalfedgeData<std::vector<double>> halfedgeIsoValues = getHalfEdgeIsoValues(geometry, stripeValues, stripesIndices, period);

  //first generate all the points
  for (Halfedge h : mesh.halfedges()){

    if (h.face().isBoundaryLoop() || stripesIndices[h.face()] != 0) continue;//skip singular faces and boundary faces
    std::vector<double> heSet = halfedgeIsoValues[h];

    for (int i = 0; i < heSet.size(); i++){
      double bary;
      if (halfedgeContainsLevelSet(stripeValues[h.corner()], stripeValues[h.next().corner()], bary, heSet[i])){//this will always be true
        //interpolating from tip to tail
        Vector3 pos = (bary * geometry.vertexPositions[h.tailVertex()] +
                       (1 - bary) * geometry.vertexPositions[h.tipVertex()]);
        
        PolyLinePoint currPoint;
        currPoint.position = pos;
        currPoint.f = h.face();
        currPoint.e = h.edge();
        currPoint.isoval = heSet[i];

        isoline_points.push_back(currPoint);
      }
    }
  }

  std::vector<Vector3> points;
  std::vector<std::array<int, 2>> edges;

  for (int i = 0; i < isoline_points.size(); i++){
    points.push_back(isoline_points[i].position);
  }

  //now make the polyline 
  for (int i = 0; i < isoline_points.size(); i++){
    
    PolyLinePoint p1 = isoline_points[i];

    //search for the point p1 should connect to
    for (int j = 0; j < isoline_points.size(); j++){

      PolyLinePoint p2 = isoline_points[j];
      if (p1.f == p2.f && std::abs(p1.isoval - p2.isoval) < pert && i != j){
        std::array<int, 2> edge = {i, j};
        edges.push_back(edge);
      }
    }
  }

  return std::tie(points, edges);
}

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>> removeCurveNetworkDuplicatedVertices(VertexPositionGeometry& globalGeometry, 
                                                                                                          std::vector<Vector3>& vertices, std::vector<std::array<int, 2>>& edges){


  SurfaceMesh& globalMesh = globalGeometry.mesh;
  // Map to store unique vertices and their indices
  std::unordered_map<Vector3, size_t> uniqueVertices;
  std::vector<size_t> vertexMapping(vertices.size(), -1);
  std::vector<Vector3> newVertices;

  // Identify unique vertices and build the mapping
  for (size_t i = 0; i < vertices.size(); ++i) {
      Vector3 pos = vertices[i];

      if (uniqueVertices.find(pos) == uniqueVertices.end()) {
          // New unique vertex
          uniqueVertices[pos] = newVertices.size();
          vertexMapping[i] = newVertices.size();
          newVertices.push_back(pos);
        } else {
          // Duplicate vertex, map to the unique one
          vertexMapping[i] = uniqueVertices[pos];
        }
    }

    // Create new edges with updated vertex indices
    std::vector<std::array<int, 2>> newEdges;
    for (auto e : edges) {
        int v1 = vertexMapping[e[0]];
        int v2 = vertexMapping[e[1]];

        // Avoid duplicate edges
        if (v1 != v2) {
            if (v1 > v2) std::swap(v1, v2); // Ensure consistent ordering
            newEdges.emplace_back(std::array<int, 2>{v1, v2});
        }
    }

    // Remove duplicate edges
    std::sort(newEdges.begin(), newEdges.end());
    newEdges.erase(std::unique(newEdges.begin(), newEdges.end()), newEdges.end());

    return std::tie(newVertices, newEdges);

}

//find connected components of a curve network
void findCurveNetworkConnectedComponents(VertexPositionGeometry& globalGeometry,
                            std::vector<Vector3>& vertices, std::vector<std::array<int, 2>>& edges){
  
  SurfaceMesh& globalMesh = globalGeometry.mesh;

  // Initialize disjoint set for vertices
  DisjointSets vertexSets(vertices.size());

  // Union-find: process each edge to connect the vertices
  for (auto e : edges) {
      auto v1 = e[0];
      auto v2 = e[1];
      vertexSets.merge(v1, v2);
  }

    // Map each set leader to its component
    std::unordered_map<size_t, std::vector<size_t>> components;
    for (size_t v = 0; v < vertices.size(); ++v) {
        size_t leader = vertexSets.find(v);
        components[leader].push_back(v);
    }

    // Print the connected components
    size_t componentIndex = 0;
    for (const auto& comp : components) {
        std::cout << "Component " << componentIndex++ << ": ";
        for (size_t v : comp.second) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

}

//-------------------------re-impleminting Knoppel's stripes-----------------------//

// Compute the 1-form \omega_{ij} such as defined in eq.7 of [Knoppel et al. 2015]
//this returns omega per edge in the glued mesh setting 
EdgeData<double> computeOmega(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                    std::map<int, int>& globalToGluedEdgeMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, int direction, FaceData<Vector3>& gradient){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    globalGeometry.requireFaceNormals();
    EdgeData<double> globalOmega(globalMesh);
    EdgeData<double> gluedOmega(gluedMesh);
    FaceData<Vector3> faceGradientsCopy = gradient;
    //if we're computing the matching 1-form in the wale direction, rotate all the gradients 
    if (direction == 1){
        for (Face f : globalMesh.faces()){
          faceGradientsCopy[f] = faceGradientsCopy[f].rotateAround(globalGeometry.faceNormals[f], PI/2.);
        }
    }

    //create a map from the mapped edges
    std::map<int, int> edgeMap;
    for (std::pair<int, int> pair : edgeMappingsPairs){
        edgeMap.insert({pair.first, pair.second});
    }

    //edges which we've handles already
    std::map<int, bool> seenEdges;
    for (Edge e : globalMesh.edges()){
        seenEdges.insert({e.getIndex(), false});
    }
    for (Edge e : globalMesh.edges()){
        if (seenEdges[e.getIndex()]) continue;
        if (e.halfedge().twin().isInterior()){//found an interior halfedge
            Vector3 eVector = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
            globalOmega[e] = 0.5 * dot((faceGradientsCopy[e.halfedge().face()] + faceGradientsCopy[e.halfedge().twin().face()]),
                                    eVector);
            seenEdges[e.getIndex()] = true;
        }
        else{//found a boundary halfedge
            if (edgeMap.find(e.getIndex()) != edgeMap.end()){//found a stitched together edge
                Vector3 faceGradientsCopy1 = faceGradientsCopy[e.halfedge().face()];
                Vector3 faceGradientsCopy2 = faceGradientsCopy[globalMesh.edge(edgeMap.at(e.getIndex())).halfedge().face()];
                Vector3 e1 = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
                globalOmega[e] = 0.5 * dot((faceGradientsCopy1 + faceGradientsCopy2), e1);
                globalOmega[globalMesh.edge(edgeMap.at(e.getIndex()))] = 0.5 * dot((faceGradientsCopy1 + (faceGradientsCopy2)), e1);
                seenEdges[e.getIndex()] = true;
                seenEdges[edgeMap.at(e.getIndex())] = true;
            }
            else{//found a boundary edge that's not stitched to anything
                faceGradientsCopy[e.halfedge().face()] = faceGradientsCopy[e.halfedge().face()];
                Vector3 eVector = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
                globalOmega[e] = dot(faceGradientsCopy[e.halfedge().face()], eVector);
                seenEdges[e.getIndex()] = true;
            }
        }
    }
    
    gluedOmega = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, globalOmega, globalToGluedEdgeMap);

    return gluedOmega;

}

// Build a Laplace-like matrix with double entries (necessary to represent complex conjugation)
SparseMatrix<double> buildVertexEnergyMatrix(EdgeLengthGeometry& gluedGeometry, const FaceData<int>& branchIndices, const EdgeData<double>& omega){

  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  gluedGeometry.requireVertexIndices();
  gluedGeometry.requireHalfedgeCotanWeights();

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(12 * gluedMesh.nEdges());
  for (Edge e : gluedMesh.edges()) {
    // compute the discrete 1-form
    double omegaIJ = omega[e];

    // compute the cotan weight
    double w = 0;
    if (branchIndices[e.halfedge().face()] == 0) {
      w += gluedGeometry.halfedgeCotanWeights[e.halfedge()];
    }
    if (!e.isBoundary() && branchIndices[e.halfedge().twin().face()] == 0) {
      w += gluedGeometry.halfedgeCotanWeights[e.halfedge().twin()];
    }

    int i = 2 * gluedGeometry.vertexIndices[e.halfedge().vertex()];
    int j = 2 * gluedGeometry.vertexIndices[e.halfedge().twin().vertex()];

    // add the diagonal terms
    triplets.emplace_back(i + 0, i + 0, w);
    triplets.emplace_back(i + 1, i + 1, w);

    triplets.emplace_back(j + 0, j + 0, w);
    triplets.emplace_back(j + 1, j + 1, w);

    // compute the new transport coefficient
    Vector2 rij = w * Vector2::fromAngle(omegaIJ);

    // these terms are the same in both cases
    triplets.emplace_back(i + 0, j + 0, -rij.x);
    triplets.emplace_back(i + 1, j + 0, rij.y);

    triplets.emplace_back(j + 0, i + 0, -rij.x);
    triplets.emplace_back(j + 0, i + 1, rij.y);

    triplets.emplace_back(i + 0, j + 1, -rij.y);
    triplets.emplace_back(i + 1, j + 1, -rij.x);

    triplets.emplace_back(j + 1, i + 0, -rij.y);
    triplets.emplace_back(j + 1, i + 1, -rij.x);

  }

  // assemble matrix from triplets
  SparseMatrix<double> vertexEnergyMatrix(2 * gluedMesh.nVertices(), 2 * gluedMesh.nVertices());
  vertexEnergyMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // Shift to avoid singularity
  SparseMatrix<double> eye(2 * gluedMesh.nVertices(), 2 * gluedMesh.nVertices());
  eye.setIdentity();
  vertexEnergyMatrix += 1e-4 * eye;
  
  return vertexEnergyMatrix;
}

// Build a lumped mass matrix with double entries
SparseMatrix<double> computeRealVertexMassMatrix(IntrinsicGeometryInterface& geometry) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireVertexDualAreas();

  std::vector<Eigen::Triplet<double>> triplets(2 * mesh.nVertices());
  for (size_t i = 0; i < mesh.nVertices(); ++i) {
    double area = geometry.vertexDualAreas[i];
    triplets[2 * i] = Eigen::Triplet<double>(2 * i, 2 * i, area);
    triplets[2 * i + 1] = Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area);
  }

  // assemble matrix from triplets
  SparseMatrix<double> realVertexMassMatrix(2 * mesh.nVertices(), 2 * mesh.nVertices());
  realVertexMassMatrix.setFromTriplets(triplets.begin(), triplets.end());

  return realVertexMassMatrix;

}

// Build a lumped mass matrix with double entries
SparseMatrix<double> computeRealVertexMassMatrix(EdgeLengthGeometry& gluedGeometry){

  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  gluedGeometry.requireVertexDualAreas();

  std::vector<Eigen::Triplet<double>> triplets(2 * gluedMesh.nVertices());
  for (size_t i = 0; i < gluedMesh.nVertices(); ++i) {
    double area = gluedGeometry.vertexDualAreas[i];
    triplets[2 * i] = Eigen::Triplet<double>(2 * i, 2 * i, area);
    triplets[2 * i + 1] = Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area);
  }

  // assemble matrix from triplets
  SparseMatrix<double> realVertexMassMatrix(2 * gluedMesh.nVertices(), 2 * gluedMesh.nVertices());
  realVertexMassMatrix.setFromTriplets(triplets.begin(), triplets.end());

  return realVertexMassMatrix;

}

// Solve the generalized eigenvalue problem in equation 9 [Knoppel et al. 2015]
VertexData<Vector2> computeParameterization(EdgeLengthGeometry& gluedGeometry,
                                            const FaceData<int>& branchIndices, const EdgeData<double>& omega){

  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  // Compute vertex energy matrix A and mass matrix B
  SparseMatrix<double> energyMatrix = buildVertexEnergyMatrix(gluedGeometry, branchIndices, omega);
  SparseMatrix<double> massMatrix = computeRealVertexMassMatrix(gluedGeometry);


  // Find the smallest eigenvector
  Vector<double> solution = smallestEigenvectorPositiveDefinite(energyMatrix, massMatrix, 20);

  // Copy the result to a VertexData vector
  VertexData<Vector2> toReturn(gluedMesh);
  for (size_t i = 0; i < gluedMesh.nVertices(); ++i) {
    toReturn[i].x = solution(2 * i);
    toReturn[i].y = solution(2 * i + 1);
    toReturn[i] = toReturn[i].normalize();
  }
  return toReturn;
}

// extract the final texture coordinates from the parameterization
std::tuple<CornerData<double>, FaceData<int>> computeTextureCoordinates(EdgeLengthGeometry& gluedGeometry,
                                                                        const EdgeData<double>& omega,
                                                                        const VertexData<Vector2>& parameterization){


  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  CornerData<double> textureCoordinates(gluedMesh);
  FaceData<int> paramIndices(gluedMesh);

  gluedGeometry.requireTransportVectorsAlongHalfedge();

  for (Face f : gluedMesh.faces()) {
    textureCoordinates[f.halfedge().corner()] = parameterization[f.halfedge().vertex()].arg();

    for (Halfedge he : f.adjacentHalfedges()) {
      // grab the parameter values at vertices
      Vector2 psiI = parameterization[he.vertex()];
      Vector2 psiJ = parameterization[he.next().vertex()];

      // is each halfedge canonical?
      double cIJ = (he.edge().halfedge() != he ? -1 : 1);

      // grab the connection coeffients
      double omegaIJ = cIJ * omega[he.edge()];

      // construct complex transport coefficients
      Vector2 rij = Vector2::fromAngle(omegaIJ);

      if (he.next() != f.halfedge()) {
        textureCoordinates[he.next().corner()] = textureCoordinates[he.corner()] + omegaIJ - (rij * psiI / psiJ).arg();
      } else {
        double alpha = textureCoordinates[he.corner()] + omegaIJ - (rij * psiI / psiJ).arg();
        paramIndices[f] = std::round((alpha - textureCoordinates[he.next().corner()]) / (2 * PI));
      }
    }
  }

  return std::tie(textureCoordinates, paramIndices);
}

// Isolines of this function are stripes perpendicular to the direction field spaced according to the target frequencies
std::tuple<CornerData<double>, FaceData<int>>
computeStripePattern(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                    std::map<int, int>& globalToGluedEdgeMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, int direction, FaceData<Vector3>& gradient, 
                    polyscope::SurfaceMesh& psMesh, std::vector<bool>& orientations){

  SurfaceMesh& globalMesh = globalGeometry.mesh;
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;
  EdgeData<double> omega(gluedMesh);
  EdgeData<double> omegaViz(globalMesh);

  //unnormalized gradient vector field integrable everywhere
  FaceData<int> branchIndices(gluedMesh, 0);

  omega = computeOmega(globalGeometry, gluedGeometry,
                    globalToGluedEdgeMap, edgeMappingsPairs, 1, gradient); 
  
  omegaViz = convertGluedToGlobalEdgeFunction(globalGeometry, gluedGeometry, omega, globalToGluedEdgeMap);
  
  psMesh.addOneFormTangentVectorQuantity("rotated 1-form", omegaViz, orientations);
 
  // find singularities of the direction field
  //FaceData<int> branchIndices = computeFaceIndex(gluedGeometry, directionField, 2);

  // solve the eigenvalue problem (multiply by 2pi to get the right frequencies)
  VertexData<Vector2> parameterization =
      computeParameterization(gluedGeometry, branchIndices, omega);
  
  // compute the final corner-based values, along with singularities of the stripe pattern
  CornerData<double> textureCoordinates;
  FaceData<int> zeroIndices;
  std::tie(textureCoordinates, zeroIndices) =
      computeTextureCoordinates(gluedGeometry, omega, parameterization);

  return std::tie(textureCoordinates, zeroIndices);

}