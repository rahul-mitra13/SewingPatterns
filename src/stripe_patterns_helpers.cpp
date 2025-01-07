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
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                              HalfedgeData<double>& sigmaTilde, float period){

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
  sigma_mod[startVertex] = 0;
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
  
  if (std::abs(val1 - val2) < pert) return false; // this eliminates halfedges that are exactly parallel to the isoline
  
  //just check if the isovalue lies in between the two values 
  if (currIsoVal > (std::min(val1, val2) - pert) && currIsoVal < (std::max(val1, val2) + pert)){
    bary = (currIsoVal - val2) / (val1 - val2);
    return true;
  }

  return false;

}

//find the set of edges that a given isovalue of the stripes crosses
EdgeData<int> stripeIsoValEdges(EmbeddedGeometryInterface& geometry,
                                  const CornerData<double>& stripeValues,
                                  const FaceData<int>& stripesIndices, double period, polyscope::SurfaceMesh& psMesh){
  
  SurfaceMesh& mesh = geometry.mesh;
  HalfedgeData<std::vector<double>> halfedgeIsoValues = getHalfEdgeIsoValues(geometry, stripeValues, stripesIndices, period);
  //decalare a custom comparator for floating point values 
  auto customComparator = [](const double a, const double b){
    double epsilon = 1e-9;
    return ((a < b) && (std::fabs(a - b) >= epsilon));
  };
  std::set<double, decltype(customComparator)> uniqueIsoVals(customComparator);
  //std::set<double> uniqueIsoVals;
  EdgeData<int> result(mesh, 0);
  
  //find all the faces that a given isovalue passes through
  for (Halfedge h : mesh.halfedges()){
    if (result[h.edge()]) continue;
    if (h.face().isBoundaryLoop() || stripesIndices[h.face()] != 0) continue;//skip singular faces and boundary faces
    std::vector<double> heSet = halfedgeIsoValues[h];
    for (int i = 0; i < heSet.size(); i++){
      result[h.edge()] = 1;
      uniqueIsoVals.insert(heSet[i]);
    }
  }

  //a visited edges map where the key is 
  //(leader edge, isoval) -> set of edges on that isovalue visited from that leader
  std::map<std::pair<int, double>, std::vector<int>> visitedEdges;
  //a used edges map where the key is 
  //(edge index, isoval) -> whether or not that edge has been visited by that isovalue
  std::map<std::pair<int, double>, int> usedEdges;
  double bary;
  EdgeData<int> leaders(mesh, 0);

  //this approach will have some repitions since isolines take on 
  //different values depending on the edges they're on 
  //but if the isoline hits the same edges it shouldn't matter in the long run
  for (double isoVal : uniqueIsoVals){
    std::cout << "unique isoval = " << isoVal << std::endl;
    for (Halfedge he : mesh.halfedges()){
      if (halfedgeContainsLevelSet(stripeValues[he.corner()], stripeValues[he.next().corner()], bary, isoVal) && 
        usedEdges.find(std::make_pair(he.edge().getIndex(), isoVal)) == usedEdges.end()){//this edge has not been visited for isoVal
        std::cout << "leader for isovalue " << isoVal << " is edge " << he.edge().getIndex() << std::endl;
        leaders[he.edge()] = 1;
        usedEdges[std::make_pair(he.edge().getIndex(), isoVal)] = 1;
        visitedEdges[std::make_pair(he.edge().getIndex(), isoVal)] = std::vector<int>{static_cast<int>(he.edge().getIndex())};

        //walk along this isovalue starting from the leader
        Halfedge startHe = he;
        Halfedge currHe = startHe;
        while(true){
          Halfedge jk = currHe.next();
          Halfedge ki = currHe.next().next();
          if (halfedgeContainsLevelSet(stripeValues[jk.corner()], stripeValues[jk.next().corner()], bary, isoVal) && 
          usedEdges.find(std::make_pair(jk.edge().getIndex(), isoVal)) == usedEdges.end()){//this edge has not been visited for this isoval
            
            usedEdges[std::make_pair(jk.edge().getIndex(), isoVal)] = 1;
            visitedEdges[std::make_pair(he.edge().getIndex(), isoVal)].push_back(jk.edge().getIndex());
            currHe = jk;
          }

          else if (halfedgeContainsLevelSet(stripeValues[ki.corner()], stripeValues[ki.next().corner()], bary, isoVal) && 
          usedEdges.find(std::make_pair(ki.edge().getIndex(), isoVal)) == usedEdges.end()){//this edge has not been visited for this isoval
            
            usedEdges[std::make_pair(ki.edge().getIndex(), isoVal)] = 1;
            visitedEdges[std::make_pair(he.edge().getIndex(), isoVal)].push_back(ki.edge().getIndex());
            currHe = ki;
          }

          else{//isoline ends
            std::cout << "iso line is ending! " << std::endl;
            break;
          }

          if (currHe == startHe){//made a closed loop
            std::cout << "isoline made a closed! " << std::endl;
            break;
          }

          //otherwise keep walking
          currHe = currHe.twin();

        }
      }
    }
  }

  std::cout << "size of visited edges = " << visitedEdges.size() << std::endl;

  for (auto entry : visitedEdges){
    std::cout << "for pair " << entry.first.first << ", " << entry.first.second << std::endl;
    for (int i = 0; i < entry.second.size(); i++){
      std::cout << "edge it hits = " << entry.second[i] << std::endl;
    }
    std::cout << "-----------------------------" << std::endl;
  }

  psMesh.addEdgeScalarQuantity("edge leaders", leaders);

  

  return result;
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
        currPoint.he = h;
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

// lists all the k's \in \mathbb{Z} such that val1 < 2 k \pi < val2 or val2 < 2 k \pi < val1
std::vector<double> crossingsModuloPeriod(double val1, double val2, double period) {
  std::vector<double> barys;
  if (val1 == val2) return barys;

  int maxCrossings = std::ceil(std::abs(val1 - val2) / (period));

  if (val1 < val2) {
    for (int i = 0; i < maxCrossings; ++i) {
      int k = std::ceil(val1 / (period)) + i;
      double isoval = 2 * PI * k;

      if (isoval < val2) {
        barys.push_back((isoval - val2) / (val1 - val2));
      }
    }
  } else {
    for (int i = maxCrossings - 1; i >= 0; --i) {
      int k = std::ceil(val2 / (period)) + i;
      double isoval = 2 * PI * k;

      if (isoval < val1) {
        barys.push_back((isoval - val2) / (val1 - val2));
      }
    }
  }

  return barys;
}

// Matches crossings based on a strategy proposed in "Navigating intrinsic triangulations" [Sharp et al. 2019].
// See https://github.com/nmwsharp/geometry-central/pull/89#issuecomment-936150222 for more details
std::vector<std::array<int, 2>> matchCrossings(const std::vector<std::vector<int>>& crossings) {
  assert(crossings.size() == 3);

  int idxIJ = 2;
  if (crossings[0].size() >= crossings[1].size() && crossings[0].size() >= crossings[2].size()) {
    idxIJ = 0;
  } else if (crossings[1].size() >= crossings[2].size() && crossings[1].size() >= crossings[0].size()) {
    idxIJ = 1;
  }
  int idxJK = (idxIJ + 1) % 3;
  int idxKI = (idxIJ + 2) % 3;

  const std::vector<int>& IJ = crossings[idxIJ];
  const std::vector<int>& JK = crossings[idxJK];
  const std::vector<int>& KI = crossings[idxKI];

  int nIJ = IJ.size();
  int nJK = JK.size();
  int nKI = KI.size();

  assert(nIJ >= nJK && nIJ >= nKI);
  assert(nIJ <= nJK + nKI);
  assert((nIJ + nJK + nKI) % 2 == 0);

  std::vector<std::array<int, 2>> matchings;
  if (nIJ == nJK + nKI) { // Case 1: all edges intersecting ijk cross a common edge ij
    // match IJ with IK
    for (int m = 0; m < nKI; ++m) {
      matchings.push_back({IJ[m], KI[nKI - m - 1]});
    }
    // match IJ with KJ
    for (int m = 0; m < nJK; ++m) {
      matchings.push_back({IJ[nKI + m], JK[nJK - m - 1]});
    }
  } else { // Case 2: there is no common edge
    int nRemainingCrossings = (nIJ + nJK + nKI) / 2;
    int m = 0;
    while (nRemainingCrossings > nJK) {
      matchings.push_back({IJ[m], KI[nKI - m - 1]});
      ++m;
      nRemainingCrossings -= 1;
    }

    int l = 0;
    while (nRemainingCrossings > nKI - m) {
      matchings.push_back({IJ[nIJ - 1 - l], JK[l]});
      nRemainingCrossings -= 1;
      ++l;
    }

    int p = 0;
    while (nRemainingCrossings > 0) {
      matchings.push_back({JK[nJK - 1 - p], KI[p]});
      ++p;
      nRemainingCrossings -= 1;
    }
  }

  return matchings;
}


//utility functions to find level sets of stripes
bool existsInVertexList(std::vector<Vector3>& vertices, Vector3& newVertex){

  for (Vector3 v : vertices){
    if (norm(v - newVertex) < 1e-8)
      return true;
  }
  return false;

}
bool existsInVertexList(std::vector<PolyLinePoint>& vertices, PolyLinePoint& p){

  for (PolyLinePoint v : vertices){
    if (norm(v.position - p.position) < 1e-8)
      return true;
  }
  return false;

}

int find(std::vector<Vector3>& vertices, Vector3& newVertex){

  for (int i = 0; i < vertices.size(); i++){
    if (norm(vertices[i] - newVertex) < 1e-8){
      return i;
    }
  }

  return -1;

}

int find(std::vector<PolyLinePoint>& vertices, PolyLinePoint& p){

  for (int i = 0; i < vertices.size(); i++){
    if (norm(vertices[i].position - p.position) < 1e-8){
      return i;
    }
  }

  return -1;

}


std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>> removeCurveNetworkDuplicatedVertices(VertexPositionGeometry& globalGeometry, 
                                                                        std::vector<Vector3>& vertices, std::vector<std::array<int, 2>>& edges){

                                                                                  
  // SurfaceMesh& globalMesh = globalGeometry.mesh;
  // // Map to store unique vertices and their indices
  // std::unordered_map<Vector3, size_t> uniqueVertices;
  // std::vector<size_t> vertexMapping(vertices.size(), -1);
  // std::vector<Vector3> newVertices;
  // double eps = 1e-8;

  // // Identify unique vertices and build the mapping
  // for (size_t i = 0; i < vertices.size(); ++i) {
  //   Vector3 pos = vertices[i];
  //   if (uniqueVertices.find(pos) == uniqueVertices.end()) {
  //       // New unique vertex
  //       uniqueVertices[pos] = newVertices.size();
  //       vertexMapping[i] = newVertices.size();
  //       newVertices.push_back(pos);
  //     } else {
  //       // Duplicate vertex, map to the unique one
  //       vertexMapping[i] = uniqueVertices[pos];
  //     }
  // }

  SurfaceMesh& globalMesh = globalGeometry.mesh;
  std::vector<Vector3> uniqueVertices;
  std::vector<size_t> vertexMapping(vertices.size(), -1);
  std::vector<Vector3> newVertices;
  double eps = 1e-8;

  // Identify unique vertices and build the mapping
  for (size_t i = 0; i < vertices.size(); ++i) {
    Vector3 pos = vertices[i];
    if (!existsInVertexList(uniqueVertices, pos)){
      //New unique vertex
      uniqueVertices.push_back(pos);
      vertexMapping[i] = newVertices.size();
      newVertices.push_back(pos);
    } 
    else {
      // Duplicate vertex, map to the unique one
      //vertexMapping[i] = uniqueVertices[pos];
      vertexMapping[i] = find(uniqueVertices, pos);
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
    std::cout << "number of components = " << components.size() << std::endl;
    // for (const auto& comp : components) {
    //     std::cout << "Component " << componentIndex++ << ": ";
    //     for (size_t v : comp.second) {
    //         std::cout << v << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>> findStripeConnectedComponents(VertexPositionGeometry& globalGeometry, 
                                  EdgeLengthGeometry& gluedGeometry, 
                                  const CornerData<double>& stripeValues,
                                  const FaceData<int>& stripesIndices, double period,
                                  std::map<int, int>& edgeMap,
                                  std::unordered_map<size_t, std::vector<PolyLinePoint>>& components){
  
  SurfaceMesh& globalMesh = globalGeometry.mesh; 
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;
  //clear out the map first 
  components.clear();

  globalGeometry.requireFaceIndices();
  globalGeometry.requireVertexPositions();

  //pertubation
  double pert = 1e-8;

  std::vector<PolyLinePoint> isoLinePoints;
  FaceData<std::vector<int>> isoLinePointsPerFace(globalMesh);
  HalfedgeData<std::vector<int>> isoLinePointsPerHalfedge(globalMesh);

  std::vector<Vector3> isoLinePointCoords;

  // For each halfedge, find all isolines traversing it
  HalfedgeData<std::vector<double>> halfedgeIsoValues = getHalfEdgeIsoValues(globalGeometry, stripeValues, stripesIndices, period);

  //first generate all the points
  for (Halfedge h : globalMesh.halfedges()){

    if (h.face().isBoundaryLoop() || stripesIndices[h.face()] != 0) continue;//skip singular faces and boundary faces
    std::vector<double> heSet = halfedgeIsoValues[h]; // isolines traversing h

    for (int i = 0; i < heSet.size(); i++){
      double bary;
      if (halfedgeContainsLevelSet(stripeValues[h.corner()], stripeValues[h.next().corner()], bary, heSet[i])){//this will always be true, except if the halfedge is parallel to the level set
        //interpolating from tip to tail
        Vector3 pos = (bary * globalGeometry.vertexPositions[h.tailVertex()] +
                       (1 - bary) * globalGeometry.vertexPositions[h.tipVertex()]);
        
        PolyLinePoint currPoint;
        currPoint.position = pos;
        currPoint.f = h.face();
        currPoint.e = h.edge();
        currPoint.he = h;
        currPoint.isoval = heSet[i];

        isoLinePointsPerFace[h.face()].push_back(isoLinePoints.size());
        isoLinePointsPerHalfedge[h].push_back(isoLinePoints.size());
        isoLinePoints.push_back(currPoint);
        isoLinePointCoords.push_back(pos);

      }
    }
  }

  // polyscope::registerPointCloud("isoLinePoints", isoLinePointCoords);

  std::vector<std::array<int, 2>> edges;
  //now make the polyline 
  for (Face f : globalMesh.faces()) {
    std::vector<int> facePoints = isoLinePointsPerFace[f];
    for (int i : facePoints) {
      PolyLinePoint p1 = isoLinePoints[i];
      //search for the point p1 should connect to
      for (int j : facePoints) {
        PolyLinePoint p2 = isoLinePoints[j];
        if (p1.f == p2.f && std::abs(p1.isoval - p2.isoval) < pert && i != j){
          std::array<int, 2> edge = {i, j};
          edges.push_back(edge);
        }
      }
    }
  }

  // Build mapping for vertices at same location, using a union find data structure
  UF clusters(isoLinePoints.size());
  for (Halfedge he : globalMesh.halfedges()) {
    for (int i : isoLinePointsPerHalfedge[he]) { // for each point on halfedge
      PolyLinePoint p = isoLinePoints[i];
      for (int j : isoLinePointsPerHalfedge[he.twin()]) {
        PolyLinePoint q = isoLinePoints[j];
        if (norm(p.position - q.position) < 1e-8) // epsilon choice seems OK
          clusters.merge(i, j);
      }
      // Also search on the face. This is useful when an isoline passes exactly through a mesh vertex
      for (int j : isoLinePointsPerFace[he.face()]) {
        PolyLinePoint q = isoLinePoints[j];
        if (norm(p.position - q.position) < 1e-8) // epsilon choice seems OK
          clusters.merge(i, j);
      }
    }
  }

  //now remove all the duplicate points 
  std::vector<PolyLinePoint> uniqueVertices;
  std::vector<int> vertexMapping(isoLinePoints.size(), -1);
  std::vector<PolyLinePoint> newVertices; // what's the difference with uniqueVertices?
  
  // Find leaders
  for (int i = 0; i < isoLinePoints.size(); i++) {
    if (clusters.find(i) == i) { // if this is a leader
      vertexMapping[i] = uniqueVertices.size();
      uniqueVertices.push_back(isoLinePoints[i]);
      newVertices.push_back(isoLinePoints[i]);
    }
  }
  // Assign the remaining mappings
  for (int i = 0; i < isoLinePoints.size(); i++)
    vertexMapping[i] = vertexMapping[clusters.find(i)];

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

  //finally find the connected components
  // Initialize disjoint set for vertices
  DisjointSets vertexSets(newVertices.size());

  // Union-find: process each edge to connect the vertices
  for (auto e : newEdges) {
    auto v1 = e[0];
    auto v2 = e[1];
    vertexSets.merge(v1, v2);
  }

  // Map each set leader to its component
  std::unordered_map<size_t, std::vector<size_t>> indexedComponents;

  std::vector<Vector3> newVertexPositions;
  for (size_t v = 0; v < newVertices.size(); ++v) {
    newVertexPositions.push_back(newVertices[v].position);
    size_t leader = vertexSets.find(v);
    indexedComponents[leader].push_back(v);
  }

  //reset the component index
  size_t componentIndex = 0;
  // std::cout << "number of components = " << indexedComponents.size() << std::endl;
  for (const auto& comp : indexedComponents) {
      //std::cout << "Component " << componentIndex++ << ": ";
      componentIndex++;
      for (size_t v : comp.second) {
      //    std::cout << v << " ";
          components[componentIndex].push_back(newVertices[v]);
      }
      //std::cout << std::endl;
  }

  return std::tie(newVertexPositions, newEdges);

}

//find a set of edge singularity pairs on the same isoline as the stripe patterns
std::vector<std::pair<int, int>> findEdgeSingularityPairFromStripeIsoVals(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                                            EdgeData<double>& edgeCurl, std::map<int, int>& edgeMap,
                                                            std::unordered_map<size_t, std::vector<PolyLinePoint>>& components, int numPairs){
  
  //causing some malloc errors
  //edge curl is in the global setting (maybe should do it in the glued setting)
  //return top 3 pairs here?
  SurfaceMesh& globalMesh = globalGeometry.mesh;
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;
  double eps = 1e-8;

  //a map from edge curl to component id
  std::vector<double> componentCurl(components.size());
  std::map<int, size_t> edgeCurlToComponentId;

  for (const auto &c : components){
    std::vector<PolyLinePoint> vertices = c.second;
    double curlSum = 0.0;
    for (const auto &v : vertices){
      curlSum += std::fabs(edgeCurl[v.e]);
    }
    double avgCurl = curlSum / vertices.size();
    edgeCurlToComponentId[hashFloatQuantized(curlSum)] = c.first;
    componentCurl[c.first - 1] = avgCurl;
    //componentCurl[c.first - 1] = curlSum;
  }

  std::sort(componentCurl.begin(), componentCurl.end(), std::greater<double>()); // Sort in descending order
  componentCurl.resize(numPairs); // Keep only the top numPairs

  std::vector<std::pair<int, int>> edgeSingularityPairs;
  std::cout << "numPairs = " << numPairs << std::endl;
  std::cout << "size of component curl = " << componentCurl.size() << std::endl;
  std::cout << "number of components = " << components.size() << std::endl;
  std::cout << "size of edge curl to component Id = " << edgeCurlToComponentId.size() << std::endl;

  //now find the max and min curl edges on those components
  for (const auto &c : componentCurl){
    auto vertices = components[edgeCurlToComponentId[hashFloatQuantized(c)]];
    double maxCurl = -DBL_MAX;
    double minCurl = DBL_MAX;
    int maxEdge, minEdge = -1;
    for (const auto &v : vertices){
      if (edgeCurl[v.e] > maxCurl){
        maxCurl = edgeCurl[v.e];
        maxEdge = v.e.getIndex();
      }
      if (edgeCurl[v.e] < minCurl){
        minCurl = edgeCurl[v.e];
        minEdge = v.e.getIndex();
      }
    }
    edgeSingularityPairs.push_back(std::make_pair(maxEdge, minEdge));
  }
  return edgeSingularityPairs; 
}


//given a stripe pattern as corner data, compute an integer value around a specified path
std::vector<int> computeIntegralValueAlongPaths(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                CornerData<double>& stripeValues, std::vector<std::vector<double>>& paths, double period){
  
  SurfaceMesh& globalMesh = globalGeometry.mesh;
  SurfaceMesh& gluedMesh = gluedGeometry.mesh;

  std::vector<int> integerVariables(paths.size());
  std::cout << "size of paths = " << paths.size() << std::endl;

  for (int i = 0; i < paths.size(); i++){
    std::vector<double> path = paths[i];
    double sum = 0;
    int integerVal = 0;
    double pathIntegral = 0;
  
    for (int j = 0; j < gluedMesh.nEdges(); j++){
      if (std::fabs(path[j]) > 1e-10){
        sum += stripeValues[gluedMesh.edge(j).halfedge().corner()];
      }
    }
    std::cout << "period = " << period << std::endl;
    std::cout << "sum for boundary " << i << " = " << sum << std::endl;
    integerVal = std::round(sum  / period);
    std::cout << "integer val = " << integerVal << std::endl;
    integerVariables[i] = integerVal;
  }

  return integerVariables;

}

//-------------------------re-implementing some of  Knoppel's stripes-----------------------//

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