#include "stripe_patterns_helpers.h"

std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& geometry, EdgeData<double>& sigma, float period, 
                                                                            std::vector<std::pair<int, int>>& vertexMappingsPairs, std::vector<std::pair<int, int>> edgeMappingsPairs){
  
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
  for (std::pair<int, int> p : edgeMappingsPairs){
    for (Halfedge he : mesh.halfedges()){
      Edge e2 = mesh.edge(p.second);
      sigma_halfedges[e2.halfedge()] = -sigma_halfedges[e2.halfedge()];
      sigma_halfedges[e2.halfedge().twin()] = -sigma_halfedges[e2.halfedge().twin()];
    }
    Edge e1 = mesh.edge(p.first);
    Edge e2 = mesh.edge(p.second);
    sigma_halfedges[e2.halfedge().twin()] = sigma_halfedges[e1.halfedge()];
    sigma_halfedges[e2.halfedge()] = sigma_halfedges[e1.halfedge().twin()];
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

  // this is the standard BFS
  // std::queue<Vertex> Q;
  // Q.push(startVertex);
  // // Floating point thing for knit graph
  // sigma_mod[startVertex] = 0.0 + 1e-16;
  // visited[startVertex] = true;
  // while(!Q.empty()){
  //   Vertex vi = Q.front(); Q.pop();
  //   Halfedge h = vi.halfedge();
  //   do{
  //     Vertex vj = h.twin().vertex();
  //     if (!visited[vj]){
  //       sigma_mod[vj] = fmod(sigma_mod[vi] + sigma_halfedges[h], period);
  //       visited[vj] = true;
  //       Q.push(vj);
  //     }
  //     h = h.twin().next();
  //   }
  //   while( h != vi.halfedge());
  // }

  std::queue<Vertex> Q;
  Q.push(startVertex);
  sigma_mod[startVertex] = 0.0;
  visited[startVertex] = true;
  while(!Q.empty()){
    Vertex vi = Q.front(); Q.pop();
    std::cout << "index of vi " << vi.getIndex() << std::endl;
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
    for (std::pair<int, int> p : vertexMappingsPairs){
      if (vi.getIndex() == p.first){
        Vertex vj = mesh.vertex(p.second);//jump to "stitched" vertex
        Halfedge h = vj.halfedge();
        do{
          Vertex vk = h.twin().vertex();
          if (!visited[vk]){
            sigma_mod[vk] = fmod(sigma_mod[vj] + sigma_halfedges[h], period);
            visited[vk] = true;
            Q.push(vk);
          }
          h = h.twin().next();
        }
        while (h != vj.halfedge());
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

  FaceData<std::vector<double>> faceIsoValues = getFaceIsoValues(geometry, stripeValues, stripesIndices, period);
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


//integrate the 1-form but use an intrinsic geometry representation to handle the stitched edges 
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneFormGluedMesh(IntrinsicGeometryInterface& geometry, EdgeData<double>& sigma, float period, 
                                            std::vector<std::pair<int, int>>& vertexMappingsPairs, std::vector<std::pair<int, int>> edgeMappingsPairs){
        
  SurfaceMesh& mesh = geometry.mesh;
  //store a map from the original index to index in the edge length geometry
  std::map<int, int> originalIndexToELGIndex;
  //number of "unique vertices" i.e., only consider one vertex per stitch
  int numUniqueVertices = 0;
   
  for (Vertex v : mesh.vertices()){
    size_t iV = v.getIndex();
    //iterate over the mappings 
    for (auto p : vertexMappingsPairs){
      if (p.first == iV || p.second == iV){
        if (originalIndexToELGIndex.find(p.first) == originalIndexToELGIndex.end()
        && originalIndexToELGIndex.find(p.second) == originalIndexToELGIndex.end()){
          originalIndexToELGIndex.insert({p.first, numUniqueVertices});
          originalIndexToELGIndex.insert({p.second, numUniqueVertices});
          numUniqueVertices++;
        }
        if (originalIndexToELGIndex.find(p.first) != originalIndexToELGIndex.end()&&
          originalIndexToELGIndex.find(p.second) == originalIndexToELGIndex.end()){
          originalIndexToELGIndex.insert({p.second, originalIndexToELGIndex.at(p.first)});
        }
        if (originalIndexToELGIndex.find(p.second) != originalIndexToELGIndex.end() &&
          originalIndexToELGIndex.find(p.first) == originalIndexToELGIndex.end()){
          originalIndexToELGIndex.insert({p.first, originalIndexToELGIndex.at(p.second)});
        }
      }
    }
    if (originalIndexToELGIndex.find(iV) == originalIndexToELGIndex.end()){
        originalIndexToELGIndex.insert({iV, numUniqueVertices});
        numUniqueVertices++;
    }
  }

  std::vector<std::vector<size_t>> polygons;
  geometry.requireEdgeLengths();
  for (Face f : mesh.faces()){
    if (f.isBoundaryLoop()) continue;
    std::vector<size_t> currPolygon;
    int i = originalIndexToELGIndex[f.halfedge().tailVertex().getIndex()];
    int j = originalIndexToELGIndex[f.halfedge().next().tailVertex().getIndex()];
    int k = originalIndexToELGIndex[f.halfedge().next().next().tailVertex().getIndex()];
    currPolygon.push_back(i);
    currPolygon.push_back(j);
    currPolygon.push_back(k);
    polygons.emplace_back(currPolygon);
  }
  //for some reason this is not manifold and has duplicate edges? 
  SurfaceMesh * gluedMesh = new SurfaceMesh(polygons);
  EdgeData<double> sigma_edges(*gluedMesh);
  for (Edge e1 : mesh.edges()){
    double sigmaFromOriginalMesh = sigma[e1];
    for (Edge e2 : gluedMesh -> edges()){
      //find the corresponding edge in the glued mesh
      if (e2.halfedge().tailVertex().getIndex() == originalIndexToELGIndex[e1.halfedge().tailVertex().getIndex()]
      && e2.halfedge().tipVertex().getIndex() == originalIndexToELGIndex[e1.halfedge().tipVertex().getIndex()]){
        sigma_edges[e2] = sigmaFromOriginalMesh;
      }
    }
  }

  //store the sigma values per halfedge
  HalfedgeData<double> sigma_halfedges(*gluedMesh);
  
  //assign halfedge sigma values taking orientation into account
  for (Halfedge he : gluedMesh -> halfedges()){
  	if (he.orientation()){
  		sigma_halfedges[he] = sigma_edges[he.edge()];
  	}
  	else{
  		sigma_halfedges[he] = -1.0 * sigma_edges[he.edge()];
  	}
  }

  //integrate up the sigma values now onto vertices
  VertexData<double> sigma_mod(*gluedMesh); 
  VertexData<bool> visited(*gluedMesh, false);
  Vertex startVertex;
  for (Vertex v : gluedMesh -> vertices()){
    if (v.isBoundary()){
      startVertex = v;
      break;
    }
  }
  //now integrate the 1-form on the glued mesh
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

  HalfedgeData<double> sigmaHalfedgesOriginalMesh(mesh);
  VertexData<double> sigmaModOriginalMesh(mesh);

  //first copy over the sigma mod values;
  for (Vertex v : mesh.vertices()){
    sigmaModOriginalMesh[v] = sigma_mod[gluedMesh->vertex(originalIndexToELGIndex[v.getIndex()])];
  }
  //then assign values to the halfedges
  for (Halfedge he : mesh.halfedges()){
    if (he.orientation()){
  		sigmaHalfedgesOriginalMesh[he] = sigma[he.edge()];
  	}
  	else{
  		sigmaHalfedgesOriginalMesh[he] = -1.0 * sigma[he.edge()];
  	}
  }

  //update the 1-form for edges that are "stitched together"
  //not convinced I need to do this. In essence just trying to make the values same across edges that are stitched together
  for (std::pair<int, int> p : edgeMappingsPairs){
    Edge e1 = mesh.edge(p.first);
    Edge e2 = mesh.edge(p.second);
    sigmaHalfedgesOriginalMesh[e2.halfedge().twin()] = sigmaHalfedgesOriginalMesh[e1.halfedge()];
    sigmaHalfedgesOriginalMesh[e2.halfedge()] = sigmaHalfedgesOriginalMesh[e1.halfedge().twin()];
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
    double sigma_ij = sigmaHalfedgesOriginalMesh[hij];
    double sigma_jk = sigmaHalfedgesOriginalMesh[hjk];
    double sigma_ki = sigmaHalfedgesOriginalMesh[hki];

    //compute alpha values at triangle corners
    double alphaI = sigmaModOriginalMesh[f.halfedge().vertex()];
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