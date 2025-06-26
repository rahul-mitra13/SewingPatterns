#include "KnitGraph.h"

using namespace std;

KnitGraph::KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, double coursePeriod, double walePeriod,
                      CornerData<double>& courseOneForm, EdgeData<double>& courseSingularEdges, CornerData<double>& waleOneForm, EdgeData<double>& waleSingularEdges,
                      std::map<int, int>& globalToGluedEdgeMap){

    this->globalGeometry = &globalGeometry; 
    this->gluedGeometry = &gluedGeometry;
    this->psMesh = &psMesh; 
    this->coursePeriod = coursePeriod;
    this->walePeriod = walePeriod;
    this->courseOneForm = courseOneForm;
    this->courseSingularEdges = courseSingularEdges;
    this->waleSingularEdges = waleSingularEdges;
    this->waleOneForm = waleOneForm;
    this->globalToGluedEdgeMap = &globalToGluedEdgeMap;
}

void KnitGraph::buildGraph(){

    std::cout << "Building knit graph..." << std::endl;
    SurfaceMesh& gluedMesh = this->gluedGeometry->mesh;
    //convert global edges to glued edge indices 
    EdgeData<double> courseSingularEdgesGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedGeometry, courseSingularEdges, *globalToGluedEdgeMap);
    EdgeData<double> waleSingularEdgesGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedGeometry, waleSingularEdges, *globalToGluedEdgeMap);

    // Flag edges that were glued together
    isGlued = EdgeData<bool>(gluedGeometry->mesh, false);
    for (auto &[globalEdgeID, gluedEdgeID] : *globalToGluedEdgeMap) {
        Edge globalEdge = globalGeometry->mesh.edge(globalEdgeID);
        Edge gluedEdge = gluedGeometry->mesh.edge(gluedEdgeID);
        if (globalEdge.isBoundary() && !gluedEdge.isBoundary())
            isGlued[gluedEdge] = true;
    }

    for (Face f : globalGeometry->mesh.faces()){
        //make the connections on a non-singular face in both directions
        handleCourseNonSingularFaceWaleNonSingularFace(f);
    }

    // Plot vertices
    std::vector<Vector3> vertexPositions;
    for (knitGraphVertex &v : vertices)
        vertexPositions.push_back(v.position);
    polyscope::registerPointCloud("knit graph vertices before merge", vertexPositions)->setPointRadius(0.001)->setEnabled(false);

    //handle the merge in the intrinsic setting 
    intrinsicMerge();

    //reorder vertex indices
    makeRealVertices();

    // Plot real vertices
    std::vector<Vector3> realVertexPositions;
    for (knitGraphVertex &v : realVertices)
        realVertexPositions.push_back(v.position);
    polyscope::registerPointCloud("knit graph vertices (real)", realVertexPositions)->setPointRadius(0.001)->setEnabled(false);

    // std::vector<Vector3> realVertexPositions;
    // for (knitGraphVertex &v : realVertices)
    //     realVertexPositions.push_back(v.position);
    // polyscope::registerPointCloud("real vertices", realVertexPositions)->setPointRadius(0.001);


    // tag the increases and decreases 
    //tagIncreasesDecreases();
    
    //render the knit graph
    renderGraph();

    //trace the short rows 
    traceShortRows();

    //sanity check the graph
    sanityCheck();

    std::cout << "Completed building graph..." << std::endl;


}

void KnitGraph::handleCourseNonSingularFaceWaleNonSingularFace(Face &f){

    double alphaI, alphaJ, alphaK, alpha_min, alpha_max, alpha_start, alpha_end, betaI, betaJ, betaK, beta_min, beta_max, beta_start, beta_end;

    double eps = 1e-12;//for error tolerance

    std::vector<knitGraphVertex> faceVertices;

    //looping stuff and linear system solving stuff
    double RHS_alpha, RHS_beta, detA, a1, b1, c1, a2, b2, c2, j, k, bi, bj, bk;

    //grab the vertices on the face
    Vertex vI = f.halfedge().vertex();
    Vertex vJ = f.halfedge().next().vertex();
    Vertex vK = f.halfedge().next().next().vertex();

    Vector3 pI = globalGeometry->vertexPositions[vI];
    Vector3 pJ = globalGeometry->vertexPositions[vJ];
    Vector3 pK = globalGeometry->vertexPositions[vK];

    Vector3 e1 = pJ - pI; // edge (i,j)
    Vector3 e2 = pK - pI; // edge (i,k)
    Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
    double area = n.norm(); // triangle area

    //grab the alpha values
    alphaI = courseOneForm[f.halfedge().corner()];
    alphaJ = courseOneForm[f.halfedge().next().corner()];
    alphaK = courseOneForm[f.halfedge().next().next().corner()];
    alpha_min = std::min({alphaI, alphaJ, alphaK});
    alpha_max = std::max({alphaI, alphaJ, alphaK});
    a1 = alphaI - alphaK;
    b1 = alphaJ - alphaK;
    c1 = alphaK;
    Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);

    //grab the beta values
    betaI = waleOneForm[f.halfedge().corner()];
    betaJ = waleOneForm[f.halfedge().next().corner()];
    betaK = waleOneForm[f.halfedge().next().next().corner()];
    beta_min = std::min({betaI, betaJ, betaK});
    beta_max = std::max({betaI, betaJ, betaK});
    a2 = betaI - betaK;
    b2 = betaJ - betaK;
    c2 = betaK;
    Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);

    //trace the middle of the stripes
    alpha_start = (std::ceil((alpha_min - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    alpha_end = (std::floor((alpha_max - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    beta_start = (std::ceil((beta_min - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
    beta_end = (std::floor((beta_max - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;

    //trace the stripes
    // alpha_start = (std::ceil((alpha_min - coursePeriod/2.)/coursePeriod/2.) * coursePeriod) + coursePeriod/2.;
    // alpha_end = (std::floor((alpha_max - coursePeriod/2.)/coursePeriod/2.) * coursePeriod) + coursePeriod/2.;
    // beta_start = (std::ceil((beta_min - walePeriod/2.)/walePeriod/2.) * walePeriod) + walePeriod/2.;
    // beta_end = (std::floor((beta_max - walePeriod/2.)/walePeriod/2.) * walePeriod) + walePeriod/2.;

    Eigen::Matrix2f A, A_inv;
    Eigen::Vector2f x, b;

    //basically solving a 2 by 2 linear system over every face
    //create real vertices
   	//shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//step alpha
   	    for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){//step beta
   			//set up the linear system as Ax = b
   			//set RHS to constant
   			RHS_alpha = j - c1;
   			RHS_beta = k - c2;
   			A << a1, b1, a2, b2;
   			//finding the determinant of A
   			detA = (a1 * b2) - (b1 * a2);
   			//find the inverse of a
   			A_inv << (b2 / detA), (-b1/detA), (-a2/detA), (a1/detA);
   			b << RHS_alpha, RHS_beta;
   			x = A_inv * b;
   			bi = x(0);
   			bj = x(1);
   			bk = 1.0 - bi - bj;
   			if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                knitGraphVertex v;						
   			   	v.baryCoords = Vector3{bi, bj, bk};
                v.id = vertexID++;
                v.alpha_tag = j;
                v.beta_tag = k;
                v.face = f;
                v.isVirtual = false;
                v.isAlphaVirtual = false;
                v.isBetaVirtual = false;
                faceVertices.push_back(v);
                //set the position from the bary coords
                v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
                vertices.push_back(v);
                //vertexInfoMap[v.id] = v;
            }
        }
    }

    //interpolate the betas
    //shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//fix alpha
        //ij edge
        bi = (j - alphaJ) / (alphaI - alphaJ);
        bj = 1.0 - bi;
        bk = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            k = bi * betaI + bj * betaJ + bk * betaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }

        //jk edge
        bj = (j - alphaK) / (alphaJ - alphaK);
        bk = 1.0 - bj;
        bi = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            k = bi * betaI + bj * betaJ + bk * betaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }

        //ki edge
        bi = (j - alphaK) / (alphaI - alphaK);
        bk = 1.0 - bi;
        bj = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            k = bi * betaI + bj * betaJ + bk * betaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }
    }

    //interpolate the alphas
    //shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){        
        //ij edge
        bi = (k - betaJ) / (betaI - betaJ);
        bj = 1.0 - bi;
        bk = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            j = bi * alphaI + bj * alphaJ + bk * alphaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }

        //jk edge
        bj = (k - betaK) / (betaJ - betaK);
        bk = 1.0 - bj;
        bi = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            j = bi * alphaI + bj * alphaJ + bk * alphaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }

        //ki edge
        bi = (k - betaK) / (betaI - betaK);
        bk = 1.0 - bi;
        bj = 0.0;
        if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
            knitGraphVertex v;						
   			v.baryCoords = Vector3{bi, bj, bk};
            j = bi * alphaI + bj * alphaJ + bk * alphaK;
            v.id = vertexID++;
            v.alpha_tag = j;
            v.beta_tag = k;
            v.face = f;
            //edges and halfedges are in the glued mesh setting 
            v.edge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next().edge();
            v.halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            v.position = v.baryCoords[0] * pI + v.baryCoords[1] * pJ + v.baryCoords[2] * pK;
            vertices.push_back(v);
            //vertexInfoMap[v.id] = v;
        }
    }

    if (dot(cross(gradAlpha, gradBeta), n) < 0) {
        for (knitGraphVertex &v : faceVertices) {
            v.alpha_tag = -v.alpha_tag;
            v.beta_tag = -v.beta_tag;
        }
    }

    connectOnSmoothFace(faceVertices);
}

//update connections on the vertices of a smooth face 
void KnitGraph::connectOnSmoothFace(std::vector<knitGraphVertex>& faceVertices){

    //3 cases: 
    //a. Connecting real vertices to real vertices
    //b. Connecting real vertices to virtual vertices 
    //c. Connecting virtual vertices to virtual vertices

    std::vector<double> uniqueAlphas;
    std::vector<double> uniqueBetas;
    std::map<double, int> currAlphaRow;
    std::map<double, int> currBetaCol;
    knitGraphVertex* currVertex;
    knitGraphVertex* nextVertex;
    double eps = 1e-8;

    // Collect unique alphas and betas
    for (knitGraphVertex &v : faceVertices) {
        if (!v.isBetaVirtual) {
            bool isUnique = true;
            for (double alpha : uniqueAlphas)
                if (std::abs(alpha - v.alpha_tag) < eps) {
                    isUnique = false;
                    break;
                }
            if (isUnique)
                uniqueAlphas.push_back(v.alpha_tag);
        }

        if (!v.isAlphaVirtual) {
            bool isUnique = true;
            for (double beta : uniqueBetas)
                if (std::abs(beta - v.beta_tag) < eps) {
                    isUnique = false;
                    break;
                }
            if (isUnique)
                uniqueBetas.push_back(v.beta_tag);
        }
    }
   
    
    for (double currAlphaVal : uniqueAlphas){
        currAlphaRow.clear();
        //make an ordered map to store beta values for the current alpha value
        for (knitGraphVertex &v : faceVertices){
            if (std::fabs(v.alpha_tag - currAlphaVal) < eps){
                currAlphaRow[v.beta_tag] = v.id;
            }
        }
        //set the course connections using the map
        for (auto it = currAlphaRow.begin(); it != std::prev(currAlphaRow.end()); it++){
            currVertex = &vertices[it->second];
            nextVertex = &vertices[std::next(it)->second]; 
            currVertex->row_out = nextVertex->id;
            nextVertex->row_in = currVertex->id;
        }
    }
    
    for (double currBetaVal : uniqueBetas){
        //make an ordered map to store the alpha tags for the current beta value
        currBetaCol.clear();
        for (knitGraphVertex &v : faceVertices){
            if (std::fabs(v.beta_tag - currBetaVal) < eps){
                currBetaCol[v.alpha_tag] = v.id;
            }
        }
        //set the wale connections using the map
        for (auto it = currBetaCol.begin(); it != std::prev(currBetaCol.end()); it++){
            currVertex = &vertices[it->second];
            nextVertex = &vertices[std::next(it)->second];
            currVertex->col_out[0] = nextVertex->id;
            nextVertex->col_in[0] = currVertex->id;
        }
    }
}

void KnitGraph::renderGraph(){
    
    //visualize the knit graph vertices with the real connections
    for (auto &v : realVertices){
        // if (v.isVirtual) continue;
        // if (vertices[v.row_out].isVirtual) continue;
        // if (vertices[v.col_out[0]].isVirtual) continue;
        vertexPositions.push_back(v.position);
        if (v.row_out != -1)
            edges.push_back({v.id, v.row_out});
        if (v.col_out[0] != -1)
            edges.push_back({v.id, v.col_out[0]});
        if (v.col_out[1] != -1)
            edges.push_back({v.id, v.col_out[1]});
        if (v.col_in[0] != -1)
            edges.push_back({v.id, v.col_in[0]});
        if (v.col_in[1] != -1)
            edges.push_back({v.id, v.col_in[1]});
    }
    auto graphReal = polyscope::registerCurveNetwork("knit graph with real connections", vertexPositions, edges);
    graphReal -> setRadius(0.001);
    graphReal -> setEnabled(true);

    //visualize the knit graph with virtual connections
    // for (auto &v : vertices){ 
    //     //if (v.isVirtual) continue;
    //     if (v.row_out != -1 ){
    //         vertexPositions.push_back(v.position);
    //         vertexPositions.push_back(vertices[v.row_out].position);
    //     }
    //     if (v.col_out[0] != -1){
    //         vertexPositions.push_back(v.position);
    //         vertexPositions.push_back(vertices[v.col_out[0]].position);
    //     }
    // }

    // for (int i = 0; i < vertexPositions.size(); i+=2){
    //     edges.push_back({i, i + 1});
    // }

    // auto graph = polyscope::registerCurveNetwork("knit graph with virtual connections", vertexPositions, edges);
    // graph -> setRadius(0.001);
    // graph -> setEnabled(false);

}

//find the sign of a value
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


//carry out a merging in the intrinsic setting
void KnitGraph::intrinsicMerge(){

    //convert the singular edges to the glued mesh setting 
    EdgeData<double> courseSingularEdgesGlued(gluedGeometry->mesh, 0.0);
    EdgeData<double> waleSingularEdgesGlued(gluedGeometry->mesh, 0.0);
    courseSingularEdgesGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedGeometry, courseSingularEdges, *globalToGluedEdgeMap);
    waleSingularEdgesGlued = convertGlobalToGluedEdgeFunction(*globalGeometry, *gluedGeometry, waleSingularEdges, *globalToGluedEdgeMap);


    // Store these vertices in order of the "direction of the halfedge"
    std::map<Halfedge, std::vector<knitGraphVertex>> halfedgeCourseVertices;
    std::map<Halfedge, std::vector<knitGraphVertex>> halfedgeWaleVertices;
    for (const knitGraphVertex& v : vertices) if (v.isVirtual) {

        // // Prevent merging on singular edges
        // if (v.isAlphaVirtual && std::fabs(courseSingularEdgesGlued[v.edge.value()]) > 0) continue; // end of course stripe: skip
        // if (v.isBetaVirtual  && std::fabs(waleSingularEdgesGlued[v.edge.value()]) > 0)   continue; // end of wale stripe: skip

        if (v.isAlphaVirtual) halfedgeCourseVertices[v.halfedge.value()].push_back(v);
        if (v.isBetaVirtual) halfedgeWaleVertices[v.halfedge.value()].push_back(v);
    }

    // // Check that each edge has an even number of virtual vertices
    // for (Edge e : (gluedGeometry->mesh).edges()){
    //     if (e.isBoundary()) continue;
    //     if (numVirtualVertices[e] % 2 != 0){
    //         std::cout << "something wrong on edge " << e << std::endl;
    //         std::cout << "non singular edge has an uneven number of virtual vertices " << std::endl;
    //         polyscope::show();
    //         exit(1);    
    //     }
    // }

    // Custom comparator function
    auto compareByJK = [](const knitGraphVertex &a, const knitGraphVertex &b) -> bool {
        //sort in ascending order 
        return a.baryCoords[1] > b.baryCoords[1];
    };

    // Custom comparator function
    auto compareByKI = [](const knitGraphVertex &a, const knitGraphVertex &b) -> bool {
        //sort in ascending order 
        return a.baryCoords[2] > b.baryCoords[2];
    };

    // Custom comparator function
    auto compareByIJ = [](const knitGraphVertex &a, const knitGraphVertex &b) -> bool {
        //sort in ascending order 
        return a.baryCoords[0] > b.baryCoords[0];
    };

    //now sort the halfedge vertices based on their barycoords (course + wale)
    for (auto &[he, verts] : halfedgeCourseVertices) {
        if (std::fabs(verts[0].baryCoords[0]) < 1e-8){//this is a halfedge on the jk edge
            std::sort(verts.begin(), verts.end(), compareByJK);
        }
        else if (std::fabs(verts[0].baryCoords[1]) < 1e-8){//this is a halfedge on the ki edge
            std::sort(verts.begin(), verts.end(), compareByKI);
        }
        else{//this is a halfedge on the ij edge
            std::sort(verts.begin(), verts.end(), compareByIJ);
        }
    }
    for (auto &[he, verts] : halfedgeWaleVertices) {
        if (std::fabs(verts[0].baryCoords[0]) < 1e-8){//this is a halfedge on the jk edge
            std::sort(verts.begin(), verts.end(), compareByJK);
        }
        else if (std::fabs(verts[0].baryCoords[1]) < 1e-8){//this is a halfedge on the ki edge
            std::sort(verts.begin(), verts.end(), compareByKI);
        }
        else{//this is a halfedge on the ij edge
            std::sort(verts.begin(), verts.end(), compareByIJ);
        }
    }

    // Make connections across regular course edges
    for (Edge e : (gluedGeometry->mesh).edges()) if (!e.isBoundary() && courseSingularEdgesGlued[e] == 0) {

        std::vector<knitGraphVertex> he1CourseVertices = halfedgeCourseVertices[e.halfedge()];
        std::vector<knitGraphVertex> he2CourseVertices = halfedgeCourseVertices[e.halfedge().twin()];

        // Trivial matchings
        vector<pair<int, int>> matchings;
        for (int i = 0; i < he1CourseVertices.size(); i++) {
            matchings.push_back({i, (he1CourseVertices.size() - i) - 1});
        }

        // Connect all matchings
        for (auto [i1, i2] : matchings) {
            knitGraphVertex v1 = he1CourseVertices[i1];
            knitGraphVertex v2 = he2CourseVertices[i2];

            // if (v1.col_in[0] == -1 && v2.col_in[0] == -1 && v1.col_out[0] == -1 && v2.col_out[0] == -1){//this is a row merge
            if (v1.row_in == -1 && v2.row_in != -1){
                vertices[v2.id].row_out = v1.id;
                vertices[v1.id].row_in  = v2.id;
            }
            if (v2.row_in == -1 && v1.row_in != -1){
                vertices[v1.id].row_out = v2.id;
                vertices[v2.id].row_in  = v1.id;                    
            }
            // }
        }

    }

    // Now we need to connect across course singular edges.
    // For each positive sing, we choose an arbitrary stripe to start a short row (TODO: make a better choice here).
    // We then propagate it until reaching another sing edge. If that sing edge matches the first one, we're done.
    // If not, we connect it across in a way that respects the ordering.

    map<int, int> matchings; // -1 means short row
    vector<Edge> processedEdges; // for debugging

    // Order positive edges in decreasing order of time
    vector<Edge> orderedPosEdges;
    map<int, Edge, greater<int>> posEdgesByTime;
    for (Edge edge : (gluedGeometry->mesh).edges()) if (!edge.isBoundary() && courseSingularEdgesGlued[edge] > 0)
        posEdgesByTime[courseSingularEdgesGlued[edge]] = edge;
    for (auto &[time,edge] : posEdgesByTime)
        orderedPosEdges.push_back(edge);

    for (Edge startEdge : orderedPosEdges) {
        // Loop on positive course sings

        processedEdges.push_back(startEdge);

        int startEdgeOrder = round(courseSingularEdgesGlued[startEdge]);

        std::vector<knitGraphVertex> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
        std::vector<knitGraphVertex> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];
        if (he1Vertices.size() > he2Vertices.size())
            swap(he1Vertices, he2Vertices);

        ensure(he2Vertices.size() - he1Vertices.size() == 1);

        bool success = false;
        while (!success) {

            // As a starting vertex, pick the one after the last matched vertex on he2Vertices
            knitGraphVertex startVertex = he2Vertices[0];
            for (int i = 0; i < he2Vertices.size()-1; i++) {
                if (matchings.count(he2Vertices[i].id))
                    startVertex = he2Vertices[i+1];
            }
            ensure(!matchings.count(startVertex.id)); // if this fails we're screwed lol

            matchings[startVertex.id] = -1;

            // Trace short row. Now that we go in the direction of row_in (== right)!
            knitGraphVertex walker = startVertex;
            ensure(walker.row_in != -1);
            while (true) {

                if (walker.row_in == -1) {
                    // We hit a singular edge
                    ensure(walker.isVirtual); // walker must be virtual

                    // if (matchings.count(walker.id)) {
                    //     polyscope::registerPointCloud("matched walker", vector<Vector3>{startVertex.position, walker.position});
                    //     showEdges("processedEdges", processedEdges, *globalGeometry);
                    //     H(matchings[walker.id]);
                    //     H(startVertex.id);
                    //     polyscope::show();
                    // }

                    ensure(!matchings.count(walker.id)); // walker must be unmatched
                    ensure(walker.edge.has_value());

                    Edge edge = walker.edge.value();
                    int edgeOrder = round(courseSingularEdgesGlued[edge]);

                    if (edgeOrder == -startEdgeOrder) {
                        // It's a match! We're done
                        matchings[walker.id] = -1;
                        success = true;
                        break;
                    } else {
                        // We need to cross this edge. Above or below?
                        std::vector<knitGraphVertex> leftVertices = halfedgeCourseVertices[edge.halfedge()];
                        std::vector<knitGraphVertex> rightVertices = halfedgeCourseVertices[edge.halfedge().twin()];

                        if (sgn((int)rightVertices.size() - (int)leftVertices.size()) != sgn(edgeOrder))
                            swap(leftVertices, rightVertices);
                        ensure(sgn((int)rightVertices.size() - (int)leftVertices.size()) == sgn(edgeOrder)); // make sure the sides are correct!

                        // Find index along half-edge
                        int indexAlongHalfedge = -1;
                        for (int i = 0; i < leftVertices.size(); i++) {
                            if (leftVertices[i].id == walker.id) {
                                indexAlongHalfedge = i;
                                break;
                            }
                        }

                        knitGraphVertex connectTo;
                        if (abs(edgeOrder) > startEdgeOrder) {
                            // Go below
                            connectTo = rightVertices[rightVertices.size()-1-indexAlongHalfedge];
                        } else {
                            // Go above (also happens if we looped to the starting edge)
                            connectTo = rightVertices[leftVertices.size()-1-indexAlongHalfedge];
                        }
                        matchings[walker.id] = connectTo.id;
                        matchings[connectTo.id] = walker.id;
                        walker = connectTo;
                        if (edgeOrder == startEdgeOrder) {
                            // We looped around! Break and continue onto the next short row candidate
                            break;
                        }
                    }
                } else {
                    walker = vertices[walker.row_in];
                }

            }

        }
    }



    // Connect the remaining unmatched virtual vertices across singular edges
    for (Edge startEdge : (gluedGeometry->mesh).edges()) if (!startEdge.isBoundary() && courseSingularEdgesGlued[startEdge] != 0) {
        std::vector<knitGraphVertex> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
        std::vector<knitGraphVertex> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];

        int i2 = he2Vertices.size()-1;
        for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
            if (!matchings.count(he1Vertices[i1].id)) {
                while (matchings.count(he2Vertices[i2].id)) {
                    i2--;
                    if (i2 < 0) {
                        showEdges("i2 < 0", {startEdge}, *globalGeometry);
                        // polyscope::registerPointCloud("i2 < 0", vector<Vector3> {he2Vertices[i2].position});
                        polyscope::show();
                    }
                    ensure (i2 >= 0);
                }
                matchings[he1Vertices[i1].id] = he2Vertices[i2].id;
                matchings[he2Vertices[i2].id] = he1Vertices[i1].id;
                i2--;
            }
        }
    }

    // Now actually connect all course matchings
    for (auto [i1, i2] : matchings) {
        if (i1 == -1 || i2 == -1) continue;
        knitGraphVertex v1 = vertices[i1];
        knitGraphVertex v2 = vertices[i2];
        if (v1.row_in == -1 && v2.row_in != -1){
            vertices[v2.id].row_out = v1.id;
            vertices[v1.id].row_in  = v2.id;
        }
        if (v2.row_in == -1 && v1.row_in != -1){
            vertices[v1.id].row_out = v2.id;
            vertices[v2.id].row_in  = v1.id;                    
        }
    }


    // Connect wale vertices 
    for (Edge e : (gluedGeometry->mesh).edges()) {
        if (e.isBoundary()) continue;//don't need to handle boundary vertices

        std::vector<knitGraphVertex> he1Vertices = halfedgeWaleVertices[e.halfedge()];
        std::vector<knitGraphVertex> he2Vertices = halfedgeWaleVertices[e.halfedge().twin()];

        // classical singularity with 1 stripe being born/dying: nothing to do
        if (he1Vertices.size() + he2Vertices.size() == 1)
            continue;

        vector<pair<int, int>> matchings;

        // more funky singularity: we need to match
        if (he1Vertices.size() != he2Vertices.size()) { 

            // Check which side of the triangles he1 and he2 are located
            // TODO: find a way to do this without epsilons
            int side1, side2;
            for (int i = 0; i < 3; i++) {
                if (abs(he1Vertices[0].baryCoords[i]) < 1e-8)
                    side1 = (i+1)%3;
                if (abs(he2Vertices[0].baryCoords[i]) < 1e-8)
                    side2 = (i+1)%3;
            }

            // Fetch coordinate along edge
            vector<double> coordsAlongHe1, coordsAlongHe2;
            for (knitGraphVertex &v : he1Vertices)
                coordsAlongHe1.push_back(v.baryCoords[side1]);
            for (knitGraphVertex &v : he2Vertices)
                coordsAlongHe2.push_back(1 - v.baryCoords[side2]); // need to invert to be in the same basis

            // Find out best matchings (greedy approach)
            if (he1Vertices.size() < he2Vertices.size()) { // match every vertex of he1 to closest vertex of he2
                for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
                    int i2closest = -1;
                    for (int i2 = 0; i2 < he2Vertices.size(); i2++)
                        if (i2closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1] - coordsAlongHe2[i2closest]))
                            i2closest = i2;
                    if (i2closest == -1) {
                        cout << "Error: No match found for vertex " << he1Vertices[i1].id << endl;
                        polyscope::show();
                    } else {
                        matchings.push_back({i1, i2closest});
                    }
                }
            } else { // match every vertex of he2 to closest vertex of he1
                for (int i2 = 0; i2 < he2Vertices.size(); i2++) {
                    int i1closest = -1;
                    for (int i1 = 0; i1 < he1Vertices.size(); i1++)
                        if (i1closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1closest] - coordsAlongHe2[i2]))
                            i1closest = i1;
                    if (i1closest == -1) {
                        cout << "Error: No match found for vertex " << he2Vertices[i2].id << endl;
                        polyscope::show();
                    } else {
                        matchings.push_back({i1closest, i2});
                    }
                }
            }
        } else {
            // Edge is regular: matching is trivial
            for (int i = 0; i < he1Vertices.size(); i++) {
                matchings.push_back({i, (he1Vertices.size() - i) - 1});
            }
        }

        // Connect the wale matchings we found
        for (auto [i1, i2] : matchings) {
            knitGraphVertex v1 = vertices[he1Vertices[i1].id];
            knitGraphVertex v2 = vertices[he2Vertices[i2].id];
            if (v1.col_in[0] == -1 && v2.col_in[0] != -1){
                vertices[v1.id].col_in[0]  = v2.id;
                vertices[v2.id].col_out[0] = v1.id;
            }
            if (v2.col_in[0] == -1 && v1.col_in[0] != -1){
                vertices[v1.id].col_out[0] = v2.id;
                vertices[v2.id].col_in[0]  = v1.id;
            }
        }

    }

    // Now connect real vertices to one another
    for (knitGraphVertex v0 : vertices) if (!v0.isVirtual) {
        knitGraphVertex v, vprev;
        bool isGluedPath;
        
        // Find row_in and row_out
        ensure(v0.row_out != -1);
        v = vertices[v0.row_out];
        isGluedPath = false;
        while(v.isVirtual) {
            if (isGlued[v.edge.value()]) isGluedPath = true;
            if (v.row_out == -1) break; // sing or bdy edge
            v = vertices[v.row_out]; 
        }
        if (v.isVirtual) // ended up on a sing or bdy edge
            vertices[v0.id].row_out = -1;
        else { // v is real: new connection
            vertices[v0.id].row_out = v.id;
            vertices[v.id].row_in = v0.id;
            if (isGluedPath) stitchedVertices.push_back({v0.id, v.id});
        }

        // Find col_in[0] and col_out[0]
        ensure(v0.col_out[0] != -1);
        v = vertices[v0.col_out[0]];
        isGluedPath = false;
        while (v.isVirtual) {
            if (isGlued[v.edge.value()]) isGluedPath = true;
            if (v.col_out[0] == -1) break; // sing or bdy edge
            vprev = v; // save previous
            v = vertices[v.col_out[0]];
        }
        if (v.isVirtual)
            vertices[v0.id].col_out[0] = -1;
        else { // v is real: new connection
            vertices[v0.id].col_out[0] = v.id;
            vertices[v.id].col_in[0] = v0.id;
            if (isGluedPath) stitchedVertices.push_back({v0.id, v.id});
        }        
    }

    //remove any lingering connections from real vertices to virtual vertices 
    for (knitGraphVertex &v : vertices){
        if (vertices[v.row_out].isVirtual) v.row_out = -1;
        if (vertices[v.row_in].isVirtual) v.row_in = -1;
        if (vertices[v.col_in[0]].isVirtual) v.col_in[0] = -1;
        if (vertices[v.col_out[0]].isVirtual) v.col_out[0] = -1;
    }

    // // Print out real vertices
    // for (knitGraphVertex& v: vertices) if (!v.isVirtual) {
    //     std::cout << "id = " << v.id << std::endl;
    //     std::cout << "row in = " << v.row_in << std::endl;
    //     std::cout << "row out = " << v.row_out << std::endl;
    //     std::cout << "col_in[0] = " << v.col_in[0] << std::endl;
    //     std::cout << "col_out[0] = " << v.col_out[0] << std::endl;
    //     std::cout << "---------------" << std::endl;
    // }
}

//perform epsilon merging to 
//merge connections across faces
void KnitGraph::epsilonMerging(){

    std::vector<Vector3> badVirtuals;
    for (knitGraphVertex &v : vertices){
        if (!v.isVirtual) continue;
        if ((v.row_in == -1 && v.row_out == -1) 
            && (v.col_in[0] == -1 && v.col_out[0] == -1)){
                badVirtuals.push_back(v.position);
        }
    }
    auto PC = polyscope::registerPointCloud("bad virtuals", badVirtuals);
    //PC->setRadius(0.001);

    //epsilon ball around vertex we hope to merge
    const double eps = 1e-8;
    for (knitGraphVertex &v : vertices) {
        //skip singular edges in the intial merging
        // if (v.edge.has_value() && (std::fabs(courseSingularEdges[v.edge.value()]) > 0 || std::fabs(waleSingularEdges[v.edge.value()]) > 0)){ 
        //     //std::cout << "skipping singular edge " << v.edge.value() << std::endl;
        //     continue;
        // }

        // find clusters for all vertices that have not been handled
        if (!v.hasBeenHandled) {
            auto vCluster = findCluster(v, eps);
            if (vCluster.size() == 1) {
                v.hasBeenHandled = true;
                continue;
            }
            if (vCluster.size() != 2){
                std::cout << "size of cluster is not 2 " << std::endl;
                //exit(1);
            }
            mergeCluster(vCluster, eps);
        }
    }

    //remove any lingering connections from real vertices to virtual vertices 
    for (knitGraphVertex &v : vertices){
        if (vertices[v.row_out].isVirtual) v.row_out = -1;
        if (vertices[v.row_in].isVirtual) v.row_in = -1;
        if (vertices[v.col_in[0]].isVirtual) v.col_in[0] = -1;
        if (vertices[v.col_out[0]].isVirtual) v.col_out[0] = -1;
    }
}

//find a cluster of vertices to merge
std::vector<knitGraphVertex> KnitGraph::findCluster(const knitGraphVertex &v, double eps){
    std::vector<knitGraphVertex> vCluster;
    const auto vPos = v.position;
    for (knitGraphVertex &w : vertices) {
        const auto wPos = w.position;
        // ensure 'w' is within epsilon of 'v'
        if (norm(vPos - wPos) < eps && !w.hasBeenHandled){
            w.hasBeenHandled = true;
            vCluster.push_back(w);
        }
    }
    return vCluster;
}

//merge a cluster
void KnitGraph::mergeCluster(std::vector<knitGraphVertex>& vCluster, double eps){

    //std::cout << "size of cluster = " << vCluster.size() << std::endl;
    bool hasRealVertex = false; 
    int id_to_keep = -1;
    //check if this cluster has a real vertex 
    //in that case we need to preserve the real vertex in this cluster.
    for (auto &v : vCluster) {
        if (!v.isVirtual) {
            hasRealVertex = true;
            id_to_keep = v.id;
            break;
        }
    }

    if (hasRealVertex){//if the cluster has a real vertex

        int global_row_in = -1;
        int global_row_out = -1;
        int global_col_in = -1;
        int global_col_out = -1;

        // find the row in/out and col in/outs
        for (auto &vi : vCluster){
            for (auto &vj : vertices){
                if (norm(vi.position - vj.position) > eps && (vi.row_in == vj.id)){
                    global_row_in = vj.id;
                }
                if (norm(vi.position - vj.position) > eps && (vi.row_out == vj.id)){
                    global_row_out = vj.id;
                }
                if (norm(vi.position - vj.position) > eps && (vi.col_in[0] == vj.id)){
                    global_col_in = vj.id;
                }
                if (norm(vi.position - vj.position) > eps && (vi.col_out[0] == vj.id)){
                    global_col_out = vj.id;
                }
            }
        }

        //set all the vertices in this cluster to virtual 
        //except for the vertex we want to keep
        // for (auto &vi : vCluster){
        //     if (vi.id != id_to_keep){
        //         vi.isVirtual = true;
        //     }
        // }

        // auto &realVertex = vertices[id_to_keep];

        // std::cout << "global row in = " << global_row_in << std::endl;
        // std::cout << "global row out = " << global_row_out << std::endl;
        // std::cout << "global col in = " << global_col_in << std::endl;
        // std::cout << "global col out = " << global_col_out << std::endl;
        // std::cout << "----------------------------" << std::endl;
 
        // realVertex.row_in = global_row_in;
        // vertices[global_row_in].row_out = realVertex.id;

        // realVertex.row_out = global_row_out;
        // vertices[global_row_out].row_in = realVertex.id;

        // realVertex.col_in[0] = global_col_in;
        // vertices[global_col_in].col_out[0] = realVertex.id;

        // realVertex.col_out[0] = global_col_out;
        // vertices[global_col_out].col_in[0] = realVertex.id;

        if (global_row_in >=0 && global_row_in < vertices.size()){
            if (global_row_out >=0 && global_row_out < vertices.size()){//bounds checking
                vertices[id_to_keep].row_in = global_row_in;
                vertices[global_row_in].row_out = id_to_keep;

                vertices[id_to_keep].row_out = global_row_out;
                vertices[global_row_out].row_in = id_to_keep;
            }
        }

        if (global_col_in >= 0 && global_col_in < vertices.size()){
            if (global_col_out >= 0 && global_col_out < vertices.size()){//bounds checking
                vertices[id_to_keep].col_in[0] = global_col_in;
                vertices[global_col_in].col_out[0] = id_to_keep;

                vertices[id_to_keep].col_out[0] = global_col_out;
                vertices[global_col_out].col_in[0] = id_to_keep;
            }
        }


    }

    else{
        
        int global_row_in = -1;
        int global_row_out = -1;
        int global_col_in = -1;
        int global_col_out = -1;

        for (auto &vi : vCluster){ // for every vertex in the current cluster
            // for every vertex outside of the current cluster
            for (auto &vj : vertices) if (norm(vi.position - vj.position) > 2*eps) {
                if (vi.row_in == vj.id){
                    //std::cout << "here for global row in!! " << std::endl;
                    global_row_in = vj.id;
                }
                if (vi.row_out == vj.id){
                    //std::cout << "here for global row out!! " << std::endl;
                    global_row_out = vj.id;
                }
                if (vi.col_in[0] == vj.id){
                    //std::cout << "here for global col in!! " << std::endl;
                    global_col_in = vj.id;
                }
                if (vi.col_out[0] == vj.id){
                    //std::cout << "here for global col out!! " << std::endl;
                    global_col_out = vj.id;
                }
            }
        }

        //unset the connections for the virtual vertices in the cluster
        // for (auto &vi : vCluster){
        //     vi.row_in = -1;
        //     vi.row_out = -1;
        //     vi.col_in[0] = -1;
        //     vi.col_in[1] = -1;
        //     vi.col_out[0] =  -1;
        //     vi.col_out[1] = -1;
        // }

        // std::cout << "global row in = " << global_row_in << std::endl;
        // std::cout << "global row out = " << global_row_out << std::endl;
        // std::cout << "global col in = " << global_col_in << std::endl;
        // std::cout << "global col out = " << global_col_out << std::endl;
        // std::cout << "----------------------------" << std::endl;

        //shouldn't really matter the order in which merging is happening
        if (global_row_in >= 0 && global_row_in < vertices.size()){
            if (global_row_out >= 0 && global_row_out < vertices.size()){//bounds checking
                vertices[global_row_in].row_out = global_row_out;
                vertices[global_row_out].row_in = global_row_in;
            }
        }
        if (global_col_in >= 0 && global_col_in < vertices.size()){
            if (global_col_out >= 0 && global_col_out < vertices.size()){//bounds checking
                vertices[global_col_in].col_out[0] = global_col_out;
                vertices[global_col_out].col_in[0] = global_col_in;
            }
        }
        
    }

}

void KnitGraph::makeRealVertices(){

    //first re-order indices in the graph
    //convert graph indices to autoknit txt format 
    std::map<int, int> mp;
    int i = 0;
    //std::cout << "pre merge info: " << std::endl;
    for (knitGraphVertex& v: vertices){
        if (v.isVirtual) continue;//only handle real vertices
        mp.insert({v.id, i++});
        // std::cout << "id = " << v.id << std::endl;
        // std::cout << "row in = " << v.row_in << std::endl;
        // std::cout << "row out = " << v.row_out << std::endl;
        // std::cout << "col_in[0] = " << v.col_in[0] << std::endl;
        // std::cout << "col_out[0] = " << v.col_out[0] << std::endl;
        // std::cout << "---------------" << std::endl;
    }

    int id = -1;
    int row_in = -1;
    int row_out = -1;
    int col_in_1 = -1;
    int col_in_2 = -1;
    int col_out_1 = -1;
    int col_out_2 = -1;

    //resize the real vertices
    realVertices.resize(mp.size());

    std::cout << "size of real vertices = " << realVertices.size() << std::endl;

    //update the ids in the knit graph
    for (knitGraphVertex& v : vertices){

        if (v.isVirtual) continue;//only handle real vertices

        int id = mp[v.id];

        if (v.row_in == -1){//row_in is unset
            row_in = -1;
        }
        else{
            // it = mp.find(v.row_in);//find the row_in from the original graph
            // row_in = it -> second;
            row_in = mp[v.row_in];
        }

        if (v.row_out == -1){//row_out is unset
            row_out = -1;
        }
        else{
            // it = mp.find(v.row_out);//find the row out from the original graph
            // row_out = it -> second;
            row_out = mp[v.row_out];
        }

        if (v.col_in[0] == -1){//col_in_1 is unset
            col_in_1 = -1;
        }
        else{
            // it = mp.find(v.col_in[0]);//find the col_in_1 from the orginal graph
            // col_in_1 = it -> second;
            col_in_1 = mp[v.col_in[0]];
        }

        if (v.col_in[1] == -1){//col_in_2 is unset
            col_in_2 = -1;
        }
        else{
            // it = mp.find(v.col_in[1]);//find the col_in_1 from the orginal graph
            // col_in_2 = it -> second;
            col_in_2 = mp[v.col_in[1]];
        }

        if (v.col_out[0] == -1){//col_out_1 is unset
            col_out_1 = -1;
        }
        else{
            // it = mp.find(v.col_out[0]);//find the col_out_1 from the orginal graph
            // col_out_1 = it -> second;
            col_out_1 = mp[v.col_out[0]]; 
        }

        if (v.col_out[1] == -1){//col_out_2 is unset
            col_out_2 = -1;
        }
        else{
            // it = mp.find(v.col_out[1]);//find the col_out_1 from the orginal graph
            // col_out_2 = it -> second;
            col_out_2 = mp[v.col_out[1]];
        }

        //update the info
        v.id = id;
        v.row_in = row_in;
        v.row_out = row_out;
        v.col_in[0] = col_in_1;
        v.col_in[1] = col_in_2;
        v.col_out[0] = col_out_1;
        v.col_out[1] = col_out_2;
        //add the updated id vertex to the vector of real vertices
        realVertices[v.id] = v;
    }

    // Map ID's of stitched vertices pairs
    for (std::pair<int,int> &p : stitchedVertices) {
        p.first  = mp[p.first];
        p.second = mp[p.second];
    }

    // for (auto &v : realVertices){
    //     std::cout << "id = " << v.id << std::endl;
    //     std::cout << "row in = " << v.row_in << std::endl;
    //     std::cout << "row out = " << v.row_out << std::endl;
    //     std::cout << "col_in[0] = " << v.col_in[0] << std::endl;
    //     std::cout << "col_out[0] = " << v.col_out[0] << std::endl;
    //     std::cout << "---------------" << std::endl;
    // }
}


void KnitGraph::tagIncreasesDecreases(){
               
    //handle increases
    for (auto &v : realVertices){
        if (v.col_in[0] == -1  && v.col_in[1] == -1){
            if ((realVertices[v.row_out].col_in[0]) != -1 && (realVertices[v.row_in].col_in[0]) != -1){//we're not at the bottom-most course row (assuming can't have short-rows on the bottom row)
                std::cout << "In increase case for vertex " << v.id << std::endl;
                auto &candidate1 = realVertices[realVertices[v.row_out].col_in[0]];
                auto &candidate2 = realVertices[realVertices[v.row_in].col_in[0]];
                if (norm(v.position - candidate1.position) < norm(v.position - candidate2.position)){
                    candidate1.col_out[1] = v.id;
                    v.col_in[0] = candidate1.id;
                }
                else{
                    candidate2.col_out[1] = v.id;
                    v.col_in[0] = candidate2.id;
                }
            }
        }
    }

    //handle decrease
    for (auto &v : realVertices){
        if (v.col_out[0] == -1 && v.col_out[1] == -1){
            if ((realVertices[v.row_out].col_out[0]) != -1 && (realVertices[v.row_in].col_out[0]) != -1){//we're not at the top-most course row (assuming can't have short-rows on the top row)
            std::cout << "In decrease case for vertex " << v.id << std::endl;
            auto &candidate1 = realVertices[realVertices[v.row_out].col_out[0]];
            auto &candidate2 = realVertices[realVertices[v.row_in].col_out[0]];
            if (norm(v.position - candidate1.position) < norm(v.position - candidate2.position)){
                candidate1.col_in[1] = v.id;
                v.col_out[0] = candidate1.id;
            }
            else{
                candidate2.col_in[1] = v.id;
                v.col_out[0] = candidate2.id;
                }
            }
        }
    }

    //check if col_out[0] and col_out[1] need flipping for increases i.e., col_out[0] should always have a row_out that points to col_out[1]
    for (auto &v : realVertices){
        if (v.col_out[0] != -1 && v.col_out[1] != -1){//found a vertex with an increase
            if (realVertices[v.col_out[1]].row_out == realVertices[v.col_out[0]].id){//col_out[1] has a row_out to col_out[0] - need to flip
                std::cout << "flipping col outs for increase vertex " << v.id << std::endl;
                int x = v.col_out[0];
                int y = v.col_out[1];
                v.col_out[0] = y;
                v.col_out[1] = x;
            }
        }
    }

    //check if col_in[0] and col_in[1] need flipping for decreases i.e., col_in[0] should always have a row_out that points to col_in[1]
    for (auto &v : realVertices){
        if (v.col_in[0] != -1 && v.col_in[1] != -1){//found a vertex with a decrease
            if (realVertices[v.col_in[1]].row_out == realVertices[v.col_in[0]].id){//col_in[1] has a row_out to col_in[0] - need to flip
                std::cout << "flipping col ins for decrease vertex " << v.id << std::endl;
                int x = v.col_in[0];
                int y = v.col_in[1];
                v.col_in[0] = y;
                v.col_in[1] = x;
            }
        }
    }

}


//----------------------helper functions----------------------//
//hash function for a floating point number
int KnitGraph::hashFloat(double value) {
    // Reinterpret the double as an integer and hash it
    std::uint64_t intValue = *reinterpret_cast<std::uint64_t*>(&value);
    return std::hash<std::uint64_t>()(intValue);
}

//get a knitgraph vertex by id
knitGraphVertex& KnitGraph::get(int id){
        
   return vertices[id];

}

//sanity check the graph
void KnitGraph::sanityCheck(){

    double eps = 1e-8;

    std::vector<Vector3> buggy_vertices;

    //make sanity check similar to last project 
    int num_errors = 0;
    for (auto &v : realVertices){
        if (v.row_in != -1 && realVertices[v.row_in].row_out != v.id){
            buggy_vertices.push_back(v.position);
            num_errors++;
            std::cout << "Row mismatch at vertex " << v.id << std::endl;
        }
        if (v.row_out != -1 && realVertices[v.row_out].row_in != v.id){
            buggy_vertices.push_back(v.position);
            num_errors++;
            std::cout << "Row mismatch at vertex " << v.id << std::endl;
        }
        if (v.col_in[0] != -1){
            if (realVertices[v.col_in[0]].col_out[0] != v.id && realVertices[v.col_in[0]].col_out[1] != v.id){
                buggy_vertices.push_back(v.position);
                num_errors++;
                std::cout << "Column 0 mismatch at vertex " << v.id << std::endl;
                H(v.id);
                H(v.col_in[0]);
                H(v.col_in[1]);
                H(v.col_out[0]);
                H(v.col_out[1]);
                H(realVertices[v.col_in[0]].col_out[0]);
                H(realVertices[v.col_in[0]].col_out[1]);
            }
        }
        if (v.col_out[0] != -1){
            if (realVertices[v.col_out[0]].col_in[0] != v.id && realVertices[v.col_out[0]].col_in[1] != v.id){
                buggy_vertices.push_back(v.position);
                num_errors++;
                std::cout << "Column 0 mismatch at vertex " << v.id << std::endl;
                H(v.id);
                H(v.col_in[0]);
                H(v.col_in[1]);
                H(v.col_out[0]);
                H(v.col_out[1]);
                H(realVertices[v.col_out[0]].col_in[0]);
                H(realVertices[v.col_out[0]].col_in[1]);
            }
        }

        if (v.col_in[1] != -1){
            if (realVertices[v.col_in[1]].col_out[0] != v.id && realVertices[v.col_in[1]].col_out[1] != v.id){
                buggy_vertices.push_back(v.position);
                num_errors++;
                std::cout << "Column 1 mismatch at vertex " << v.id << std::endl;
            }
        }
        if (v.col_out[1] != -1){
            if (realVertices[v.col_out[1]].col_in[0] != v.id && realVertices[v.col_out[1]].col_in[1] != v.id){
                buggy_vertices.push_back(v.position);
                num_errors++;
                std::cout << "Column 1 mismatch at vertex " << v.id << std::endl;
            }
        }
    }

    // // This is super expensive, disable for now
    // for (auto const &vi : realVertices){
    //     for (auto &vj : realVertices){
    //         if (norm(vi.position - vj.position) < eps && vi.id != vj.id){
    //             std::cout << "2 real vertices are in the same place " << std::endl;
    //         } 
    //     }
    // }

    if (num_errors == 0){
        std::cout << "Sanity check passed " << std::endl << std::endl;
    }
    else{
        std::cout << "Sanity check failed " << std::endl << std::endl;
    }

    polyscope::registerPointCloud("buggy vertices ", buggy_vertices);

}



//write knit graph to txt file 
void KnitGraph::writeKnitGraphToTxtFile(const std::string& file_name){
    
    std::string obj_name = file_name.c_str();

    //write knit graph to txt file 
    std::ofstream file(obj_name.erase(obj_name.length() - 4) + "_knitgraph.txt");
    

    for (auto &v : realVertices){
        if (v.isVirtual) continue;
        file << v.id << " " << v.position[0] << " " << v.position[1] << " " << v.position[2] << " " << v.row_in << " " << v.row_out << 
        " " << v.col_in[0] << " " << v.col_in[1] << " " << v.col_out[0] << " " << v.col_out[1] << "\n";
    }

    // Write pairs of stitched vertices
    // for (const auto &[v1,v2] : stitchedVertices)
    //     file << "s " << v1 << " " << v2 << "\n";

    // file.close();
    std::cout << "wrote knit graph to txt file " << std::endl;
}

//write knit graph as line element 
void KnitGraph::writeKnitGraphLineElement(){

   std::ofstream outfile("lineElement.obj");
    for (auto &v : realVertices){
        outfile << " v " << v.position.x << " " << v.position.y << " " << v.position.z << std::endl;
    }

    for (auto &v : realVertices){
        if (v.row_out != -1) outfile << "l " << v.id + 1 << " " << v.row_out + 1 << std::endl;
        if (v.col_out[0] != -1) outfile << "l " << v.id + 1  << " " << v.col_out[0] + 1 << std::endl;
        if (v.col_out[1] != -1) outfile << "l " << v.id + 1 << " " << v.col_out[1] + 1 << std::endl;
    }

    outfile.close();
    std::cout << "wrote knit graph text file " << std::endl;
}

//trace the short rows in the graph 
void KnitGraph::traceShortRows(){   

    int ctr = 0;
    for (auto &v : realVertices){
        if (v.row_in == -1){
            std::vector<Vector3> pos;
            std::vector<std::array<int, 2>> edges;
            auto walker = v;
            while(walker.row_out != -1){
                pos.push_back(walker.position);
                walker = realVertices[walker.row_out];
            }
            for (int i = 0; i < (int)pos.size() - 1; i++){
                edges.push_back(std::array{i, i + 1});
            }
            if (pos.size() > 1000) return;
            polyscope::registerCurveNetwork("traced short row " + std::to_string(ctr), pos, edges)->setRadius(0.00125);
            ctr++;
        }
    }
}


