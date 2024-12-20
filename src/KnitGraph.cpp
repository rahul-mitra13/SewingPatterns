#include "KnitGraph.h"

KnitGraph::KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, double period, 
                      CornerData<double>& courseOneForm, EdgeData<double>& courseSingularEdges, CornerData<double>& waleOneForm, EdgeData<double>& waleSingularEdges,
                      std::map<int, int>& globalToGluedEdgeMap){

    this->globalGeometry = &globalGeometry; 
    this->gluedGeometry = &gluedGeometry;
    this->psMesh = &psMesh; 
    this->period = period; 
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

    for (Face f : globalGeometry->mesh.faces()){
        //make the connections on a non-singular face in both directions
        handleCourseNonSingularFaceWaleNonSingularFace(f);
    }

    //merge vertices together 
    epsilonMerging();

    //edgeMerging();

    //merge the virtual vertices 
    //mergeVirtual();

    //reorder vertex indices
    makeRealVertices();
    
    //render the knit graph
    renderGraph();

    //sanity check the graph
    //sanityCheck();


    //make the obj
    //makeObj();

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
    alpha_start = (std::ceil((alpha_min - period/4.)/period) * period) + period/4.;
    alpha_end = (std::floor((alpha_max - period/4.)/period) * period) + period/4.;
    beta_start = (std::ceil((beta_min - period/4.)/period) * period) + period/4.;
    beta_end = (std::floor((beta_max - period/4.)/period) * period) + period/4.;

    // alpha_start = std::ceil(alpha_min/period) * period;
    // alpha_end = std::floor(alpha_max/period) * period;
    // beta_start = std::ceil(beta_min/period) * period;
    // beta_end = std::floor(beta_max/period) * period;


    Eigen::Matrix2f A, A_inv;
    Eigen::Vector2f x, b;

    //basically solving a 2 by 2 linear system over every face
    //create real vertices
   	//shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (j = alpha_start - eps; j < alpha_end + eps; j += period){//step alpha
   	    for (k = beta_start - eps; k < beta_end + eps; k += period){//step beta
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
   	for (j = alpha_start - eps; j < alpha_end + eps; j += period){//fix alpha
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
            v.edge = f.halfedge().edge();
            v.halfedge = f.halfedge();
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
            v.edge = f.halfedge().next().edge();
            v.halfedge = f.halfedge().next();
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
            v.edge = f.halfedge().next().next().edge();
            v.halfedge = f.halfedge().next().next();
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
   	for (k = beta_start - eps; k < beta_end + eps; k += period){        
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
            v.edge = f.halfedge().edge();
            v.halfedge = f.halfedge();
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
            v.edge = f.halfedge().next().edge();
            v.halfedge = f.halfedge().next();
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
            v.edge = f.halfedge().next().next().edge();
            v.halfedge = f.halfedge().next().next();
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

void KnitGraph::mergeVirtual(){

    double eps = 1e-8;
    int numCourseErrors = 0;
    int numWaleErrors = 0;

    //merged vertices
    std::map<int, knitGraphVertex> mergedVertices;
    //final map 
    std::map<int, knitGraphVertex> finalMap;


    for (auto vertexI : vertexInfoMap){
        knitGraphVertex v = vertexI.second;
        if (!v.isVirtual) continue;//skip real vertices for now 
        if (v.edge.has_value() && (std::fabs(courseSingularEdges[v.edge.value()]) > eps || std::fabs(waleSingularEdges[v.edge.value()]) > eps
        || v.edge->isBoundary()))  continue;//skip singular and boundary edges for now
        int courseMerge = 0;
        int waleMerge = 0;
        for (auto vertexJ : vertexInfoMap){
            knitGraphVertex k = vertexJ.second;
            if (!k.isVirtual) continue;//skip real vertices for now
            if (v.isAlphaVirtual){
                if (norm(v.position - k.position) < eps && k.isAlphaVirtual && v.id != k.id){
                    courseMerge++;
                    //make a new vertex that represents the merged vertex
                    knitGraphVertex vMerged;
                    vMerged.id = v.id;
                    vMerged.position = v.position;
                    //handle the row_outs first
                    if (v.row_out == -1 && k.row_out == -1)
                        vMerged.row_out = -1;
                    else if (v.row_out != -1 && k.row_out == -1)
                        vMerged.row_out = v.row_out;
                    else if (v.row_out == -1 && k.row_out != -1)
                        vMerged.row_out = k.row_out;
                
                    //handle the row ins 
                    if (v.row_in == -1 && k.row_in == -1)
                        vMerged.row_in = -1;
                    else if (v.row_in != -1 && k.row_in == -1)
                        vMerged.row_in = v.row_in;
                    else if (v.row_in == -1 && k.row_in != -1)
                        vMerged.row_in = k.row_in;
                    
                    mergedVertices[v.id] = vMerged;
                    
                }
            }
            else if (v.isBetaVirtual){
                if (norm(v.position - k.position) < eps && k.isBetaVirtual && v.id != k.id){
                    waleMerge++;
                    //make a new vertex that represents the merged vertex
                    knitGraphVertex vMerged;
                    vMerged.id = v.id;
                    vMerged.position = v.position;
                    //handle the col_outs first
                    if (v.col_out[0] == -1 && k.col_out[0] == -1)
                        vMerged.col_out[0] = -1;
                    else if (v.col_out[0] != -1 && k.col_out[0] == -1)
                        vMerged.col_out[0] = v.col_out[0];
                    else if (v.col_out[0] == -1 && k.col_out[0] != -1)
                        vMerged.col_out[0] = k.col_out[0];
                    //handle the col_ins 
                    if (v.col_in[0] == -1 && k.col_in[0] == -1)
                        vMerged.col_in[0] = -1;
                    else if (v.col_in[0] != -1 && k.col_in[0] == -1)
                        vMerged.col_in[0] = v.col_in[0];
                    else if (v.col_in[0] == -1 && k.col_in[0] != -1)
                        vMerged.col_in[0] = k.col_in[0];
                    
                    mergedVertices[v.id] = vMerged;
                }
            }
        }
        if (courseMerge != 1 && waleMerge != 1){
            std::cout << "error in merge " << std::endl;
            std::cout << "course merge ctr = " << courseMerge << std::endl;
            std::cout << "wale merge ctr = " << waleMerge << std::endl;
            std::cout << "---------------------" << std::endl;
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

//merge virtual vertices across an edge
void KnitGraph::edgeMerging(){

    double eps = 1e-8;
    //map virtual vertices to edges first 
    HalfedgeData<std::vector<knitGraphVertex>> virtualHalfedgeVertices(globalGeometry->mesh, std::vector<knitGraphVertex>{});
    for (knitGraphVertex v : vertices){
        if (!v.isVirtual) continue;
        if (v.edge.value().isBoundary()) continue;
        virtualHalfedgeVertices[v.halfedge.value()].push_back(v);
    }

    for (Edge e : globalGeometry -> mesh.edges()){
        if (e.isBoundary()) continue;
        int ctr = 0;
        std::vector<knitGraphVertex> vertices1 = virtualHalfedgeVertices[e.halfedge()];
        std::vector<knitGraphVertex> vertices2 = virtualHalfedgeVertices[e.halfedge().twin()];
        for (auto v1 : vertices1){
            for (auto v2 : vertices2){
                if (norm(v1.position - v2.position) < eps){
                    ctr++;
                }
            }
        }
        if (ctr == 0){
            std::cout << "ctr = " << ctr << std::endl;
            std::cout << "is edge boundary? " << e.halfedge().isInterior() << std::endl;
            std::cout << "is twin boundary? " << e.halfedge().twin().isInterior() << std::endl;
            std::cout << "-----------------" << std::endl;
        }
    }

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
            //v.hasBeenHandled = true;
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

    std::cout << "size of cluster = " << vCluster.size() << std::endl;
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
                    std::cout << "here for global row in!! " << std::endl;
                    global_row_in = vj.id;
                }
                if (vi.row_out == vj.id){
                    std::cout << "here for global row out!! " << std::endl;
                    global_row_out = vj.id;
                }
                if (vi.col_in[0] == vj.id){
                    std::cout << "here for global col in!! " << std::endl;
                    global_col_in = vj.id;
                }
                if (vi.col_out[0] == vj.id){
                    std::cout << "here for global col out!! " << std::endl;
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

        std::cout << "global row in = " << global_row_in << std::endl;
        std::cout << "global row out = " << global_row_out << std::endl;
        std::cout << "global col in = " << global_col_in << std::endl;
        std::cout << "global col out = " << global_col_out << std::endl;
        std::cout << "----------------------------" << std::endl;

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
    for (knitGraphVertex& v: vertices){
        if (v.isVirtual) continue;//only handle real vertices
        mp.insert({v.id, i++});
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

    // for (auto &v : realVertices){
    //     std::cout << "id = " << v.id << std::endl;
    //     std::cout << "row in = " << v.row_in << std::endl;
    //     std::cout << "row out = " << v.row_out << std::endl;
    //     std::cout << "col_in[0] = " << v.col_in[0] << std::endl;
    //     std::cout << "col_out[0] = " << v.col_out[0] << std::endl;
    //     std::cout << "---------------" << std::endl;
    // }
}


//make obj for yarn-level rendering
void KnitGraph::makeObj(){

    std::vector<Vector3> positions;
    std::vector<std::vector<int>> faces;
    std::vector<std::vector<int>> edges;
    for (knitGraphVertex v : realVertices){
        positions.push_back(v.position);
        //skip short-rows for now
        if (v.row_out == -1) continue;
        //skip increas/decreases for now
        if (v.col_out[0] == -1) continue;
        //this is your standard quad face
        std::vector<int> currFace = std::vector<int>{v.id + 1, v.row_in + 1, realVertices[v.row_in].col_out[0] + 1, v.col_out[0] + 1};
        edges.push_back(std::vector<int>{0, 1, 2, 1});
        faces.push_back(currFace);
    }

    std::ofstream outfile("render.obj");
    
    for (auto p : positions){
        outfile << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }
    for (auto f : faces){
        outfile << "f " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;
    }
    for (auto e : edges){
        outfile << "e " << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << std::endl;            
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
            }
        }
        if (v.col_out[0] != -1){
            if (realVertices[v.col_out[0]].col_in[0] != v.id && realVertices[v.col_out[0]].col_in[1] != v.id){
                buggy_vertices.push_back(v.position);
                num_errors++;
                std::cout << "Column 0 mismatch at vertex " << v.id << std::endl;
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
                std::cout << "Column 0 mismatch at vertex " << v.id << std::endl;
            }
        }

        if (!realVertices[v.col_in[0]].col_out[0] == v.id || realVertices[v.col_in[0]].col_out[1] == v.id){
            buggy_vertices.push_back(v.position);
            num_errors++;
            std::cout << "Column mismatch at vertex (autoknit assertion) " << v.id << std::endl;
        }
    }

    for (auto const &vi : realVertices){
        for (auto &vj : realVertices){
            if (norm(vi.position - vj.position) < eps && vi.id != vj.id){
                std::cout << "2 real vertices are in the same place " << std::endl;
            } 
        }
    }
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

    file.close();
    std::cout << "wrote knit graph to txt file " << std::endl;
}


