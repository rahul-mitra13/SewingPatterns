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

    renderGraph();

}

void KnitGraph::handleCourseNonSingularFaceWaleNonSingularFace(Face &f){
    
    double alphaI, alphaJ, alphaK, alpha_min, alpha_max, alpha_start, alpha_end, betaI, betaJ, betaK, beta_min, beta_max, beta_start, beta_end;

    double eps = 1e-8 * period;//for error tolerance

    std::vector<knitGraphVertex> faceVertices;

    //looping stuff and linear system solving stuff
    double RHS_alpha, RHS_beta, detA, a1, b1, c1, a2, b2, c2, j, k, bi, bj, bk;

    //grab the vertices on the face
    Vertex vI = f.halfedge().vertex();
    Vertex vJ = f.halfedge().next().vertex();
    Vertex vK = f.halfedge().next().next().vertex();

    //grab the alpha values
    alphaI = courseOneForm[f.halfedge().corner()];
    alphaJ = courseOneForm[f.halfedge().next().corner()];
    alphaK = courseOneForm[f.halfedge().next().next().corner()];
    alpha_min = std::min({alphaI, alphaJ, alphaK});
    alpha_max = std::max({alphaI, alphaJ, alphaK});
    a1 = alphaI - alphaK;
    b1 = alphaJ - alphaK;
    c1 = alphaK;

    //grab the beta values
    betaI = waleOneForm[f.halfedge().corner()];
    betaJ = waleOneForm[f.halfedge().next().corner()];
    betaK = waleOneForm[f.halfedge().next().next().corner()];
    beta_min = std::min({betaI, betaJ, betaK});
    beta_max = std::max({betaI, betaJ, betaK});
    a2 = betaI - betaK;
    b2 = betaJ - betaK;
    c2 = betaK;

    //trace the middle of the stripes to avoid floating point error
    
    alpha_start = (std::ceil((alpha_min - period/4.)/period) * period) + period/4.;
    alpha_end = (std::floor((alpha_max - period/4.)/period) * period) + period/4.;

    beta_start = (std::ceil((beta_min - period/4.)/period) * period) + period/4.;
    beta_end = (std::floor((beta_max - period/4.)/period) * period) + period/4.;

    Eigen::Matrix2f A, A_inv;
    Eigen::Vector2f x, b;

    //basically solving a 2 by 2 linear system over every face
    //create real vertices
   	//shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (j = alpha_start; j < alpha_end + eps; j += period){//step alpha
   	    for (k = beta_start; k < beta_end + eps; k += period){//step beta
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
                faceVertices.push_back(v);
                //set the position from the bary coords
                Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
                Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
                Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
                v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
                vertices.push_back(v);
            }
        }
    }

    //interpolate the betas
    //shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (j = alpha_start; j < alpha_end + eps; j += period){//fix alpha
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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
        }
    }

    //interpolate the alphas
    //shift up by some small epsilon so that we don't break when Z_start = Z_end
   	for (k = beta_start; k < beta_end + eps; k += period){        
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
            //set the position from the bary coords
            Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
            Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
            Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
            v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
            vertices.push_back(v);
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

    std::set<int> uniqueAlphas;
    std::set<int> uniqueBetas;

    for (knitGraphVertex v : faceVertices){
        if (!v.isBetaVirtual) uniqueAlphas.insert(hashFloatQuantized(v.alpha_tag));
        if (!v.isAlphaVirtual) uniqueBetas.insert(hashFloatQuantized(v.beta_tag));
    }
    for (int currHashedAlphaVal : uniqueAlphas){
        //make an ordered map to store beta values for the current alpha value
        std::map<int, int> currAlphaRow;
        for (knitGraphVertex v : faceVertices){
            if (hashFloatQuantized(v.alpha_tag) == currHashedAlphaVal){
                //if a vertex has the same alpha_tag as the current alpha val, place them in currAlphaRow and order them according to their beta
                currAlphaRow[hashFloatQuantized(v.beta_tag)] = v.id;
            }
        }
        //set the course connections using the map (skip the last element)
        for (auto it = currAlphaRow.begin(); it != std::prev(currAlphaRow.end()); it++){
            knitGraphVertex& currVertex = vertices[it->second];
            knitGraphVertex& nextVertex = vertices[std::next(it)->second]; 

            currVertex.row_out = nextVertex.id;
            nextVertex.row_in = currVertex.id;
        }
    }


    for (int currHashedBetaVal : uniqueBetas){
        //make an ordered map to store alpha values for the current beta value
        std::map<int, int> currBetaCol;
        for (knitGraphVertex v : faceVertices){
            if (hashFloatQuantized(v.beta_tag) == currHashedBetaVal){
                //if a vertex has the same beta_tag as the current beta val, place them in currBetaCol and order them according to their alpha
                currBetaCol[hashFloatQuantized(v.alpha_tag)] = v.id;
            }
        }
        //set the wale connections using the map
        for (auto it = currBetaCol.begin(); it != std::prev(currBetaCol.end()); it++){
            knitGraphVertex& currVertex = vertices[it->second];
            knitGraphVertex& nextVertex = vertices[std::next(it)->second];

            currVertex.col_out[0] = nextVertex.id;
            nextVertex.col_in[0] = currVertex.id;
        }
    }

}

void KnitGraph::renderGraph(){

    std::vector<std::array<int, 2>> edgesVirtual;
    std::vector<Vector3> positionsVirtual;
    //query the vertex positons from the barycentric coordinates
    for (knitGraphVertex &v : vertices){
        //if (v.isVirtual) continue;
        positionsVirtual.push_back(v.position);
        
        if (v.row_out != -1)//row out is set
        edgesVirtual.push_back(std::array<int, 2>{v.id, v.row_out});
        if (v.col_out[0] != -1)//col out is set 
        edgesVirtual.push_back(std::array<int, 2>{v.id, v.col_out[0]});
    }

    auto graph = polyscope::registerCurveNetwork("Virtual knit graph vertices", positionsVirtual, edgesVirtual);
    graph -> setRadius(0.001);
    graph -> setEnabled(false);

    std::vector<Vector3> positionsReal;
    std::vector<std::array<int,2>> edgesReal;
    std::vector<int> faces;

    //gonna be a lot of repetition but it's okay I guess
    for (auto &v : vertices){
        if (v.isVirtual) continue;
        if (v.row_out != -1 && vertices[v.row_out].row_in == v.id && !vertices[v.row_out].isVirtual){
            positionsReal.push_back(v.position);
            positionsReal.push_back(vertices[v.row_out].position);
            
        }
        if (v.col_out[0] != -1 && vertices[v.col_out[0]].col_in[0] == v.id && !vertices[v.col_out[0]].isVirtual){
            positionsReal.push_back(v.position);
            positionsReal.push_back(vertices[v.col_out[0]].position);
           
        }
    }
    for (int i = 0; i < positionsReal.size(); i+=2){
        edgesReal.push_back({i, i + 1});
    }

    polyscope::registerCurveNetwork("Real knit graph vertices ", positionsReal, edgesReal);
}

//perform epsilon merging to 
//merge connections across faces
void KnitGraph::epsilonMerging(){

    //epsilon ball around vertex we hope to merge
    const double eps = 1e-5 * period;
    for (knitGraphVertex &v : vertices) {
        //skip singular edges in the intial merging
        if (v.edge.has_value() && (std::fabs(courseSingularEdges[v.edge.value()]) > 0 || std::fabs(waleSingularEdges[v.edge.value()]) > 0)){ 
            std::cout << "skipping singular edge " << v.edge.value() << std::endl;
            continue;
        }
        // find clusters for all vertices that have not been handled
        if (!v.hasBeenHandled) {
            auto vCluster = findCluster(v, eps);
            if (vCluster.size() == 1) {
                continue;
            }
            v.hasBeenHandled = true;
            mergeCluster(vCluster, eps);
        }
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
        for (auto &vi : vCluster){
            if (vi.id != id_to_keep){
                vi.isVirtual = true;
            }
        }

        knitGraphVertex& realVertex = vertices[id_to_keep];

        realVertex.row_in = global_row_in;
        vertices[global_row_in].row_out = realVertex.id;

        realVertex.row_out = global_row_out;
        vertices[global_row_out].row_in = realVertex.id;

        realVertex.col_in[0] = global_col_in;
        vertices[global_col_in].col_out[0] = realVertex.id;

        realVertex.col_out[0] = global_col_out;
        vertices[global_col_out].col_in[0] = realVertex.id;
    }

    else{
        
        int global_row_in = -1;
        int global_row_out = -1;
        int global_col_in = -1;
        int global_col_out = -1;

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

        //unset the connections for the virtual vertices in the cluster
        for (auto &vi : vCluster){
            vi.row_in = -1;
            vi.row_out = -1;
            vi.col_in[0] = -1;
            vi.col_in[1] = -1;
            vi.col_out[0] =  -1;
            vi.col_out[1] = -1;
        }

        knitGraphVertex& g_row_in = vertices[global_row_in];
        knitGraphVertex& g_row_out = vertices[global_row_out];
        knitGraphVertex& g_col_in = vertices[global_col_in];
        knitGraphVertex& g_col_out = vertices[global_col_out];

        g_row_in.row_out = g_row_out.id;
        g_row_out.row_in = g_row_in.id;
        g_col_in.col_out[0] = g_col_out.id;
        g_col_out.col_in[0] = g_col_in.id;

    }

}

//make obj for yarn-level rendering
void KnitGraph::makeObj(){

    std::cout << "In making obj function " << std::endl;

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
