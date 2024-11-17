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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
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
            v.isVirtual = true;
            v.isAlphaVirtual = true;
            faceVertices.push_back(v);
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            faceVertices.push_back(v);
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
            v.isVirtual = true;
            v.isBetaVirtual = true;
            //insertSmooth(v);
            faceVertices.push_back(v);
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

    std::vector<std::array<int, 2>> edges;
    std::vector<Vector3> pos;
    //query the vertex positons from the barycentric coordinates
    for (knitGraphVertex v : vertices){
        //if (v.isVirtual) continue;
        Face f = globalGeometry->mesh.face(v.face.getIndex());
        Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
        Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
        Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
        pos.push_back(v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k);
        
        if (v.row_out != -1)//row out is set
        edges.push_back(std::array<int, 2>{v.id, v.row_out});
        if (v.col_out[0] != -1)//col out is set 
        edges.push_back(std::array<int, 2>{v.id, v.col_out[0]});
    }

    auto graph = polyscope::registerCurveNetwork("knit graph vertices", pos, edges);
    graph -> setRadius(0.001);
    graph -> setEnabled(false);

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
