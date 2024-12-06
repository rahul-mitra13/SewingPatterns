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

    //merge the virtual vertices 
    //mergeVirtual();

    //reorder vertex indices
    makeRealVertices();
    
    //render the knit graph
    renderGraph();


    //make the obj
    //makeObj();

    std::cout << "Completed building graph..." << std::endl;

}

void KnitGraph::handleCourseNonSingularFaceWaleNonSingularFace(Face &f){
    
    double alphaI, alphaJ, alphaK, alpha_min, alpha_max, alpha_start, alpha_end, betaI, betaJ, betaK, beta_min, beta_max, beta_start, beta_end;

    double eps = 1e-8;//for error tolerance

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
                v.isAlphaVirtual = false;
                v.isBetaVirtual = false;
                faceVertices.push_back(v);
                //set the position from the bary coords
                Vector3 i = globalGeometry->vertexPositions[f.halfedge().corner().vertex()];
                Vector3 j = globalGeometry->vertexPositions[f.halfedge().next().corner().vertex()];
                Vector3 k = globalGeometry->vertexPositions[f.halfedge().next().next().corner().vertex()];
                v.position = v.baryCoords[0] * i + v.baryCoords[1] * j + v.baryCoords[2] * k;
                vertices.push_back(v);
                vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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
            vertexInfoMap[v.id] = v;
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

void KnitGraph::mergeVirtual(){
    double eps = 1e-8;
    int numCourseErrors = 0;
    int numWaleErrors = 0;
    for (knitGraphVertex &v : vertices){
        int courseCtr = 0;
        int waleCtr = 0;
        if (!v.isVirtual) continue;//skip real vertices for now 
        if (v.edge.has_value() && (std::fabs(courseSingularEdges[v.edge.value()]) > eps || std::fabs(waleSingularEdges[v.edge.value()]) > eps
        || v.edge->isBoundary()))  continue;//skip singular and boundary edges for now
        for (knitGraphVertex &k : vertices){
            if (v.isAlphaVirtual){
                // if (v.row_in == -1 && v.row_out == -1){
                //     std::cout << "alpha virtual vertex " << v.id << " has no row in or out set " << std::endl;
                //     std::cout << "row_in = " << v.row_in << std::endl;
                //     std::cout << "row_out = " << v.row_out << std::endl;
                //     std::cout << "col_in[0] = " << v.col_in[0] << std::endl;
                //     std::cout << "col_out[0] = " << v.col_out[0] << std::endl;
                // }
                if (norm(v.position - k.position) < eps && k.isAlphaVirtual && v.id != k.id){
                    if (v.row_out == -1 && k.row_out != -1){
                        std::cout << "v.row_out is unset but k.row_out is set " << std::endl;
                        vertices[v.row_in].row_out = k.row_out;
                        vertices[k.row_out].row_in = v.row_in;
                    }
                    else if (k.row_out == -1 && v.row_out != -1){
                        std::cout << "k.row_out is unset but v.row_out is set " << std::endl;
                        vertices[k.row_in].row_out = v.row_out;
                        vertices[v.row_out].row_in = k.row_in;
                    }
                    else if (v.row_in == -1 && k.row_in != -1){
                        std::cout << "v.row_in is unset but k.row_in is set " << std::endl;
                        vertices[v.row_out].row_in = k.row_in;
                        vertices[k.row_in].row_out = v.row_out;
                    }
                    else if (k.row_in == -1 && v.row_in != -1){
                        std::cout << "k.row_in is unset but v.row_in is set " << std::endl;
                        vertices[k.row_out].row_in = v.row_in;
                        vertices[v.row_in].row_out = k.row_out;
                    }
                    else if (k.row_out == -1 && v.row_out == -1
                        && k.row_in == -1 && v.row_in == -1){
                        std::cout << "Both row in/out unset in merge " << std::endl;
                        numCourseErrors++;
                    }
                    courseCtr++;
                }
            }
            else if (v.isBetaVirtual){
                if (norm(v.position - k.position) < eps && k.isBetaVirtual && v.id != k.id){
                    if (v.col_out[0] == -1 && k.col_out[0] != -1){
                        std::cout << "v.col_out[0] is unset but k.col_out[0] is set " << std::endl;
                        vertices[v.col_in[0]].col_out[0] = k.col_out[0];
                        vertices[k.col_out[0]].col_in[0] = v.col_in[0];
                    }
                    else if (k.col_out[0] == -1 && v.col_out[0] != -1){
                        std::cout << "k.col_out[0] is unset but v.col_out[0] is set " << std::endl;
                        vertices[k.col_in[0]].col_out[0] = v.col_out[0];
                        vertices[v.col_out[0]].col_in[0] = k.col_in[0];
                    }
                    else if (v.col_in[0] == -1 && k.col_in[0] != -1){
                        std::cout << "v.col_in[0] is unset but k.col_in[0] is set " << std::endl;
                        vertices[k.col_in[0]].col_out[0] = v.col_out[0];
                        vertices[v.col_out[0]].col_in[0] = k.col_in[0];
                    }
                    else if (k.col_in[0] == -1 && v.col_in[0] != -1){
                        std::cout << "k.col_in[0] is unset but v.col_in[0] is set " << std::endl;
                        vertices[v.col_in[0]].col_out[0] = k.col_out[0];
                        vertices[k.col_out[0]].col_in[0] = v.col_in[0];
                    }
                    else if (k.col_out[0] == -1 && v.col_out[0] == -1
                        && k.col_in[0] == -1 && v.col_in[0] == -1){
                        std::cout << "Both row in/out unset in merge " << std::endl;
                        numWaleErrors++;
                    }
                    waleCtr++;
                }
            }
        }
        if (courseCtr != 1 && v.isAlphaVirtual){
            std::cout << "number of course vertices being merged at edge " << v.edge.value() << " is " << courseCtr << std::endl;
            std::cout << "Something's gone wrong with the knit graph construction! " << std::endl;
        }
        if (waleCtr != 1 && v.isBetaVirtual){
            std::cout << "number of wale vertices being merged at edge " << v.edge.value() << " is " << waleCtr << std::endl;
            std::cout << "Something's gone wrong with knit graph construction! " << std::endl;
        }
    }
    std::cout << "number of course errors =  " << numCourseErrors << std::endl;
    std::cout << "number of wale errors = " << numWaleErrors << std::endl;
}

void KnitGraph::renderGraph(){
    
    for (auto &v : realVertices){
        if (v.row_out != -1 && realVertices[v.row_out].row_in == v.id){
            vertexPositions.push_back(v.position);
            vertexPositions.push_back(realVertices[v.row_out].position);
        }
        if (v.col_out[0] != -1 && realVertices[v.col_out[0]].col_in[0] == v.id){
            vertexPositions.push_back(v.position);
            vertexPositions.push_back(realVertices[v.col_out[0]].position);
        }

        if (v.col_out[1] != -1 && realVertices[v.col_out[1]].col_in[1] == v.id){
            vertexPositions.push_back(v.position);
            vertexPositions.push_back(realVertices[v.col_out[1]].position);

        }
    }

    for (int i = 0; i < vertexPositions.size(); i+=2){
        edges.push_back({i, i + 1});
    }
    auto graphReal = polyscope::registerCurveNetwork("knit graph ", vertexPositions, edges);
    graphReal -> setRadius(0.001);
    graphReal -> setEnabled(false);

    //visualize the knit graph vertices with the virtual connections
    // for (auto &v : vertices){
    //     vertexPositions.push_back(v.position);
    //     if (v.row_out != -1)
    //         edges.push_back({v.id, v.row_out});
    //     if (v.col_out[0] != -1)
    //         edges.push_back({v.id, v.col_out[0]});
    // }
    // auto graphVirtual = polyscope::registerCurveNetwork("knit graph with virtual connections", vertexPositions, edges);
    // graphVirtual -> setRadius(0.001);
    // graphVirtual -> setEnabled(false);

    
    // std::vector<Vector3> preMergingVertices;
    // for (auto v : vertices){
    //     if (v.isVirtual) continue;
    //     preMergingVertices.push_back(v.position);
    // }

    // auto preMergVerts = polyscope::registerPointCloud("real vertices", preMergingVertices);
    // preMergVerts->setEnabled(false);
}

//perform epsilon merging to 
//merge connections across faces
void KnitGraph::epsilonMerging(){

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
                continue;
            }
            v.hasBeenHandled = true;
            mergeCluster(vCluster, eps);
        }
    }
    //finally, removing any lingering connections to virtual vertices for the real vertices 
    // for (auto &v : vertices){
    //     if (!v.isVirtual){
    //         if (vertices[v.row_out].isVirtual) v.row_out = -1;
    //         if (vertices[v.row_in].isVirtual) v.row_in = -1;
    //         if (vertices[v.col_out[0]].isVirtual) v.col_out[0] = -1;
    //         if (vertices[v.col_in[0]].isVirtual) v.col_in[0] = -1;
    //     }
    // }
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
        // for (auto &vi : vCluster){
        //     vi.row_in = -1;
        //     vi.row_out = -1;
        //     vi.col_in[0] = -1;
        //     vi.col_in[1] = -1;
        //     vi.col_out[0] =  -1;
        //     vi.col_out[1] = -1;
        // }

        knitGraphVertex* g_row_in;
        knitGraphVertex* g_row_out; 
        knitGraphVertex* g_col_in;
        knitGraphVertex* g_col_out;
        if (global_row_in != -1 && global_row_out != -1){
            g_row_in = &vertices[global_row_in];
            g_row_out = &vertices[global_row_out];
            g_row_in->row_out = g_row_out->id;
            g_row_out->row_in = g_row_in->id;
        }
        if (global_col_in != -1 && global_col_out != -1){
            g_col_in = &vertices[global_col_in];
            g_col_out = &vertices[global_col_out];
            g_col_in->col_out[0] = g_col_out->id;
            g_col_out->col_in[0] = g_col_in->id;
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
        //insert real vertices into the map
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

    //update the ids in the knit graph
    for (knitGraphVertex& v : vertices){
        if (v.isVirtual) continue;//only handle real vertices
        auto it = mp.find(v.id);
        id = it -> second;//new id

        if (v.row_in == -1){//row_in is unset
            row_in = -1;
        }
        else{
            it = mp.find(v.row_in);//find the row_in from the original graph
            row_in = it -> second;
        }

        if (v.row_out == -1){//row_out is unset
            row_out = -1;
        }
        else{
            it = mp.find(v.row_out);//find the row out from the original graph
            row_out = it -> second;
        }

        if (v.col_in[0] == -1){//col_in_1 is unset
            col_in_1 = -1;
        }
        else{
            it = mp.find(v.col_in[0]);//find the col_in_1 from the orginal graph
            col_in_1 = it -> second;
        }

        if (v.col_in[1] == -1){//col_in_2 is unset
            col_in_2 = -1;
        }
        else{
            it = mp.find(v.col_in[1]);//find the col_in_1 from the orginal graph
            col_in_2 = it -> second;
        }

        if (v.col_out[0] == -1){//col_out_1 is unset
            col_out_1 = -1;
        }
        else{
            it = mp.find(v.col_out[0]);//find the col_out_1 from the orginal graph
            col_out_1 = it -> second;
        }

        if (v.col_out[1] == -1){//col_out_2 is unset
            col_out_2 = -1;
        }
        else{
            it = mp.find(v.col_out[1]);//find the col_out_1 from the orginal graph
            col_out_2 = it -> second;
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
