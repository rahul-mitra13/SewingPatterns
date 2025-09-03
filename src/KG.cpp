#include "KG.h"


KG::KG(VertexPositionGeometry& globalGeometry,
       EdgeLengthGeometry& gluedGeometry,
       polyscope::SurfaceMesh& psMesh,
       double coursePeriod, double walePeriod,
       CornerData<double>& courseOneForm,
       EdgeData<double>& courseSingularEdges,
       CornerData<double>& waleOneForm,
       EdgeData<double>& waleSingularEdges,
       std::map<int,int>& globalToGluedEdgeMap)
    : vertexID(0)
    , globalGeometry(&globalGeometry)
    , gluedGeometry(&gluedGeometry)
    , psMesh(&psMesh)
    , coursePeriod(coursePeriod)
    , walePeriod(walePeriod)
    , isGlued(gluedGeometry.mesh, false)      
    , courseOneForm(courseOneForm)
    , courseSingularEdges(courseSingularEdges)
    , courseSingularEdgesGlued(gluedGeometry.mesh)
    , waleOneForm(waleOneForm)
    , waleSingularEdges(waleSingularEdges)
    , waleSingularEdgesGlued(gluedGeometry.mesh) 
    , globalToGluedEdgeMap(&globalToGluedEdgeMap)
    , allVertices()
    , finalVertices()
    , stitchedVertices()
    , faceCourseLevelSets(gluedGeometry.mesh)     
    , faceWaleLevelSets(gluedGeometry.mesh)   
    , faceKGVertices(gluedGeometry.mesh){           

        // Bring everything into the glued setting first
        courseSingularEdgesGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, courseSingularEdges, globalToGluedEdgeMap);
        waleSingularEdgesGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, waleSingularEdges, globalToGluedEdgeMap);

        // Flag edges that were glued together
        isGlued = EdgeData<bool>(gluedGeometry.mesh, false);
        for (auto &[globalEdgeID, gluedEdgeID] : globalToGluedEdgeMap) {
            Edge globalEdge = globalGeometry.mesh.edge(globalEdgeID);
            Edge gluedEdge = gluedGeometry.mesh.edge(gluedEdgeID);
            if (globalEdge.isBoundary() && !gluedEdge.isBoundary())
                isGlued[gluedEdge] = true;
        }    
}

void KG::buildGraph(){
    makeCourseVirtualVertices();
    makeWaleVirtualVertices();
    makeRealVertices();
    makeFaceConnections();
    intrinsicMerge();

    //render the graph with all the vertices 
    std::vector<Vector3> allVs;
    std::vector<std::array<int, 2>> edges;
    for (auto &v : allVertices){
        //now compute the position of v in R^3 
        int fIndex = v->halfedge->face().getIndex();
        //grab the global face 
        Face f = globalGeometry->mesh.face(fIndex);
        //grab the vertices on the face
        Vertex vI = f.halfedge().vertex();
        Vertex vJ = f.halfedge().next().vertex();
        Vertex vK = f.halfedge().next().next().vertex();
        //grab the positions
        Vector3 pI = globalGeometry->vertexPositions[vI];
        Vector3 pJ = globalGeometry->vertexPositions[vJ];
        Vector3 pK = globalGeometry->vertexPositions[vK];
        v->position = v->baryCoords[0] * pI + v->baryCoords[1] * pJ + v->baryCoords[2] * pK;
        allVs.emplace_back(v->position);
        if (v->row_out_vertex != nullptr){
            edges.push_back({v->id, v->row_out_vertex->id});
        }
        if (v->col_out_vertex[0] != nullptr){
            edges.push_back({v->id, v->col_out_vertex[0]->id});
        }
    }
    //polyscope::registerCurveNetwork("All vertices knit graph", allVs, edges);

}

//Makes virtual vertices in the course direction
void KG::makeCourseVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeVirtualVerticesOnBorder(f, true);
    }

    std::vector<Vector3> virtualCourseSingularVertices;
}

//Makes virtual vertices in the wale direction
void KG::makeWaleVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeVirtualVerticesOnBorder(f, false);
    }
}

//Make virtual vertices on the border of all faces
void KG::makeVirtualVerticesOnBorder(Face& f, bool isCourseDirection){
        
    //grab the alpha values
    double alphaI = courseOneForm[f.halfedge().corner()];
    double alphaJ = courseOneForm[f.halfedge().next().corner()];
    double alphaK = courseOneForm[f.halfedge().next().next().corner()];
    double alpha_min = std::min({alphaI, alphaJ, alphaK});
    double alpha_max = std::max({alphaI, alphaJ, alphaK});

    //grab the beta values
    double betaI = waleOneForm[f.halfedge().corner()];
    double betaJ = waleOneForm[f.halfedge().next().corner()];
    double betaK = waleOneForm[f.halfedge().next().next().corner()];
    double beta_min = std::min({betaI, betaJ, betaK});
    double beta_max = std::max({betaI, betaJ, betaK});

    //for floating point tolerance 
    double eps = 1e-12;

    //trace the middle of the stripes
    double alpha_start = (std::ceil((alpha_min - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    double alpha_end = (std::floor((alpha_max - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    double beta_start = (std::ceil((beta_min - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
    double beta_end = (std::floor((beta_max - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;

    double bi, bj, bk, j, k;
    
    //query position information to fix alpha_beta tags 
    int fIndex = f.getIndex();
    //grab the global face 
    Face fGlobal = globalGeometry->mesh.face(fIndex);
    //grab the vertices on the face
    Vertex vI = fGlobal.halfedge().vertex();
    Vertex vJ = fGlobal.halfedge().next().vertex();
    Vertex vK = fGlobal.halfedge().next().next().vertex();
    //grab the positions
    Vector3 pI = globalGeometry->vertexPositions[vI];
    Vector3 pJ = globalGeometry->vertexPositions[vJ];
    Vector3 pK = globalGeometry->vertexPositions[vK];
    Vector3 e1 = pJ - pI; // edge (i,j)
    Vector3 e2 = pK - pI; // edge (i,k)
    Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
    double area = n.norm(); // triangle area
    Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);
    Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);

    if (isCourseDirection){//course direction
        //shift by small epsilon to account for floating point error
        for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//fix alpha
            //ij edge
            bi = (j - alphaJ) / (alphaI - alphaJ);
            bj = 1.0 - bi;
            bk = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();		
                KGVertex* raw = v.get();				
   			    raw->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
                raw->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }

            //jk edge
            bj = (j - alphaK) / (alphaJ - alphaK);
            bk = 1.0 - bj;
            bi = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();	
                KGVertex* raw = v.get();		
   			    raw->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
                raw->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }

            //ki edge
            bi = (j - alphaK) / (alphaI - alphaK);
            bk = 1.0 - bi;
            bj = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();
                KGVertex* raw = v.get();							
   			    raw->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
                raw->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }
        }
    }else{//wale direction

        //shift by small epsilon to account for floating point error
   	    for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){        
            //ij edge
            bi = (k - betaJ) / (betaI - betaJ);
            bj = 1.0 - bi;
            bk = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();	
                KGVertex* raw = v.get();					
   			    raw->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
                raw->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }

            //jk edge
            bj = (k - betaK) / (betaJ - betaK);
            bk = 1.0 - bj;
            bi = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();
                KGVertex* raw = v.get();						
   			    raw->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
                raw->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }

            //ki edge
            bi = (k - betaK) / (betaI - betaK);
            bk = 1.0 - bi;
            bj = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();	
                KGVertex* raw = v.get();					
   			    raw->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
                raw->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }
        }
    }

    if (dot(cross(gradAlpha, gradBeta), n) < 0) {
        for (auto &v : faceKGVertices[f]) {
            v->alpha_tag = -v->alpha_tag;
            v->beta_tag = -v->beta_tag;
        }
    }
}

//make real vertices
void KG::makeRealVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeRealVerticesOnInterior(f);
    }
}

//compute the real vertices on triangle interior
void KG::makeRealVerticesOnInterior(Face& f){


    //grab the alpha values
    double alphaI = courseOneForm[f.halfedge().corner()];
    double alphaJ = courseOneForm[f.halfedge().next().corner()];
    double alphaK = courseOneForm[f.halfedge().next().next().corner()];
    double alpha_min = std::min({alphaI, alphaJ, alphaK});
    double alpha_max = std::max({alphaI, alphaJ, alphaK});

    //grab the beta values
    double betaI = waleOneForm[f.halfedge().corner()];
    double betaJ = waleOneForm[f.halfedge().next().corner()];
    double betaK = waleOneForm[f.halfedge().next().next().corner()];
    double beta_min = std::min({betaI, betaJ, betaK});
    double beta_max = std::max({betaI, betaJ, betaK});

    //for floating point tolerance 
    double eps = 1e-12;

    //trace the middle of the stripes
    double alpha_start = (std::ceil((alpha_min - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    double alpha_end = (std::floor((alpha_max - coursePeriod/4.)/coursePeriod) * coursePeriod) + coursePeriod/4.;
    double beta_start = (std::ceil((beta_min - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;
    double beta_end = (std::floor((beta_max - walePeriod/4.)/walePeriod) * walePeriod) + walePeriod/4.;

    //looping variables and linear system variables
    double RHS_alpha, RHS_beta, detA, a1, b1, c1, a2, b2, c2, j, k, bi, bj, bk = 0;

    //solving the linear system below
    a1 = alphaI - alphaK;
    b1 = alphaJ - alphaK;
    c1 = alphaK;
    a2 = betaI - betaK;
    b2 = betaJ - betaK;
    c2 = betaK;
    
    //query position information to fix alpha_beta tags 
    int fIndex = f.getIndex();
    //grab the global face 
    Face fGlobal = globalGeometry->mesh.face(fIndex);
    //grab the vertices on the face
    Vertex vI = fGlobal.halfedge().vertex();
    Vertex vJ = fGlobal.halfedge().next().vertex();
    Vertex vK = fGlobal.halfedge().next().next().vertex();
    //grab the positions
    Vector3 pI = globalGeometry->vertexPositions[vI];
    Vector3 pJ = globalGeometry->vertexPositions[vJ];
    Vector3 pK = globalGeometry->vertexPositions[vK];
    Vector3 e1 = pJ - pI; // edge (i,j)
    Vector3 e2 = pK - pI; // edge (i,k)
    Vector3 n = cross(e1, e2) / 2; // triangle normal (weighted by its area)
    double area = n.norm(); // triangle area
    Vector3 gradAlpha = ((alphaJ - alphaI) * e2 - (alphaK - alphaI) * e1) / (2.0 * area);
    Vector3 gradBeta = ((betaJ - betaI) * e2 - (betaK - betaI) * e1) / (2.0 * area);

    //for the linear system 
    Eigen::Matrix2f A, A_inv;
    Eigen::Vector2f x, b;

    //basically solving a 2 by 2 linear system over every face
   	//shift by small epsilon to account for floating point error
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
                auto v = std::make_unique<KGVertex>();	
                KGVertex* raw = v.get(); 					
   			    raw->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                raw->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }
        }
    }

    if (dot(cross(gradAlpha, gradBeta), n) < 0) {
        for (auto &v : faceKGVertices[f]) {
            v->alpha_tag = -v->alpha_tag;
            v->beta_tag = -v->beta_tag;
        }
    }
}

void KG::makeFaceConnections(){

    constexpr double eps = 1e-8;

    auto approx_contains = [](const std::vector<double>& xs, double x) {
        for (double v : xs) if (std::fabs(v - x) <= eps) return true;
        return false;
    };

    for (Face f : gluedGeometry->mesh.faces()) {

        // Non-owning reference to this face's vertices
        auto& faceVertices = faceKGVertices[f]; // std::vector<KGVertex*>
    

        std::vector<double> uniqueAlphas;
        std::vector<double> uniqueBetas;
        uniqueAlphas.reserve(faceVertices.size());
        uniqueBetas.reserve(faceVertices.size());

        // Collect unique alpha/beta tags using epsilon equality
        for (KGVertex* v : faceVertices) {
            if (!v->isBetaVirtual && !approx_contains(uniqueAlphas, v->alpha_tag)) {
                uniqueAlphas.push_back(v->alpha_tag);
            }
            if (!v->isAlphaVirtual && !approx_contains(uniqueBetas, v->beta_tag)) {
                uniqueBetas.push_back(v->beta_tag);
            }
        }

        // ---- Connect along course (rows): for each ~equal alpha, order by beta and link neighbors ----
        for (double currAlphaVal : uniqueAlphas) {
            std::map<double, KGVertex*> currAlphaRow; // key: beta, val: vertex*
            for (KGVertex* v : faceVertices) {
                if (!v->isBetaVirtual && std::fabs(v->alpha_tag - currAlphaVal) <= eps) {
                    currAlphaRow[v->beta_tag] = v;
                }
            }
            for (auto it = currAlphaRow.begin(); std::next(it) != currAlphaRow.end(); ++it) {
                KGVertex* currVertex = it->second;                // already KGVertex*
                KGVertex* nextVertex = std::next(it)->second;     // KGVertex*
                //update the pointers
                currVertex->row_out_vertex = nextVertex;
                nextVertex->row_in_vertex = currVertex;
            }
        }

        // ---- Connect along wale (columns): for each ~equal beta, order by alpha and link neighbors ----
        for (double currBetaVal : uniqueBetas) {
            std::map<double, KGVertex*> currBetaCol; // key: alpha, val: vertex*
            for (KGVertex* v : faceVertices) {
                if (!v->isAlphaVirtual && std::fabs(v->beta_tag - currBetaVal) <= eps) {
                    currBetaCol[v->alpha_tag] = v;
                }
            }
            for (auto it = currBetaCol.begin(); std::next(it) != currBetaCol.end(); ++it) {
                KGVertex* currVertex = it->second;
                KGVertex* nextVertex = std::next(it)->second;
                //update the pointers 
                currVertex->col_out_vertex[0] = nextVertex;
                nextVertex->col_in_vertex[0] = currVertex;
                
            }
        }
    }
}


//intrinsic merge of the virtual vertices
//this will give us the matchings across singular edges
void KG::intrinsicMerge(){

    // Store these vertices in order of the "direction of the halfedge"
    std::map<Halfedge, std::vector<KGVertex*>> halfedgeCourseVertices;
    std::map<Halfedge, std::vector<KGVertex*>> halfedgeWaleVertices;
    for (const auto& v : allVertices){
        if (v->isAlphaVirtual) halfedgeCourseVertices[v->halfedge.value()].push_back(v.get());
        if (v->isBetaVirtual) halfedgeWaleVertices[v->halfedge.value()].push_back(v.get());
    }

    // Compute a consistent 1D parameter t in [0,1] along the oriented halfedge
    // Face’s three halfedges correspond to edges: ij (h0), jk (h1), ki (h2).
    // Parameter choice consistent with your construction:
    //   ij : t = b_j
    //   jk : t = b_k
    //   ki : t = b_i   (note: uses the face's orientation)
    auto edgeParam = [](KGVertex* v) -> double {
        Halfedge he = v->halfedge.value(); // safe: we guarded above
        Face f = he.face();
        Halfedge h0 = f.halfedge();
        Halfedge h1 = h0.next();
        Halfedge h2 = h1.next();

        if (he == h0)      return v->baryCoords[1]; // ij
        else if (he == h1) return v->baryCoords[2]; // jk
        else               return v->baryCoords[0]; // ki
    };

    auto sortByParam = [&](std::vector<KGVertex*>& vec) {
        std::sort(vec.begin(), vec.end(),
                  [&](KGVertex* a, KGVertex* b) {
                      return edgeParam(a) < edgeParam(b); // ascending along edge
                  });
    };

    // Sort each bucket along the halfedge
    for (auto& [hid, vec] : halfedgeCourseVertices) sortByParam(vec);
    for (auto& [hid, vec] : halfedgeWaleVertices)   sortByParam(vec);

    // Make connections across regular course edges
    for (Edge e : (gluedGeometry->mesh).edges()){ 
    
        if (!e.isBoundary() && courseSingularEdgesGlued[e] == 0) {

            std::vector<KGVertex*> he1CourseVertices = halfedgeCourseVertices[e.halfedge()];
            std::vector<KGVertex*> he2CourseVertices = halfedgeCourseVertices[e.halfedge().twin()];

            // Trivial matchings
            std::vector<std::pair<int, int>> matchings;
            for (int i = 0; i < he1CourseVertices.size(); i++) {
                matchings.push_back({i, (he1CourseVertices.size() - i) - 1});
            }

            // Connect the virtual matchings first
            for (auto [i1, i2] : matchings) {
                KGVertex* v1 = he1CourseVertices[i1];
                KGVertex* v2 = he2CourseVertices[i2];

                if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
                    v2->row_out_vertex = v1;
                    v1->row_in_vertex = v2;
                }
                if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
                    v1->row_out_vertex = v2;
                    v2->row_in_vertex = v1;                   
                }
            }
        }
    }

    // Now we need to connect across course singular edges.
    // For each positive sing, we choose an arbitrary stripe to start a short row (TODO: make a better choice here).
    // We then propagate it until reaching another sing edge. If that sing edge matches the first one, we're done.
    // If not, we connect it across in a way that respects the ordering.
    std::map<int, int> matchings; // -1 means short row

    // Order positive edges in decreasing order of time
    std::vector<Edge> orderedPosEdges;
    std::map<int, Edge, std::greater<int>> posEdgesByTime;
    for (Edge edge : (gluedGeometry->mesh).edges()) if (!edge.isBoundary() && courseSingularEdgesGlued[edge] > 0)
        posEdgesByTime[courseSingularEdgesGlued[edge]] = edge;
    for (auto &[time,edge] : posEdgesByTime)
        orderedPosEdges.push_back(edge);

}
