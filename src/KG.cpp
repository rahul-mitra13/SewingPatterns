#include "KG.h"


//find the sign of a value
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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
    updateSingularMatchings();
    adjustFaceLevelSets();

    makeAdjustedCourseVirtualVertices();

    //render the graph
    // 1) Collect real vertices and build a pointer->index map
    // std::vector<Vector3> realVs;
    // std::vector<std::array<int, 2>> realEdges;
    // realVs.reserve(allVertices.size());

    // std::map<KGVertex*, int> idxOf;  // which index in realVs each real vertex has

    // for (auto& up : allVertices) {
    //     KGVertex* v = up.get();
    //     if (v->isAlphaVirtual || v->isBetaVirtual) continue; // render only real vertices
    //     // compute 3D position from barycentrics on the corresponding global face
    //     int fIndex = v->halfedge->face().getIndex();
    //     Face f = globalGeometry->mesh.face(fIndex);
    //     Vertex vI = f.halfedge().vertex();
    //     Vertex vJ = f.halfedge().next().vertex();
    //     Vertex vK = f.halfedge().next().next().vertex();
    //     Vector3 pI = globalGeometry->vertexPositions[vI];
    //     Vector3 pJ = globalGeometry->vertexPositions[vJ];
    //     Vector3 pK = globalGeometry->vertexPositions[vK];
    //     v->position = v->baryCoords[0] * pI + v->baryCoords[1] * pJ + v->baryCoords[2] * pK;

    //     int i = static_cast<int>(realVs.size());
    //     idxOf[v] = i;                    // pointer -> node index
    //     realVs.emplace_back(v->position);
    // }

    // // helper to add an edge only if both endpoints are in the node set
    // auto addEdge = [&](KGVertex* a, KGVertex* b) {
    //     auto ia = idxOf.find(a);
    //     auto ib = idxOf.find(b);
    //     if (ia == idxOf.end() || ib == idxOf.end()) return; // skip if an endpoint wasn't rendered
    //     if (ia->second == ib->second) return;               // skip self-loop
    //     realEdges.push_back({ ia->second, ib->second });
    // };

    // // 2) Add edges using the pointer->index map (NOT vertex ids)
    // for (auto& up : allVertices) {
    //     KGVertex* v = up.get();
    //     if (!v) continue;
    //     if (v->isAlphaVirtual || v->isBetaVirtual) continue;

    //     addEdge(v, v->row_out_vertex);
    //     addEdge(v, v->col_out_vertex[0]);
    //     // for increases
    //     // addEdge(v, v->col_out_vertex[1]);
    // }

    // // 3) Register the network
    // polyscope::registerCurveNetwork("Real vertices knit graph", realVs, realEdges);

}

//Makes virtual vertices in the course direction
void KG::makeCourseVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeVirtualVerticesOnBorder(f, true);
    }
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
        for (KGVertex *v : faceKGVertices[f]) {
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
                //assign any halfedge on that face
                v->halfedge = f.halfedge();
                raw->id = vertexID++;
                raw->alpha_tag = j;
                raw->beta_tag = k;
                allVertices.emplace_back(std::move(v));
                faceKGVertices[f].emplace_back(raw);
            }
        }
    }

    if (dot(cross(gradAlpha, gradBeta), n) < 0) {
        for (KGVertex *v : faceKGVertices[f]) {
            v->alpha_tag = -v->alpha_tag;
            v->beta_tag = -v->beta_tag;
        }
    }
}

void KG::makeFaceConnections(){

    double eps = 1e-8;

    auto approx_contains = [eps](const std::vector<double>& xs, double x) {
        for (double v : xs) if (std::fabs(v - x) <= eps) return true;
        return false;
    };

    for (Face f : gluedGeometry->mesh.faces()) {

        // Non-owning reference to this face's vertices
        std::vector<KGVertex*> faceVertices = faceKGVertices[f];
    

        std::vector<double> uniqueAlphas;
        std::vector<double> uniqueBetas;
    
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
                if (std::fabs(v->alpha_tag - currAlphaVal) <= eps) {
                    currAlphaRow[v->beta_tag] = v;
                }
            }
            for (auto it = currAlphaRow.begin(); it != std::prev(currAlphaRow.end()); it++) {
                KGVertex* currVertex = it->second;    
                KGVertex* nextVertex = std::next(it)->second; 
                //update the pointers
                currVertex->row_out_vertex = nextVertex;
                nextVertex->row_in_vertex = currVertex;
            }
        }

        // ---- Connect along wale (columns): for each ~equal beta, order by alpha and link neighbors ----
        for (double currBetaVal : uniqueBetas) {
            std::map<double, KGVertex*> currBetaCol; // key: alpha, val: vertex*
            for (KGVertex* v : faceVertices) {
                if (std::fabs(v->beta_tag - currBetaVal) <= eps) {
                    currBetaCol[v->alpha_tag] = v;
                }
            }
            for (auto it = currBetaCol.begin(); it != std::prev(currBetaCol.end()); it++) {
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

    //before we do anything else 
    //ensure that all real vertices have connections 
    for (auto& up : allVertices) {
        KGVertex* v = up.get();
        if (v->isAlphaVirtual || v->isBetaVirtual) continue;
        ensure(v->row_in_vertex != nullptr && "real vertex doesn't have row_in set");
        ensure(v->row_out_vertex != nullptr && "real vertex doesn't have row_out set");
        ensure(v->col_in_vertex[0] != nullptr && "real vertex doesn't have col_in[0] set");
        ensure(v->col_out_vertex[0] != nullptr && "real vertex doesn't have col_out[0] set");
    }

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

    // Make virtual connections across regular course edges
    for (Edge e : (gluedGeometry->mesh).edges()){ 
    
        if (!e.isBoundary() && courseSingularEdgesGlued[e] == 0) {

            std::vector<KGVertex*> he1CourseVertices = halfedgeCourseVertices[e.halfedge()];
            std::vector<KGVertex*> he2CourseVertices = halfedgeCourseVertices[e.halfedge().twin()];

            //Matchings across regular edges
            std::vector<std::pair<int, int>> regularMatchings;
            for (int i = 0; i < he1CourseVertices.size(); i++) {
                regularMatchings.push_back({i, (he1CourseVertices.size() - i) - 1});
            }

            // Connect the virtual matchings first
            for (auto [i1, i2] : regularMatchings) {
                KGVertex* v1 = he1CourseVertices[i1];
                KGVertex* v2 = he2CourseVertices[i2];
                ensure(v1->isAlphaVirtual && "vertex on halfedge is not virtual");
                ensure(v2->isAlphaVirtual && "vertex on halfege is not virtual");
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
    std::map<KGVertex*, KGVertex*> matchings; // nullptr means short row

    // Order positive edges in decreasing order of time
    std::vector<Edge> orderedPosEdges;
    std::map<int, Edge, std::greater<int>> posEdgesByTime;
    for (Edge edge : (gluedGeometry->mesh).edges()) if (!edge.isBoundary() && courseSingularEdgesGlued[edge] > 0)
        posEdgesByTime[courseSingularEdgesGlued[edge]] = edge;
    for (auto &[time,edge] : posEdgesByTime)
        orderedPosEdges.push_back(edge);

    //loop over positive course edges
    for (Edge startEdge : orderedPosEdges){

        int startEdgeOrder = round(courseSingularEdgesGlued[startEdge]);
        ensure(!startEdge.isBoundary() && "Start edge is not a boundary edge");
        
        std::vector<KGVertex*> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
        std::vector<KGVertex*> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];

        if (he1Vertices.size() > he2Vertices.size())
            swap(he1Vertices, he2Vertices);

        ensure(he2Vertices.size() - he1Vertices.size() == 1 && "More than one stripe born/dying at singular edge");

        bool success = false;
        while(!success){

            // Pick the first unmatched vertex on he2Vertices
            KGVertex* startVertex = nullptr;
            for (KGVertex* v : he2Vertices) {
                if (!matchings.count(v)) {
                    startVertex = v;
                    break;
                }
            }

            matchings[startVertex] = nullptr;
            // Trace short row. Now that we go in the direction of row_in (== right)!
            KGVertex* walker = startVertex;
            ensure(walker->row_in_vertex != nullptr && "startVertex picked doesn't have a row_in_vertex");

            while (true){

                if (walker->row_in_vertex == nullptr){
                    
                    //we've hit a singular edge
                    ensure(walker->isAlphaVirtual && "walker hit a vertex that is not virtual");//walker must be virtual 
                    ensure(!matchings.count(walker) && "walker is already matched");//walker must be unmatched
                    ensure(walker->halfedge.has_value() && "walker's halfedge option has no value");//walker must be on a halfedge

                    Edge edge = walker->halfedge->edge();
                    int edgeOrder = round(courseSingularEdgesGlued[edge]);

                    if (edgeOrder == -startEdgeOrder){
                        //It's a match! We're done 
                        matchings[walker] = nullptr;
                        success = true;
                        break;
                    }else{
                        // We need to cross this edge. Above or below?
                        std::vector<KGVertex*> leftVertices = halfedgeCourseVertices[edge.halfedge()];
                        std::vector<KGVertex*> rightVertices = halfedgeCourseVertices[edge.halfedge().twin()];

                        if (sgn((int)rightVertices.size() - (int)leftVertices.size()) != sgn(edgeOrder))
                            swap(leftVertices, rightVertices);
                        //not sure what this assertion is doing
                        ensure(sgn((int)rightVertices.size() - (int)leftVertices.size()) == sgn(edgeOrder)); // make sure the sides are correct!

                        // Find index along half-edge
                        int indexAlongHalfedge = -1;
                        for (int i = 0; i < leftVertices.size(); i++) {
                            if (leftVertices[i] == walker) {
                                indexAlongHalfedge = i;
                                break;
                            }
                        }

                        KGVertex* connectTo;
                        if (abs(edgeOrder) > startEdgeOrder) {
                            // Go below
                            connectTo = rightVertices[rightVertices.size()-1-indexAlongHalfedge];
                        } else {
                            // Go above (also happens if we looped to the starting edge)
                            connectTo = rightVertices[leftVertices.size()-1-indexAlongHalfedge];
                        }
                        ensure(connectTo->isAlphaVirtual && "connectTo vertex is not alpha virtual");
                        matchings[walker] = connectTo;
                        matchings[connectTo] = walker;
                        walker = connectTo;
                        if (edgeOrder == startEdgeOrder) {
                            // We looped around! Break and continue onto the next short row candidate
                            break;
                        }
                    }
                }else{
                    walker = walker->row_in_vertex;
                }
            }
        }
    }

    // Connect the remaining unmatched virtual vertices across singular edges
    for (Edge startEdge : (gluedGeometry->mesh).edges()){ 
        if (!startEdge.isBoundary() && courseSingularEdgesGlued[startEdge] != 0) {
            std::vector<KGVertex*> he1Vertices = halfedgeCourseVertices[startEdge.halfedge()];
            std::vector<KGVertex*> he2Vertices = halfedgeCourseVertices[startEdge.halfedge().twin()];

            int i2 = he2Vertices.size()-1;
            for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
                if (!matchings.count(he1Vertices[i1])) {
                    while (matchings.count(he2Vertices[i2])) {
                        i2--;
                        if (i2 < 0) {
                            showEdges("i2 < 0", {startEdge}, *globalGeometry);
                            polyscope::show();
                        }
                        assert(i2 >= 0);
                    }
                    matchings[he1Vertices[i1]] = he2Vertices[i2];
                    matchings[he2Vertices[i2]] = he1Vertices[i1];
                    i2--;
                }
            }
        }
    }

    //render the matched vertices in the course direction 
    std::vector<Vector3> courseMatchingsLocations;
    for (auto [v1, v2] : matchings){
        if (v1 != nullptr && v2 != nullptr){
            v1->position = getKGPosition(v1);
            v2->position = getKGPosition(v2);
            courseMatchingsLocations.emplace_back(v1->position);
            courseMatchingsLocations.emplace_back(v2->position);
            courseMatchings.emplace_back(std::make_pair(v1, v2)); 
        }
    }
    polyscope::registerPointCloud("course matchings", courseMatchingsLocations);
    
    
    // Now actually connect all course matchings
    for (auto [i1, i2] : matchings) {
        if (i1 == nullptr || i2 == nullptr) continue;
        KGVertex* v1 = i1;
        KGVertex* v2 = i2;
        if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
            v2->row_out_vertex = v1;
            v1->row_in_vertex = v2;
        }
        if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
            v1->row_out_vertex = v2;
            v2->row_in_vertex = v1;
        }
    }

    std::vector<Vector3> waleMatchingLocations;
    // Connect wale vertices 
    for (Edge e : (gluedGeometry->mesh).edges()) {
        if (e.isBoundary()) continue;//don't need to handle boundary vertices

        std::vector<KGVertex*> he1Vertices = halfedgeWaleVertices[e.halfedge()];
        std::vector<KGVertex*> he2Vertices = halfedgeWaleVertices[e.halfedge().twin()];

        // classical singularity with 1 stripe being born/dying: nothing to do
        if (he1Vertices.size() + he2Vertices.size() == 1)
            continue;

        std::vector<std::pair<int, int>> matchings;

        // more funky singularity: we need to match
        if (he1Vertices.size() != he2Vertices.size()) { 

            // Check which side of the triangles he1 and he2 are located
            // TODO: find a way to do this without epsilons
            int side1, side2;
            for (int i = 0; i < 3; i++) {
                if (abs(he1Vertices[0]->baryCoords[i]) < 1e-8)
                    side1 = (i+1)%3;
                if (abs(he2Vertices[0]->baryCoords[i]) < 1e-8)
                    side2 = (i+1)%3;
            }

            // Fetch coordinate along edge
            std::vector<double> coordsAlongHe1, coordsAlongHe2;
            for (KGVertex *v : he1Vertices)
                coordsAlongHe1.push_back(v->baryCoords[side1]);
            for (KGVertex *v : he2Vertices)
                coordsAlongHe2.push_back(1 - v->baryCoords[side2]); // need to invert to be in the same basis

            // Find out best matchings (greedy approach)
            if (he1Vertices.size() < he2Vertices.size()) { // match every vertex of he1 to closest vertex of he2
                for (int i1 = 0; i1 < he1Vertices.size(); i1++) {
                    int i2closest = -1;
                    for (int i2 = 0; i2 < he2Vertices.size(); i2++)
                        if (i2closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1] - coordsAlongHe2[i2closest]))
                            i2closest = i2;
                    if (i2closest == -1) {
                        std::cout << "Error: No match found for vertex " << he1Vertices[i1]->id << std::endl;
                        polyscope::show();
                    } else {
                        matchings.push_back({i1, i2closest});
                        KGVertex *v1 = he1Vertices[i1];
                        KGVertex *v2 = he2Vertices[i2closest];
                        ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");
                        waleMatchingLocations.emplace_back(getKGPosition(v1));
                        waleMatchingLocations.emplace_back(getKGPosition(v2));
                        waleMatchings.emplace_back(v1, v2);
                    }
                }
            } else { // match every vertex of he2 to closest vertex of he1
                for (int i2 = 0; i2 < he2Vertices.size(); i2++) {
                    int i1closest = -1;
                    for (int i1 = 0; i1 < he1Vertices.size(); i1++)
                        if (i1closest == -1 || abs(coordsAlongHe1[i1] - coordsAlongHe2[i2]) < abs(coordsAlongHe1[i1closest] - coordsAlongHe2[i2]))
                            i1closest = i1;
                    if (i1closest == -1) {
                        std::cout << "Error: No match found for vertex " << he2Vertices[i2]->id << std::endl;
                        polyscope::show();
                    } else {
                        matchings.push_back({i1closest, i2});
                        KGVertex *v1 = he1Vertices[i1closest];
                        KGVertex *v2 = he2Vertices[i2];
                        ensure(v1 != nullptr && v2 != nullptr && "wale matching has nullptr");//a wale matching should always exist
                        waleMatchingLocations.emplace_back(getKGPosition(v1));
                        waleMatchingLocations.emplace_back(getKGPosition(v2));
                        waleMatchings.emplace_back(v1, v2);
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
            KGVertex *v1 = he1Vertices[i1];
            KGVertex *v2 = he2Vertices[i2];
            if (v1->col_in_vertex[0] == nullptr && v2->col_in_vertex[0] != nullptr){
                v1->col_in_vertex[0] = v2;
                v2->col_out_vertex[0] = v1;
            }
            if (v2->col_in_vertex[0] == nullptr && v1->col_in_vertex[0] != nullptr){
                v1->col_out_vertex[0] = v2;
                v2->col_in_vertex[0] = v1;
            }
        }
    }
    polyscope::registerPointCloud("wale matchings", waleMatchingLocations);

    // Now connect real vertices to one another
    for (auto& up : allVertices) {
        KGVertex* v0 = up.get();
        if (v0->isAlphaVirtual || v0->isBetaVirtual) continue; // only real vertices

        // --------Course--------
        {
            KGVertex* v = v0->row_out_vertex;     // start from immediate neighbor
            bool isGluedPath = false;

            while (v && v->isAlphaVirtual) {
                if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true;
                v = v->row_out_vertex;            // step
            }

            if (!v || v->isAlphaVirtual) {
                v0->row_out_vertex = nullptr;     // no real neighbor reachable
            } else {
                v0->row_out_vertex = v;           // connect reciprocally
                v->row_in_vertex   = v0;
                if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
            }
        }

        // --------Wale--------
        {
            KGVertex* v = v0->col_out_vertex[0];
            bool isGluedPath = false;

            while (v && v->isBetaVirtual) {
                if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true;
                v = v->col_out_vertex[0];
            }

            if (!v || v->isBetaVirtual) {
                v0->col_out_vertex[0] = nullptr;
            } else {
                v0->col_out_vertex[0] = v;
                v->col_in_vertex[0]   = v0;
                if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
            }
        }
    }
}

void KG::updateSingularMatchings(){

    auto edgeParam = [&](KGVertex* v) -> double {
        Halfedge he = v->halfedge.value();
        Face f = he.face();
        Halfedge h0 = f.halfedge();
        Halfedge h1 = h0.next();
        Halfedge h2 = h1.next();

        if (he == h0)      return v->baryCoords[1]; // (i->j): t = b_j
        else if (he == h1) return v->baryCoords[2]; // (j->k): t = b_k
        else               return v->baryCoords[0]; // (k->i): t = b_i
    };

    auto setBaryOnEdge = [&](KGVertex* v, double t) {
        Halfedge he = v->halfedge.value();
        Face f = he.face();
        Halfedge h0 = f.halfedge();
        Halfedge h1 = h0.next();
        Halfedge h2 = h1.next();

        if (he == h0)      v->baryCoords = Vector3{1.0 - t, t,         0.0};
        else if (he == h1) v->baryCoords = Vector3{0.0,      1.0 - t,  t  };
        else               v->baryCoords = Vector3{t,        0.0,      1.0 - t};
    };


    auto updateAlphaBetaTags = [&](KGVertex* v) {
        Halfedge he = v->halfedge.value();
        Face f = he.face();
        Halfedge h0 = f.halfedge();
        Halfedge h1 = h0.next();
        Halfedge h2 = h1.next();
        double alphaI = courseOneForm[h0.corner()];
        double alphaJ = courseOneForm[h1.corner()];
        double alphaK = courseOneForm[h2.corner()];
        double betaI = waleOneForm[h0.corner()];
        double betaJ = waleOneForm[h1.corner()];
        double betaK = waleOneForm[h2.corner()];
        v->alpha_tag = v->baryCoords[0] * alphaI + v->baryCoords[1] * alphaJ + v->baryCoords[2] * alphaK;
        v->beta_tag = v->baryCoords[0] * betaI + v->baryCoords[1] * betaJ + v->baryCoords[2] * betaK;

    };

    auto clamp01 = [](double x) { return std::max(0.0, std::min(1.0, x)); };

    std::vector<Vector3> newCoursePos;
    for (auto& [v1, v2] : courseMatchings) {
       
        Halfedge h1 = v1->halfedge.value();
        Halfedge h2 = v2->halfedge.value();
        ensure(h1 == h2.twin())//must be on the same edge

        // Measure both in h1's orientation
        double t1 = edgeParam(v1);
        double t2_raw = edgeParam(v2);
        double t2_in_h1 = 1.0 - t2_raw;

        // Average in the common frame
        double t_avg = clamp01(0.5 * (t1 + t2_in_h1));

        // Write back: each side uses its own local orientation
        setBaryOnEdge(v1, t_avg);
        setBaryOnEdge(v2, (1.0 - t_avg));
        
        //update the alpha/beta tags given the new coordinates 
        updateAlphaBetaTags(v1);
        updateAlphaBetaTags(v2);

        newCoursePos.emplace_back(getKGPosition(v1));
        newCoursePos.emplace_back(getKGPosition(v2));

    }
    polyscope::registerPointCloud("new course positions", newCoursePos);

    std::vector<Vector3> newWalePos;
    for (auto& [v1, v2] : waleMatchings) {
       
        Halfedge h1 = v1->halfedge.value();
        Halfedge h2 = v2->halfedge.value();
        ensure(h1 == h2.twin())//must be on the same edge

        // Measure both in h1's orientation
        double t1 = edgeParam(v1);
        double t2_raw = edgeParam(v2);
        double t2_in_h1 = 1.0 - t2_raw;

        // Average in the common frame
        double t_avg = clamp01(0.5 * (t1 + t2_in_h1));

        // Write back: each side uses its own local orientation
        setBaryOnEdge(v1, t_avg);
        setBaryOnEdge(v2, (1.0 - t_avg));

        //update the alpha/beta tags given the new coordinates 
        updateAlphaBetaTags(v1);
        updateAlphaBetaTags(v2);

        newWalePos.emplace_back(getKGPosition(v1));
        newWalePos.emplace_back(getKGPosition(v2));
    }
    polyscope::registerPointCloud("new wale positions", newWalePos);
}

void KG::adjustFaceLevelSets(){

    for (Face f : gluedGeometry->mesh.faces()){

        std::vector<KGVertex*> faceVertices = faceKGVertices[f];

        for (KGVertex *v : faceVertices){
            if (v->isAlphaVirtual) faceCourseLevelSets[f].emplace_back(v->alpha_tag);
            if (v->isBetaVirtual) faceWaleLevelSets[f].emplace_back(v->beta_tag);
        }
    }
}

void KG::makeAdjustedCourseVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeAdjustedVirtualVerticesOnBorder(f, true);
    }
}

void KG::makeAdjustedWaleVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeAdjustedVirtualVerticesOnBorder(f, false);
    }
}

void KG::makeAdjustedVirtualVerticesOnBorder(Face& f, bool isCourseDirection){

    std::cout << "In here " << std::endl;
}

