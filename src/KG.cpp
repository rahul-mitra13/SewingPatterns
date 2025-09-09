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
    , stitchedVertices()
    , faceKGVertices(gluedGeometry.mesh)
    , adjustedVertices()
    , adjustedFaceKGVertices(gluedGeometry.mesh)
    , courseLineSegPairs(gluedGeometry.mesh)
    , waleLineSegPairs(gluedGeometry.mesh)
      {        

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

    //initial knit graph construction
    makeCourseVirtualVertices();
    makeWaleVirtualVertices();
    makeRealVertices();
    makeFaceConnections();
    intrinsicMerge();

    //store matching information 
    findLineSegmentPairs();
    updateSingularMatchings();

    //new knit graph construction
    makeAdjustedCourseVirtualVertices();
    makeAdjustedWaleVirtualVertices();
    makeAdjustedRealVertices();
    makeAdjustedFaceConnections();
    adjustedIntrinsicMerge();

    //visualize the new (after adjustment) virtual connections to make sure we're correct 
    //debugging
    // std::vector<Vector3> Vs;
    // std::vector<std::array<int, 2>> Es;
    // for (Face f : gluedGeometry->mesh.faces()){
    //     std::vector<std::pair<KGVertex*, KGVertex*>> coursePs = courseLineSegPairs[f];
    //     for (auto &[v1, v2] : coursePs){
    //         int i0 = static_cast<int>(Vs.size());
    //         Vs.emplace_back(getKGPosition(v1));
    //         Vs.emplace_back(getKGPosition(v2));
    //         Es.push_back({i0, i0 + 1});
    //     }
    //     std::vector<std::pair<KGVertex*, KGVertex*>> walePs = waleLineSegPairs[f];
    //     for (auto &[v1, v2] : walePs){
    //         int i0 = static_cast<int>(Vs.size());
    //         Vs.emplace_back(getKGPosition(v1));
    //         Vs.emplace_back(getKGPosition(v2));
    //         Es.push_back({i0, i0 + 1});
    //     }
    // }
    // polyscope::registerCurveNetwork("New virtual connections", Vs, Es);

    //render the graph
    //1) Collect real vertices and build a pointer->index map
    std::vector<Vector3> realVs;
    std::vector<std::array<int, 2>> realEdges;
    realVs.reserve(allVertices.size());

    std::map<KGVertex*, int> idxOf;  // which index in realVs each real vertex has

    for (auto& up : adjustedVertices) {
        KGVertex* v = up.get();
        //if (v->isAlphaVirtual || v->isBetaVirtual) continue; // render only real vertices
        // compute 3D position from barycentrics on the corresponding global face
        int fIndex = v->halfedge->face().getIndex();
        Face f = globalGeometry->mesh.face(fIndex);
        Vertex vI = f.halfedge().vertex();
        Vertex vJ = f.halfedge().next().vertex();
        Vertex vK = f.halfedge().next().next().vertex();
        Vector3 pI = globalGeometry->vertexPositions[vI];
        Vector3 pJ = globalGeometry->vertexPositions[vJ];
        Vector3 pK = globalGeometry->vertexPositions[vK];
        v->position = v->baryCoords[0] * pI + v->baryCoords[1] * pJ + v->baryCoords[2] * pK;

        int i = static_cast<int>(realVs.size());
        idxOf[v] = i;                    // pointer -> node index
        realVs.emplace_back(v->position);
    }

    // helper to add an edge only if both endpoints are in the node set
    auto addEdge = [&](KGVertex* a, KGVertex* b) {
        auto ia = idxOf.find(a);
        auto ib = idxOf.find(b);
        if (ia == idxOf.end() || ib == idxOf.end()) return; // skip if an endpoint wasn't rendered
        if (ia->second == ib->second) return;               // skip self-loop
        realEdges.push_back({ ia->second, ib->second });
    };

    // 2) Add edges using the pointer->index map (NOT vertex ids)
    for (auto& up : adjustedVertices) {
        KGVertex* v = up.get();
        if (!v) continue;
        //if (v->isAlphaVirtual || v->isBetaVirtual) continue;
        addEdge(v, v->row_out_vertex);
        addEdge(v, v -> col_out_vertex[0]);
        // for increases
        // addEdge(v, v->col_out_vertex[1]);
    }

    // 3) Register the network
    polyscope::registerCurveNetwork("Adjusted vertices knit graph", realVs, realEdges);

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

    //for solving the linear system below
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
    //polyscope::registerPointCloud("course matchings", courseMatchingsLocations);
    
    
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
    //polyscope::registerPointCloud("wale matchings", waleMatchingLocations);

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

    double eps = 1e-12;
    // for (Face f : gluedGeometry->mesh.faces()) {
    //     auto& F = faceKGVertices[f];
    //     for (KGVertex* v : F) {
    //         if (v->row_out_vertex) {
    //             ensure(v->beta_tag + eps <= v->row_out_vertex->beta_tag);
    //         }
    //         if (v->col_out_vertex[0]) {
    //             ensure(v->alpha_tag + eps <= v->col_out_vertex[0]->alpha_tag);
    //         }
    //     }
    // }
}


void KG::findLineSegmentPairs(){

    double eps = 1e-12;
    struct PairHash {
        size_t operator()(const std::pair<int,int>& p) const noexcept {
            // order-independent: (min, max)
            return std::hash<long long>{}((static_cast<long long>(p.first) << 32) ^ p.second);
        }
    };

    for (Face f : gluedGeometry->mesh.faces()) {

        std::vector<KGVertex*> faceVertices = faceKGVertices[f];

        std::unordered_set<std::pair<int,int>, PairHash> seenCourse, seenWale;

        // course
        for (KGVertex* v1 : faceVertices) {
            if (!v1->isAlphaVirtual) continue;
            double a1 = v1->alpha_tag;
            int cnt = 0;
            for (KGVertex* v2 : faceVertices) {
                if (!v2->isAlphaVirtual || v1 == v2) continue;
                if (std::fabs(a1 - v2->alpha_tag) < eps) {
                    auto key = std::minmax(v1->id, v2->id);
                    if (seenCourse.insert({key.first, key.second}).second) {
                        courseLineSegPairs[f].emplace_back(v1, v2);
                    }
                    cnt++;
                }
            }
            assert(cnt == 1 && "more than 2 virtual vertices with the same alpha tag");
        }

        // wale (analogous)
        for (KGVertex* v1 : faceVertices) {
            if (!v1->isBetaVirtual) continue;
            double b1 = v1->beta_tag;
            int cnt = 0;
            for (KGVertex* v2 : faceVertices) {
                if (!v2->isBetaVirtual || v1 == v2) continue;
                if (std::fabs(b1 - v2->beta_tag) < eps) {
                    auto key = std::minmax(v1->id, v2->id);
                    if (seenWale.insert({key.first, key.second}).second) {
                        waleLineSegPairs[f].emplace_back(v1, v2);
                    }
                    cnt++;
                }
            }
            assert(cnt == 1 && "more than 2 virtual vertices with the same beta tag");
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

    //view the pairs 
    std::vector<Vector3> coursePPos;
    std::vector<Vector3> walePPos;
    for (Face f : gluedGeometry->mesh.faces()){
        std::vector<std::pair<KGVertex*, KGVertex*>> coursePs = courseLineSegPairs[f];
        std::vector<std::pair<KGVertex*, KGVertex*>> walePs = waleLineSegPairs[f];

        for (auto &[v1, v2] : coursePs){
            coursePPos.emplace_back(getKGPosition(v1));
            coursePPos.emplace_back(getKGPosition(v2));
        }

        for (auto &[v1, v2] : walePs){
            walePPos.emplace_back(getKGPosition(v1));
            walePPos.emplace_back(getKGPosition(v2));
        }
    }

    // polyscope::registerPointCloud("course line seg pairs", coursePPos);
    // polyscope::registerPointCloud("wale line seg pairs", walePPos);


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

    // grab the alpha values
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

    std::vector<KGVertex*> faceVertices = faceKGVertices[f];
    if (isCourseDirection){
        std::vector<std::pair<KGVertex*, KGVertex*>> coursePs = courseLineSegPairs[f];
        for (auto &[v1, v2] : coursePs){
            //first vertex in pair
            auto v1New = std::make_unique<KGVertex>();		
            KGVertex* raw1 = v1New.get();				
   			raw1->baryCoords = v1->baryCoords;
            raw1->id = vertexID++;
            raw1->alpha_tag = v1->alpha_tag;
            raw1->beta_tag = v1->beta_tag;
            //edges and halfedges are in the glued mesh setting 
            raw1->halfedge = v1->halfedge;
            raw1->isAlphaVirtual = true;
            adjustedVertices.emplace_back(std::move(v1New));
            adjustedFaceKGVertices[f].emplace_back(raw1);

            //second vertex in pair
            auto v2New = std::make_unique<KGVertex>();		
            KGVertex* raw2 = v2New.get();				
   			raw2->baryCoords = v2->baryCoords;
            raw2->id = vertexID++;
            raw2->alpha_tag = v2->alpha_tag;
            raw2->beta_tag = v2->beta_tag;
            //edges and halfedges are in the glued mesh setting 
            raw2->halfedge = v2->halfedge;
            raw2->isAlphaVirtual = true;
            adjustedVertices.emplace_back(std::move(v2New));
            adjustedFaceKGVertices[f].emplace_back(raw2);
        }
    }
    else{
        std::vector<std::pair<KGVertex*, KGVertex*>> walePs = waleLineSegPairs[f];
        for (auto &[v1, v2] : walePs){
            //first vertex in pair
            auto v1New = std::make_unique<KGVertex>();		
            KGVertex* raw1 = v1New.get();				
   			raw1->baryCoords = v1->baryCoords;
            raw1->id = vertexID++;
            raw1->alpha_tag = v1->alpha_tag;
            raw1->beta_tag = v1->beta_tag;
            //edges and halfedges are in the glued mesh setting 
            raw1->halfedge = v1->halfedge;
            raw1->isBetaVirtual = true;
            adjustedVertices.emplace_back(std::move(v1New));
            adjustedFaceKGVertices[f].emplace_back(raw1);

            //second vertex in pair
            auto v2New = std::make_unique<KGVertex>();		
            KGVertex* raw2 = v2New.get();				
   		    raw2->baryCoords = v2->baryCoords;
            raw2->id = vertexID++;
            raw2->alpha_tag = v2->alpha_tag;
            raw2->beta_tag = v2->beta_tag;
            //edges and halfedges are in the glued mesh setting 
            raw2->halfedge = v2->halfedge;
            raw2->isBetaVirtual = true;
            adjustedVertices.emplace_back(std::move(v2New));
            adjustedFaceKGVertices[f].emplace_back(raw2);
        }
    }

    if (dot(cross(gradAlpha, gradBeta), n) < 0) {
        for (KGVertex *v : adjustedFaceKGVertices[f]) {
            v->alpha_tag = -v->alpha_tag;
            v->beta_tag = -v->beta_tag;
        }
    }
}

void KG::makeAdjustedRealVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeAdjustedRealVerticesOnInterior(f);
    }

}

void KG::makeAdjustedRealVerticesOnInterior(Face& f){

    const double eps = 1e-12;

    // Orienters: course by beta, wale by alpha (no tie-breaks)
    auto orientCourse = [&](KGVertex*& p0, KGVertex*& p1) {
        if (p1->beta_tag < p0->beta_tag) std::swap(p0, p1);
    };
    auto orientWale = [&](KGVertex*& p0, KGVertex*& p1) {
        if (p1->alpha_tag < p0->alpha_tag) std::swap(p0, p1);
    };

    // Solve (a0 + u*(a1-a0)) = (b0 + v*(b1-b0)) in barycentrics via 2x2 systems
    auto solveUV = [&](const Vector3& a0, const Vector3& a1,
                       const Vector3& b0, const Vector3& b1,
                       double& u, double& v) -> bool {
        Vector3 A = a1 - a0, B = b1 - b0;
        auto try2 = [&](int i, int j) -> bool {
            double m00 = A[i], m01 = -B[i];
            double m10 = A[j], m11 = -B[j];
            double r0  = b0[i] - a0[i], r1 = b0[j] - a0[j];
            double det = m00 * m11 - m01 * m10;
            if (std::abs(det) < 1e-14) return false;
            u = ( r0 * m11 - m01 * r1) / det;
            v = ( m00 * r1 - r0  * m10) / det;
            return true;
        };
        return try2(0,1) || try2(1,2) || try2(0,2);
    };

    auto strictlyInside = [&](const Vector3& b) -> bool {
        return (b[0] > eps && b[1] > eps && b[2] > eps);
    };

    auto setTagsOnFace = [&](KGVertex* v, Face f) {
        Halfedge h0 = f.halfedge();
        Halfedge h1 = h0.next();
        Halfedge h2 = h1.next();
        double aI = courseOneForm[h0.corner()];
        double aJ = courseOneForm[h1.corner()];
        double aK = courseOneForm[h2.corner()];
        double bI = waleOneForm[h0.corner()];
        double bJ = waleOneForm[h1.corner()];
        double bK = waleOneForm[h2.corner()];
        v->alpha_tag = v->baryCoords[0] * aI + v->baryCoords[1] * aJ + v->baryCoords[2] * aK;
        v->beta_tag  = v->baryCoords[0] * bI + v->baryCoords[1] * bJ + v->baryCoords[2] * bK;
    };

    
    const auto& coursePairs = courseLineSegPairs[f];
    const auto& walePairs   = waleLineSegPairs[f];

    for (auto cp : coursePairs) {
        KGVertex *c0 = cp.first, *c1 = cp.second;
        orientCourse(c0, c1); // beta increases
        const Vector3 a0 = c0->baryCoords, a1 = c1->baryCoords;

        for (auto wp : walePairs) {
            KGVertex *w0 = wp.first, *w1 = wp.second;
            orientWale(w0, w1); // alpha increases
            const Vector3 b0 = w0->baryCoords, b1 = w1->baryCoords;

            double u, v;
            if (!solveUV(a0, a1, b0, b1, u, v)) continue;
            if (u < -eps || u > 1 + eps || v < -eps || v > 1 + eps) continue;

            Vector3 bary = (1.0 - u) * a0 + u * a1;
            double s = bary[0] + bary[1] + bary[2];
            if (std::abs(s - 1.0) > 1e-9) bary = (1.0 / s) * bary;

            if (!strictlyInside(bary)) continue;

            auto nv = std::make_unique<KGVertex>();
            KGVertex* raw = nv.get();
            raw->baryCoords = bary;
            raw->id         = vertexID++;
            raw->halfedge   = f.halfedge(); // any halfedge on this face
            setTagsOnFace(raw, f);

            adjustedFaceKGVertices[f].push_back(raw);
            adjustedVertices.emplace_back(std::move(nv));
        }
    }
}


//there has to be a better way of writing this
void KG::makeAdjustedFaceConnections(){

    const double tol_col = 1e-10;   // colinearity/point-on-line tolerance in barycentric space
    const double tol_t   = 1e-10;   // param range tolerance

    // Return true and t if x lies on segment p0->p1 in barycentric space.
    auto onSegmentParam = [&](const Vector3& x, const Vector3& p0, const Vector3& p1, double& t) -> bool {
        Vector3 A = p1 - p0;
        double AA = dot(A, A);
        if (AA < 1e-18) return false;                  // degenerate segment
        Vector3 r = x - p0;
        t = dot(r, A) / AA;                            // least-squares param along A
        if (t < -tol_t || t > 1.0 + tol_t) return false;
        Vector3 proj = p0 + t * A;                     // closest point on the line
        return (proj - x).norm() <= tol_col;           // near the segment line
    };

    for (Face f : gluedGeometry->mesh.faces()) {
        auto& F = adjustedFaceKGVertices[f];           // all adjusted verts (endpoints + intersections) on this face

        // reset connections on this face
        for (KGVertex* v : F) {
            v->row_in_vertex = v->row_out_vertex = nullptr;
            v->col_in_vertex[0] = v->col_out_vertex[0] = nullptr;
        }

        // ----------------- Course: connect along each course segment by increasing beta -----------------
        const auto& coursePairs = courseLineSegPairs[f];
        for (const auto& cp : coursePairs) {
            const Vector3 p0 = cp.first->baryCoords;
            const Vector3 p1 = cp.second->baryCoords;
            if ((p1 - p0).norm() < 1e-18) continue;    // skip degenerate

            std::vector<KGVertex*> stripe;
            stripe.reserve(F.size());
            for (KGVertex* v : F) {
                double t;
                if (onSegmentParam(v->baryCoords, p0, p1, t)) stripe.push_back(v);
            }

            std::sort(stripe.begin(), stripe.end(),
                      [](KGVertex* a, KGVertex* b) { return a->beta_tag < b->beta_tag; });

            for (size_t i = 0; i + 1 < stripe.size(); ++i) {
                KGVertex* a = stripe[i];
                KGVertex* b = stripe[i + 1];
                a->row_out_vertex = b;
                b->row_in_vertex  = a;
            }
        }

        // ----------------- Wale: connect along each wale segment by increasing alpha -------------------
        const auto& walePairs = waleLineSegPairs[f];
        for (const auto& wp : walePairs) {
            const Vector3 q0 = wp.first->baryCoords;
            const Vector3 q1 = wp.second->baryCoords;
            if ((q1 - q0).norm() < 1e-18) continue;

            std::vector<KGVertex*> stripe;
            stripe.reserve(F.size());
            for (KGVertex* v : F) {
                double t;
                if (onSegmentParam(v->baryCoords, q0, q1, t)) stripe.push_back(v);
            }

            std::sort(stripe.begin(), stripe.end(),
                      [](KGVertex* a, KGVertex* b) { return a->alpha_tag < b->alpha_tag; });

            for (size_t i = 0; i + 1 < stripe.size(); ++i) {
                KGVertex* a = stripe[i];
                KGVertex* b = stripe[i + 1];
                a->col_out_vertex[0] = b;
                b->col_in_vertex[0]  = a;
            }
        }
    }

    for (Face f : gluedGeometry->mesh.faces()) {
        auto& F = adjustedFaceKGVertices[f];
        for (KGVertex* v : F){
            if (v->row_out_vertex != nullptr) ensure (v -> beta_tag < v -> row_out_vertex -> beta_tag && "row ordering assertion failed on adjusted vertices");
            if (v->col_out_vertex[0] != nullptr) ensure (v -> alpha_tag < v -> col_out_vertex[0] -> alpha_tag && "column ordering assertion failed on adjusted vertices");
        }
    }
}

void KG::adjustedIntrinsicMerge(){

    //before we do anything else 
    //ensure that all real vertices have connections 
    for (auto& up : adjustedVertices) {
        KGVertex* v = up.get();
        if (v->isAlphaVirtual || v->isBetaVirtual) continue;
        ensure(v->row_in_vertex != nullptr && "real vertex doesn't have row_in set");
        ensure(v->row_out_vertex != nullptr && "real vertex doesn't have row_out set");
        ensure(v->col_in_vertex[0] != nullptr && "real vertex doesn't have col_in[0] set");
        ensure(v->col_out_vertex[0] != nullptr && "real vertex doesn't have col_out[0] set");
    }

    //also ensure all the ordering is correct
    for (Face f : gluedGeometry->mesh.faces()) {
        auto& F = adjustedFaceKGVertices[f];
        for (KGVertex* v : F){
            if (v->row_out_vertex != nullptr) ensure (v -> beta_tag < v -> row_out_vertex -> beta_tag && "row ordering assertion failed on adjusted vertices");
            if (v->col_out_vertex[0] != nullptr) ensure (v -> alpha_tag < v -> col_out_vertex[0] -> alpha_tag && "column ordering assertion failed on adjusted vertices");
        }
    }

    std::cout << "Ordering assertion passed" << std::endl;

    // const double epsT = 1e-12; // tolerance for "same place" along the edge

    // // Param along an oriented halfedge (same convention you used earlier)
    // auto edgeParam = [](KGVertex* v) -> double {
    //     Halfedge he = v->halfedge.value();
    //     Face f = he.face();
    //     Halfedge h0 = f.halfedge();
    //     Halfedge h1 = h0.next();
    //     Halfedge h2 = h1.next();
    //     if (he == h0)      return v->baryCoords[1]; // (i->j): t = b_j
    //     else if (he == h1) return v->baryCoords[2]; // (j->k): t = b_k
    //     else               return v->baryCoords[0]; // (k->i): t = b_i
    // };

    // // Store these vertices in order of the "direction of the halfedge"
    // std::map<Halfedge, std::vector<KGVertex*>> halfedgeCourseVertices;
    // std::map<Halfedge, std::vector<KGVertex*>> halfedgeWaleVertices;
    // for (const auto& v : adjustedVertices){
    //     if (v->isAlphaVirtual) halfedgeCourseVertices[v->halfedge.value()].push_back(v.get());
    //     if (v->isBetaVirtual) halfedgeWaleVertices[v->halfedge.value()].push_back(v.get());
    // }

    // // Connect course vertices 
    // for (Edge e : (gluedGeometry->mesh).edges()) {
    //     if (e.isBoundary()) continue;//don't need to handle boundary vertices

    //     std::vector<KGVertex*> he1Vertices = halfedgeCourseVertices[e.halfedge()];
    //     std::vector<KGVertex*> he2Vertices = halfedgeCourseVertices[e.halfedge().twin()];
        
    //     // Matchings across regular edges
    //     if (courseSingularEdgesGlued[e] == 0){ 
        
    //         ensure(he1Vertices.size() == he2Vertices.size() && "non course singular edge has unequal number of virtual vertices on either side of the halfedge");
    //         std::vector<std::pair<int, int>> regularMatchings;
    //         for (int i = 0; i < he1Vertices.size(); i++) {
    //             regularMatchings.push_back({i, (he1Vertices.size() - i) - 1});
    //         }

    //         // Connect the virtual matchings first
    //         for (auto [i1, i2] : regularMatchings) {
    //             KGVertex* v1 = he1Vertices[i1];
    //             KGVertex* v2 = he2Vertices[i2];
    //             ensure(v1->isAlphaVirtual && "vertex on halfedge is not virtual");
    //             ensure(v2->isAlphaVirtual && "vertex on halfege is not virtual");
    //             if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
    //                 v2->row_out_vertex = v1;
    //                 v1->row_in_vertex = v2;
    //             }
    //             if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
    //                 v1->row_out_vertex = v2;
    //                 v2->row_in_vertex = v1;                   
    //             }
    //         }
    //     }
    //     else{
            
    //         if (he1Vertices.size() > he2Vertices.size())
    //             swap(he1Vertices, he2Vertices);

    //         // he1Vertices.size() < he2Vertices.size() guaranteed by your swap above
    //         ensure(he2Vertices.size() - he1Vertices.size() == 1 && "More than one stripe born/dying at singular edge");

    //         // Build params for both sides in the SAME orientation (he1's)
    //         std::vector<std::pair<double, KGVertex*>> L, R;
    //         L.reserve(he1Vertices.size());
    //         R.reserve(he2Vertices.size());

    //         for (KGVertex* v : he1Vertices) L.emplace_back(edgeParam(v), v);
    //         for (KGVertex* v : he2Vertices) R.emplace_back(1.0 - edgeParam(v), v); // flip twin into he1 frame

    //         std::sort(L.begin(), L.end(), [](auto& a, auto& b){ return a.first < b.first; });
    //         std::sort(R.begin(), R.end(), [](auto& a, auto& b){ return a.first < b.first; });

    //         // Two-pointer: for every v1 in he1, find exactly-one v2 in he2 at the same param
    //         size_t i = 0, j = 0;
    //         while (i < L.size() && j < R.size()) {
    //             double tL = L[i].first;
    //             double tR = R[j].first;

    //             if (std::abs(tL - tR) <= epsT) {
    //                 KGVertex* v1 = L[i].second;   // on he1
    //                 KGVertex* v2 = R[j].second;   // on he2 (twin)

    //                 ensure(v1->isAlphaVirtual && v2->isAlphaVirtual && "expected alpha-virtuals on course singular edge");

    //                 if (v1->row_in_vertex == nullptr && v2->row_in_vertex != nullptr){
    //                     v2->row_out_vertex = v1;
    //                     v1->row_in_vertex = v2;
    //                 }
    //                 if (v2->row_in_vertex == nullptr && v1->row_in_vertex != nullptr){
    //                     v1->row_out_vertex = v2;
    //                     v2->row_in_vertex = v1;                   
    //                 }
    //                 ++i; ++j;
    //             } else if (tL < tR) {
    //                 ++i; // advance left to catch up
    //             } else {
    //                 ++j; // advance right to catch up
    //             }
    //         }
    //     }
    // }

    // // Connect wale vertices 
    // for (Edge e : (gluedGeometry->mesh).edges()) {
    //     if (e.isBoundary()) continue;//don't need to handle boundary vertices

    //     std::vector<KGVertex*> he1Vertices = halfedgeWaleVertices[e.halfedge()];
    //     std::vector<KGVertex*> he2Vertices = halfedgeWaleVertices[e.halfedge().twin()];
        
    //     // Matchings across regular edges
    //     if (waleSingularEdgesGlued[e] == 0){ 
        
    //         ensure(he1Vertices.size() == he2Vertices.size() && "non course singular edge has unequal number of virtual vertices on either side of the halfedge");
    //         std::vector<std::pair<int, int>> regularMatchings;
    //         for (int i = 0; i < he1Vertices.size(); i++) {
    //             regularMatchings.push_back({i, (he1Vertices.size() - i) - 1});
    //         }

    //         // Connect the virtual matchings first
    //         for (auto [i1, i2] : regularMatchings) {
    //             KGVertex* v1 = he1Vertices[i1];
    //             KGVertex* v2 = he2Vertices[i2];
    //             ensure(v1->isBetaVirtual && "vertex on halfedge is not virtual");
    //             ensure(v2->isBetaVirtual && "vertex on halfege is not virtual");
    //             if (v1->col_in_vertex[0] == nullptr && v2->col_in_vertex[0] != nullptr){
    //                 v1->col_in_vertex[0] = v2;
    //                 v2->col_out_vertex[0] = v1;
    //             }
    //             if (v2->col_in_vertex[0] == nullptr && v1->col_in_vertex[0] != nullptr){
    //                 v1->col_out_vertex[0] = v2;
    //                 v2->col_in_vertex[0] = v1;
    //             }
    //         }
    //     }
    //     else{
            
    //         if (he1Vertices.size() > he2Vertices.size())
    //             swap(he1Vertices, he2Vertices);

    //         // he1Vertices.size() < he2Vertices.size() guaranteed by your swap above
    //         ensure(he2Vertices.size() - he1Vertices.size() == 1 && "More than one stripe born/dying at singular edge");

    //         // Build params for both sides in the SAME orientation (he1's)
    //         std::vector<std::pair<double, KGVertex*>> L, R;
    //         L.reserve(he1Vertices.size());
    //         R.reserve(he2Vertices.size());

    //         for (KGVertex* v : he1Vertices) L.emplace_back(edgeParam(v), v);
    //         for (KGVertex* v : he2Vertices) R.emplace_back(1.0 - edgeParam(v), v); // flip twin into he1 frame

    //         std::sort(L.begin(), L.end(), [](auto& a, auto& b){ return a.first < b.first; });
    //         std::sort(R.begin(), R.end(), [](auto& a, auto& b){ return a.first < b.first; });

    //         // Two-pointer: for every v1 in he1, find exactly-one v2 in he2 at the same param
    //         size_t i = 0, j = 0;
    //         while (i < L.size() && j < R.size()) {
    //             double tL = L[i].first;
    //             double tR = R[j].first;

    //             if (std::abs(tL - tR) <= epsT) {
    //                 KGVertex* v1 = L[i].second;   // on he1
    //                 KGVertex* v2 = R[j].second;   // on he2 (twin)

    //                 ensure(v1->isBetaVirtual && v2->isBetaVirtual && "expected beta-virtuals on course singular edge");

    //                 if (v1->col_in_vertex[0] == nullptr && v2->col_in_vertex[0] != nullptr){
    //                     v1->col_in_vertex[0] = v2;
    //                     v2->col_out_vertex[0] = v1;
    //                 }
    //                 if (v2->col_in_vertex[0] == nullptr && v1->col_in_vertex[0] != nullptr){
    //                     v1->col_out_vertex[0] = v2;
    //                     v2->col_in_vertex[0] = v1;
    //                 }
    //                 ++i; ++j;
    //             } else if (tL < tR) {
    //                 ++i; // advance left to catch up
    //             } else {
    //                 ++j; // advance right to catch up
    //             }
    //         }
    //     }
    // }

    // // Now connect real vertices to one another
    // for (auto& up : adjustedVertices) {
    //     KGVertex* v0 = up.get();
    //     if (v0->isAlphaVirtual || v0->isBetaVirtual) continue; // only real vertices

    //     // --------Course--------
    //     {
    //         KGVertex* v = v0->row_out_vertex;     // start from immediate neighbor
    //         bool isGluedPath = false;

    //         while (v && v->isAlphaVirtual) {
    //             if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true;
    //             v = v->row_out_vertex;            // step
    //         }

    //         if (!v || v->isAlphaVirtual) {
    //             v0->row_out_vertex = nullptr;     // no real neighbor reachable
    //         } else {
    //             v0->row_out_vertex = v;           // connect reciprocally
    //             v->row_in_vertex   = v0;
    //             if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
    //         }
    //     }

    //     // --------Wale--------
    //     {
    //         KGVertex* v = v0->col_out_vertex[0];
    //         bool isGluedPath = false;

    //         while (v && v->isBetaVirtual) {
    //             if (v->halfedge && isGlued[v->halfedge->edge()]) isGluedPath = true;
    //             v = v->col_out_vertex[0];
    //         }

    //         if (!v || v->isBetaVirtual) {
    //             v0->col_out_vertex[0] = nullptr;
    //         } else {
    //             v0->col_out_vertex[0] = v;
    //             v->col_in_vertex[0]   = v0;
    //             if (isGluedPath) stitchedVertices.emplace_back(v0->id, v->id);
    //         }
    //     }
    // }

}