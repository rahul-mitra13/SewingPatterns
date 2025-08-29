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

        //bring everything into the glued setting first
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
}

//Makes virtual vertices in the course direction
void KG::makeCourseVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeVirtualVerticesOnBorder(f, true);
    }

    std::vector<Vector3> virtualCourseSingularVertices;

    for (auto &v : allVertices){
        if (!v->isAlphaVirtual) continue;
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
        virtualCourseSingularVertices.emplace_back(v->position);
    }

    polyscope::registerPointCloud("Course Singular Vertices", virtualCourseSingularVertices);
}

//Makes virtual vertices in the wale direction
void KG::makeWaleVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeVirtualVerticesOnBorder(f, false);
    }

    std::vector<Vector3> virtualWaleSingularVertices;

    for (auto &v : allVertices){
        if (!v->isBetaVirtual) continue;
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
        virtualWaleSingularVertices.emplace_back(v->position);
    }

    polyscope::registerPointCloud("Wale Singular Vertices", virtualWaleSingularVertices);
}

//make virtual vertices on the border of all faces
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
        for (auto &v : allVertices) {
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
    double RHS_alpha, RHS_beta, detA, a1, b1, c1, a2, b2, c2, j, k, bi, bj, bk;
    
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
}

void KG::makeFaceConnections(){

    constexpr double eps = 1e-8;

    auto approx_contains = [eps](const std::vector<double>& xs, double x) {
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
                currVertex->row_out = nextVertex->id;
                nextVertex->row_in  = currVertex->id;
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
                currVertex->col_out[0] = nextVertex->id;
                nextVertex->col_in[0]  = currVertex->id;
                //update the pointers 
                currVertex->col_out_vertex[0] = nextVertex;
                nextVertex->col_in_vertex[0] = currVertex;
                
            }
        }
    }

}
