#include "KG.h"


KG::KG(VertexPositionGeometry& globalGeometry,
                     EdgeLengthGeometry& gluedGeometry,
                     polyscope::SurfaceMesh& psMesh,
                     double coursePeriod, double walePeriod,
                     CornerData<double>& courseOneForm,
                     EdgeData<double>& courseSingularEdges,
                     CornerData<double>& waleOneForm,
                     EdgeData<double>& waleSingularEdges,
                     std::map<int, int>& globalToGluedEdgeMap)
    : globalGeometry(&globalGeometry)
    , gluedGeometry(&gluedGeometry)
    , psMesh(&psMesh)
    , coursePeriod(coursePeriod)
    , walePeriod(walePeriod)
    , courseOneForm(courseOneForm)
    , courseSingularEdges(courseSingularEdges)
    , waleOneForm(waleOneForm)
    , waleSingularEdges(waleSingularEdges)
    , globalToGluedEdgeMap(&globalToGluedEdgeMap){

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
}

//Makes virtual vertices in the course direction
void KG::makeCourseVirtualVertices(){

    for (Face f : gluedGeometry->mesh.faces()){
        makeSingularFaceVirtualVertices(allVertices, f, coursePeriod, true);
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
        makeSingularFaceVirtualVertices(allVertices, f, walePeriod, false);
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

//make virtual vertices on the border of a singular face
void KG::makeSingularFaceVirtualVertices(std::vector<std::unique_ptr<KGVertex>>& allVertices, Face& f, double period, 
                                     bool isCourseDirection){
        
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
    
    if (isCourseDirection){//course direction
        //shift by small epsilon to allow for floating point tolerance
        for (j = alpha_start - eps; j < alpha_end + eps; j += coursePeriod){//fix alpha
            //ij edge
            bi = (j - alphaJ) / (alphaI - alphaJ);
            bj = 1.0 - bi;
            bk = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();						
   			    v->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
                v->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }

            //jk edge
            bj = (j - alphaK) / (alphaJ - alphaK);
            bk = 1.0 - bj;
            bi = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();			
   			    v->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
                v->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }

            //ki edge
            bi = (j - alphaK) / (alphaI - alphaK);
            bk = 1.0 - bi;
            bj = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();							
   			    v->baryCoords = Vector3{bi, bj, bk};
                k = bi * betaI + bj * betaJ + bk * betaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
                v->isAlphaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }
        }
    }else{//wale direction

        //shift up by some small epsilon so that we don't break when Z_start = Z_end
   	    for (k = beta_start - eps; k < beta_end + eps; k += walePeriod){        
            //ij edge
            bi = (k - betaJ) / (betaI - betaJ);
            bj = 1.0 - bi;
            bk = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();						
   			    v->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge();
                v->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }

            //jk edge
            bj = (k - betaK) / (betaJ - betaK);
            bk = 1.0 - bj;
            bi = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();						
   			    v->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next();
                v->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }

            //ki edge
            bi = (k - betaK) / (betaI - betaK);
            bk = 1.0 - bi;
            bj = 0.0;
            if ((bi >= 0.0 - eps && bi <= 1.0 + eps) && (bj >= 0.0 - eps && bj <= 1.0 + eps) && (bk >= 0.0 - eps && bk <= 1.0 + eps)){//only store points in the triangle
                auto v = std::make_unique<KGVertex>();						
   			    v->baryCoords = Vector3{bi, bj, bk};
                j = bi * alphaI + bj * alphaJ + bk * alphaK;
                v->id = vertexID++;
                v->alpha_tag = j;
                v->beta_tag = k;
                //edges and halfedges are in the glued mesh setting 
                v->halfedge = gluedGeometry->mesh.face(f.getIndex()).halfedge().next().next();
                v->isBetaVirtual = true;
                allVertices.emplace_back(std::move(v));
            }
        }
    }
}