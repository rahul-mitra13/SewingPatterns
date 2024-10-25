#include "KnitGraph.h"

KnitGraph::KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, 
                        double period, EdgeData<double>& courseOneForm, EdgeData<double>& waleOneForm){

    this->globalGeometry = &globalGeometry; 
    this->gluedGeometry = &gluedGeometry;
    this->psMesh = &psMesh; 
    this->period = period; 
    this->courseOneForm = courseOneForm;
    this->waleOneForm = waleOneForm;
}

void KnitGraph::buildGraph(){

    SurfaceMesh& gluedMesh = this->gluedGeometry->mesh;

}