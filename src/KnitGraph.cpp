#include "KnitGraph.h"

KnitGraph::KnitGraph(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, 
                        double period, EdgeData<double>& courseOneForm, EdgeData<double>& waleOneForm){

    this->geometry = &geometry; 
    this->psMesh = &psMesh; 
    this->period = period; 
    this->courseOneForm = courseOneForm;
    this->waleOneForm = waleOneForm;
}

void KnitGraph::buildGraph(){

    SurfaceMesh& mesh = this->geometry->mesh;

}