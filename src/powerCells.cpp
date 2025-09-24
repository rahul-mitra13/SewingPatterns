#include "powerCells.h"
void computeSingularities(powerCellOptions& options){

    const EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    const SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    const VertexPositionGeometry* globalGeometry = options.globalGeometry;
    const SurfaceMesh& globalMesh = globalGeometry->mesh;

    //first compute the measure in the course direction 
    VertexData<double> courseMeasure; 
}

//compute the course curl measure in the glued setting 
VertexData<double> computeCourseCurlMeasure(powerCellOptions& options){

    SurfaceMesh& globalMesh = options.globalGeometry->mesh;
    SurfaceMesh& gluedMesh = options.gluedGeometry->mesh;
    VertexData<double> curlGlobal(globalMesh);
    VertexData<double> curlGlued(gluedMesh);
    options.globalGeometry->requireFaceAreas();
    auto vertexMap = options.vertexMap;
    auto field = options.normalizedTFGrad;
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        double area = 0.0; // area of the 1-ring of faces
        for (Halfedge he : options.gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = options.globalGeometry->vertexPositions[hjk.tipVertex()] - options.globalGeometry->vertexPositions[hjk.tailVertex()];
            field[he.face()] = field[he.face()].normalize(); //always normalize the field
            area += options.globalGeometry->faceArea(he.face());
            sum += dot(hjkVec, field[he.face()]);
        }

        if (vi.isBoundary()) curlGlobal[vi] = 0.;//set curl to 0 on boundary vertices
        else curlGlobal[vi] = sum / area;
    }

    //convert curl to glued setting 
    for (Vertex v : globalMesh.vertices()){
        curlGlued[vertexMap[v.getIndex()]] = curlGlobal[v];
    }

    return curlGlued;

}