#include "experiments.h"

std::tuple<HalfedgeData<double>, VertexData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    FaceData<Vector3>& courseFaceGradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

}

//compute the per-face gradient of a 1-form
FaceData<Vector3> computeOneFormFaceGrad(VertexPostionGeometry& globalGeometry, HalfedgeData<double>& sigmaTilde){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
}