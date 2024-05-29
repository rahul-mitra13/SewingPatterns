//eigen includes
#include <Eigen/Dense>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

//libigl includes
#include <igl/cotmatrix.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/grad.h>
#include <igl/average_onto_vertices.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//this is a utility file that will handle basic mesh processing tasks such as 
//1. solveLaplace()
//2. computeTimeFunctionGrad()

//
//solve Laplace equation over an input mesh 
//
//@param[in]    geometry        VertexPositionGeometry          input geometry
//@param[in]    zeroVertices    std::vector<Vertex>             vertices that take the value 0 as a boundary condition i.e., where knitting should start
//@paramm[in]   oneVertices     std::vector<Vertex>             vertices that take the value 1 as a boundary condition i.e., where knitting should stop 
//
//@return       timeFunction    VertexData<double>              #V by 1 scalar function that represents the knitting time function
VertexData<double> solveLaplace(VertexPositionGeometry& geometry, std::vector<Vertex>& zeroVertices, std::vector<Vertex>& oneVertices);


//compute the per face gradient of some function defined as a scalar over vertices
//
//@param[in]    geometry                VertexPositionGeometry          input geometry
//@param[in]    vertexScalarFunction    VertexData<double>              scalar function defined as a scalar over vertices
//
//@return       faceGradients           FaceData<Vector3>               #F by 3 vector per face
FaceData<Vector3> computeTimeFunctionFaceGrad(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction);

//compute a vector (in ambient space) per vertex that is aligned with the gradient of a scalar field
//
//@param[in]    geometry                VertexPositionGeometry  input geometry
//@param[in]    vertexScalarFunction    VertexData<double>      scalar function defined as a scalar over vertices
//
//@return       vertexValuedField       VertexData<Vector3>     #V by 3 vector per face     
VertexData<Vector3> computeVertexValuedField(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction);

