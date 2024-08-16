//eigen includes
#include <Eigen/Dense>

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/surface/stripe_patterns.h"

//libigl includes
//I guess we'll just put all the libigl includes here for now 
#include <igl/cotmatrix.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/grad.h>
#include <igl/average_onto_vertices.h>
#include <igl/hessian_energy.h>
#include <igl/curved_hessian_energy.h>
#include <igl/isolines.h>

//helper files
#include "helpers.h"

//includes to solve the optimization problem
#include "gurobi_c++.h"
#include "model.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


//compute time function using a vector of pairs of vertex mappings instead of the map because we miss stitches then 
//
//@param[in]    geometry        VertexPositionGeometry                  input geometry
//@param[in]    vertexMappings  std::vector<std::pair<int,int>>         map that stores global vertex mappings 
//
//@return       timeFunction    VertexData<double>                      the global time function computed over the entire mesh
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::vector<std::pair<int,int>>& vertexMappingsPairs,  globalBoundaryConditions& boundaryConditions, std::map<int, int>& indexMap);

//compute time function directly in the glued mesh setting 
VertexData<double> computeTimeFunction(EdgeLengthGeometry& gluedGeometry, globalBoundaryConditions& bdyConditions);


//compute the per face gradient of some function defined as a scalar over vertices
//
//@param[in]    geometry                VertexPositionGeometry          input geometry
//@param[in]    vertexScalarFunction    VertexData<double>              scalar function defined as a scalar over vertices
//
//@return       faceGradients           FaceData<Vector3>               #F by 3 vector per face
FaceData<Vector3> computeTimeFunctionFaceGrad(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction);

//compute the gradient of the a function defined as a scalar over vertices in the glued mesh setting 
//Don't think the formula I'm using in here is right
FaceData<Vector3> computeTimeFunctionFaceGrad(EdgeLengthGeometry& geometry, VertexData<double>& vertexScalarFunction);

//compute a vector (in ambient space) per vertex that is aligned with the gradient of a scalar field
//
//@param[in]    geometry                VertexPositionGeometry  input geometry
//@param[in]    vertexScalarFunction    VertexData<double>      scalar function defined as a scalar over vertices
//@param[in]    angle                   double                  how much to rotate the vector about the vertex normal axis in radians
//
//@return       vertexValuedField       VertexData<Vector3>     #V by 3 vector per face     
VertexData<Vector3> computeVertexValuedField(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction, double angle);


//compute a direction field per vertex from a vector field (in ambient space) per vertex
//
//@param[in]    geometry            VertexPositionGeometry geometry     input geometry
//@param[in]    vertexValuedField   VertexData<Vector3>                 #V by 3 vector per vertex
//
//@return       directionField      VertexData<Vector2>                 direction field at each vertex
VertexData<Vector2> vertexDirectionField(VertexPositionGeometry& geometry, VertexData<Vector3>& vertexValuedField);

//compute a 1-form that will be used to generate the stripes over the quad patches 
//
//@param[in]    geometry            VertexPositionGeometry              input geometry 
//@param[in]    model               Model                               The gurobi model we're interested in solving
//@param[in]    edgeMappingPairs    std::vector<std::pair<int, int>>    vector of pairs where each pair representes "stitched together" edge indices   
//
//@return       oneForm     EdgeData<double>        The course 1-form used to specify the stripes
EdgeData<double> computeOneForm(VertexPositionGeometry& geometry, Model& gbModel, polyscope::SurfaceMesh& globalPSMesh);

//compute a 1-form that will be used to generate the stripes over the pathces 
//overload from the function above
//it does the optimization in the glued mesh setting 
EdgeData<double> computeOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& vertexMap, VertexData<double>& timeFunction, polyscope::SurfaceMesh& psMesh);

//compute \omega i.e., the 1-form we're trying to match over each edge 
//
//@param[in]    geometry            VertexPositionGeometry  input geometry
//@param[in]    direction           int                     integer representing the direction for omega that we're attempting to generate (0 -> course, 1 -> wale)
//@param[in]    faceGradients       FaceData<Vector3>       the gradients of the time function we're using to solve for the matching 1-form
//
//@return       matchingOneForm     EdgeData<double>        The 1-form we're trying to match 
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3>& faceGradients);

//compute matching 1-form while taking into account "stitched together" edges
//same parameters and return type as the function above 
//additional parameter 
//@param[in]    edgeMappingsPairs       std::vector<std::pair<int, int>>    vector of pairs where each pair stores edges that are stitched together
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, int direction, FaceData<Vector3>& faceGradients, std::vector<std::pair<int, int>>& edgeMappingsPairs);

