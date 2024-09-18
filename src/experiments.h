//file includes
#include "knitting_utils.h"
#include "helpers.h"
#include "stripe_patterns_helpers.h"
#include <tuple>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//-----------------------------strategy1-------------------------------------------//
//compute the iterative 1-form and vertex singularity positions trying the strategy 1 
//strategy 1
//1. Sample level sets at intervals of P, say C_i 
//2. For each C_i, consider all the level sets it passes through, say T_i 
//3. Calc average ||\frac{\nabla h}{||\nabla h}||^2 - \grad \sigma_j} (probably area-weighted) over T_i
//4. For C_i with largest average ||\frac{\nabla h}{||\nabla h||, search max/min singularities and insert 
//5. Stop when ||\delta \sigma - \frac{\nabla h}{||\nabla h}||^2 doesn't decrease any more 
std::tuple<HalfedgeData<double>, VertexData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    FaceData<Vector3>& courseFaceGradients);

//compute the per-face gradient of a 1-form
FaceData<Vector3> computeOneFormFaceGrad(VertexPostionGeometry& globalGeometry, HalfedgeData<double>& sigmaTilde);