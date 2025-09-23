#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct powerCellOptions {
  EdgeLengthGeometry* gluedMesh; //intrinsic mesh 
  ManifoldSurfaceMesh* globalMesh; //extrinsic mesh 
  VertexPositionGeometry* globalGeometry; //extrinsic geometry 
  std::map<int, int> vertexMap;//index from vertex in original mesh to vertex in glued mesh (probably don't need this and should just do everything in the glued setting)
  std::map<int, int> edgeMap;//index from an edge in the original mesh to to an edge in the glued mesh (probably don't need this and shouldn just do everything in the glued setting)
  bool maskSaddleLoops;//should the saddle loops be masked
  bool maskBoundaries;//should boundaries be masked
};

struct powerCellResults {
  EdgeData<double> courseSingularities;//accepted couse singularities on the glued mesh
  CornerData<double> courseOneForm;//final course stripes 
  EdgeData<double> waleSingularities;//accepted wale singularities on the glued mesh 
  EdgeData<double> waleOneForm;//final wale stripes
};

void computeSingularities(powerCellOptions& options);
