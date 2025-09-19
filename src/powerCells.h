#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct powerCellOptions {
  std::unique_ptr<EdgeLengthGeometry> gluedMesh; //intrinsic mesh 
  std::unique_ptr<ManifoldSurfaceMesh> globalMesh; //extrinsic mesh 
  std::unique_ptr<VertexPositionGeometry> globalGeometry; //extrinsic geometry 
  VertexData<double> courseCurlMeasure; //curl measure in the course direction 
  VertexData<double> waleCurlMeasure; //curl measure in the wale direction 
  std::map<int, int> vertexMap;//index from vertex in original mesh to vertex in glued mesh (probably don't need this and should just do everything in the glued setting)
  std::map<int, int> edgeMap;//index from an edge in the original mesh to to an edge in the glued mesh (probably don't need this and shouldn just do everything in the glued setting)
};

struct powerCellResults {
  EdgeData<double> courseSingularities;//accepted couse singularities on the glued mesh
  CornerData<double> courseOneForm;//final course stripes 
  EdgeData<double> waleSingularities;//accepted wale singularities on the glued mesh 
  EdgeData<double> waleOneForm;//final wale stripes
};

void computeSingularities(powerCellOptions& options);
