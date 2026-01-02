#pragma once
#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct powerCellOptions {
  EdgeLengthGeometry* gluedGeometry; //intrinsic mesh 
  VertexPositionGeometry* globalGeometry; //extrinsic geometry 
  std::map<int, int> vertexMap;//index from vertex in original mesh to vertex in glued mesh (probably don't need this and should just do everything in the glued setting)
  std::map<int, int> edgeMap;//index from an edge in the original mesh to to an edge in the glued mesh (probably don't need this and shouldn just do everything in the glued setting)
  std::vector<std::vector<double>> saddleLoops;//saddle loops for this geometry
  std::vector<Vertex> heatSourceVerts; //source vertices for the heat method
  std::map<int, std::vector<Halfedge>> gluedOneRingMap;//one ring map in the glued setting
  FaceData<Vector3> normalizedTFGrad; //normalized time function gradient (in the global setting)
  bool maskSaddle = true;//masking the saddle vertices
  polyscope::SurfaceMesh* psMesh;//the polyscope mesh
  double period;//period of the stripe patterns
  VertexData<double> timeFunction;//time function over the mesh 

};

struct ComponentMesh {
    std::unique_ptr<SurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;

    std::unordered_map<Vertex, Vertex> oldToNewV;
    std::vector<Vertex> newToOldV;

    std::unordered_map<Edge, Edge> oldToNewE;
    std::vector<Edge> newToOldE;
};

// Note: surface points are on edges, not faces!
// This is a result of how projectOnIsoline() works.
std::vector<std::pair<SurfacePoint, SurfacePoint>>  computeCourseSingularities(powerCellOptions& options);

//compute the course curl measure in the glued setting 
VertexData<double> computeCourseCurlMeasure(powerCellOptions& options);

//compute the singularity positions per bucket
std::vector<std::pair<SurfacePoint, SurfacePoint>> computeBucketSingularities(powerCellOptions& options, VertexData<double>& curlMeasure, VertexData<double>& allDistance, 
                                                                              bool cylindricalParameterization);

//Create a submesh without preserving boundaries i.e., nBoundaryLoops in each submesh = 0
ComponentMesh extractComponent(SurfaceMesh& mesh, VertexPositionGeometry& geom, const std::vector<Face>& faces);

