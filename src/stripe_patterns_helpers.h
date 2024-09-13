#pragma once

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/surface/stripe_patterns.h"
#include "geometrycentral/surface/edge_length_geometry.h"

#include "polyscope/surface_mesh.h"

#include <vector>
#include <queue>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Implementation of "Stripe Patterns on Surfaces" [Knoppel et al. 2015]
struct Isoline {
  std::vector<std::pair<Halfedge, double>> barycenters;
  bool open;
};

struct PolyLinePoint{
  Vector3 position;//the position of the point in space 
  Face f;//the face this point belongs to
  double isoval;//the isoval (modulo P of this point)
};

//do the integration in the glued mesh setting itself
//this is helpful because we do not need to pass around edge maps and vertex maps
//sigma is passed over the EDGES of the mesh
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& sigma, float period);

//carry out the integration in the glued mesh setting
//sigma is defined over the HALFEDGES of the glued mesh
//sigma is passed over the HALFEDGES of the mesh
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde, float period);

//computes the stripes texture coordinates from the passed 1-form
std::tuple<CornerData<double>, FaceData<int>> computeStripeValuesFromOneForm(IntrinsicGeometryInterface& geometry, EdgeData<double>& sigma, float period, SurfaceMesh& gluedMesh, 
                                                                            std::vector<std::pair<int, int>>& vertexMappingsPairs, std::vector<std::pair<int, int>>& edgeMappingsPairs, 
                                                                            std::map<int, std::vector<Halfedge>>& gluedOneRingMap);

                                                          

//returns a vector per face that stores all the isovalues (modulo period) passing through this face
FaceData<std::vector<double>> getFaceIsoValues(IntrinsicGeometryInterface& geometry,
                                              const CornerData<double>& stripeValues,
                                              const FaceData<int>& stripesIndices, double period);


//returns a vector per halfedge that stores all the isovalues (modulo period) passing through this halfedge
HalfedgeData<std::vector<double>> getHalfEdgeIsoValues(IntrinsicGeometryInterface& geometry,
                                              const CornerData<double>& stripeValues,
                                              const FaceData<int>& stripesIndices, double period);

//extract isolines when multiple isolines may pass through a face (using face-based matching)
std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>> generateIsoLines(EmbeddedGeometryInterface& geometry,
                                                      const CornerData<double>& stripeValues,
                                                      const FaceData<int>& stripesIndices, double period);
