#pragma once

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include "polyscope/surface_mesh.h"

namespace geometrycentral {
namespace surface {

struct VoronoiResult {
  std::vector<SurfacePoint> siteLocations;               // sites at the centers of the Voronoi cells
  std::vector<VertexData<double>> siteDistributions;    // soft indicator functions for each Voronoi cell
  bool hasDistributions = false;                       // is siteDistributions populated?
  std::vector<std::vector<SurfacePoint>> steps;       // DEBUG: steps taken from each starting point
  std::vector<std::vector<VertexData<double>>> stepSiteDistribution; // DEBUG: visualize the size distrubutions during the steps
  std::vector<SurfacePoint> initialSites;             //DEBUG: location of initial sites
};

struct VoronoiOptions {
  size_t nSites = 1;                     // number of sites to place
  std::vector<SurfacePoint> initialSites; // desired locations for sites. If blank, locations are chosen randomly
  size_t iterations = 500;                 // number of iterations to run for
  double stepSize = 1.0;                    // step size for steps towards cell centers
  bool useDelaunay = true;                // solve on an intrinsic Delaunay triangulation of the input
  bool computeDistributions = false;      // return the indicator functions for each cell (`result.siteDistributions`)
  double tCoef = 1.0;                       // diffusion time for vector heat method
  size_t nSubIterations = 1;              // number of iterations to use when computing Karcher means
  double eps = 1e-6;                      // stopping criterion on the sum of Karcher mean updates
  int seed = 42;                          // random seed for site intialization
  // double epsWeights = 1e-5;               // stopping criterion on the weight optimization gradient norm
};


struct alignOptions {

  VertexData<double> timeFunction;  //time function over the mesh 
  std::vector<std::pair<SurfacePoint, SurfacePoint>> pairedSites;   //paired sites obtained from the previous optimization stored in (pos site, neg site) order 
  bool usingPosCurl = false; 
  double lambda = 1;                      //weighting term
  size_t iterations = 5;                  // number of iterations to run for
  double tCoef = 1;                       // diffusion time for vector heat method
  size_t nSubIterations = 1;              // number of iterations to use when computing Karcher means
  double stepSize = 1;                    // step size for steps towards cell centers
  bool computeDistributions = false;      // return the indicator functions for each cell (`result.siteDistributions`)
};



extern const VoronoiOptions defaultVoronoiOptions;

VoronoiResult computeGeodesicCentroidalVoronoiTessellationWithWeights(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh);

// Given a set of singularity pairs in `options`, aligns them (in place) to the midpoint of their time values.
// Based on geodesic tracing so not exact, but fast.
void alignPointsOnIsolineFast(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, alignOptions& options, polyscope::SurfaceMesh &psMesh);


//TO-DO: jointly optmize the positions of both the positive sites and the negtive site 
VoronoiResult alignPointsOnIsoline(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                  alignOptions options, VertexData<double>& posMeasure, VertexData<double>& negMeasure, polyscope::SurfaceMesh &psMesh);


} // namespace surface
} // namespace geometrycentral
