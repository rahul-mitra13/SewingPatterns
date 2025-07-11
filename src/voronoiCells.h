#pragma once

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include "polyscope/surface_mesh.h"
#include <LBFGS.h>

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
  size_t nSites = 1;                      // number of sites to place
  std::vector<SurfacePoint> initialSites; // desired locations for sites. If blank, locations are chosen randomly
  size_t iterations = 100;                // number of iterations to run for
  size_t descIter = 10;                   // number of steps of gradient descent for weights update 
  double stepSize = 1.0;                  // step size for steps towards cell centers
  bool useDelaunay = true;                // solve on an intrinsic Delaunay triangulation of the input
  bool computeDistributions = false;      // return the indicator functions for each cell (`result.siteDistributions`)
  double tCoef = 1.0;                     // diffusion time for vector heat method
  size_t nSubIterations = 1;              // number of iterations to use when computing Karcher means
  double eps = 1e-6;                      // stopping criterion on the sum of Karcher mean updates
  int seed = 42;                          // random seed for site intialization
  bool useLineSearch = false;             // whether to use line search for weights update
  VertexData<double> measure;             // the measure we're trying to quantize
  double shortTime;                       // short time parameter for heat diffusion
  // double epsWeights = 1e-5;            // stopping criterion on the weight optimization gradient norm
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


//perform line search by computing objective 
VoronoiResult computeSitesWithFunction(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                      VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh);


//defining a class to do LBFGS with F_OT as defined by the seminal work, "Stochastic Wassertein Barycenters"
class F_OT{
  private: 
    IntrinsicGeometryInterface& geom; 
    SurfaceMesh& mesh;
    VoronoiOptions& options;
    std::vector<VertexData<Vector2>> logMapPerSite;
    std::vector<VertexData<double>> heatKernel;
    std::vector<SurfacePoint> sites;

  public: 
    F_OT(IntrinsicGeometryInterface& myGeom, SurfaceMesh& myMesh, VoronoiOptions& myOptions, std::vector<SurfacePoint>& mySites):
      geom(myGeom), mesh(myMesh), options(myOptions), sites(mySites){}
    
    void requireLogMap(std::vector<VectorHeatMethodSolver>& vSolvers){
      
      #pragma omp parallel for 
      for (size_t iSite = 0; iSite < options.nSites; iSite++){
        int tid = omp_get_thread_num(); // thread ID
        SurfacePoint site = sites[iSite];
        this->logMapPerSite[iSite] = vSolvers[tid].computeLogMap(site);
      }
    }

    void requireHeatKernel(std::vector<VectorHeatMethodSolver>& vSolvers){

      auto computeRHS = [&](const std::vector<SurfacePoint>& points) -> VertexData<double> {
        VertexData<double> rhs(mesh, 0);
        for (const SurfacePoint& point : points) {
          SurfacePoint facePoint = point.inSomeFace();
          Halfedge he = facePoint.face.halfedge();

          rhs[he.vertex()] += facePoint.faceCoords.x;
          rhs[he.next().vertex()] += facePoint.faceCoords.y;
          rhs[he.next().next().vertex()] += facePoint.faceCoords.z;
        }
        return rhs;
      };

      #pragma omp parallel for 
      for (size_t iSite = 0; iSite < options.nSites; iSite++) {
        int tid = omp_get_thread_num(); // thread ID
        SurfacePoint site = sites[iSite];
        //Just do a scalar diffues
        this->heatKernel[tid] = vSolvers[tid].scalarDiffuse(computeRHS(std::vector<SurfacePoint>{site}));
      }
    }

    double operator()(const Eigen::VectorXd& phiWeights, Eigen::VectorXd& grad){ 
      
      double maxWeight = *std::max_element(phiWeights.begin(), phiWeights.end());
      std::vector<VertexData<double>> rho(options.nSites);

      //find the desired mass 
      double desiredMass = 0;
      for (Vertex v : mesh.vertices()){
        desiredMass +=  geom.vertexDualAreas[v] * options.measure[v];
      }
      desiredMass /= options.nSites;


      geom.requireVertexDualAreas();
      VertexData<double> normD(mesh, 0); // denominator for the rho's
      for (size_t iSite = 0; iSite < options.nSites; iSite++) {
        normD += exp((phiWeights[iSite] - maxWeight)/(4 * options.shortTime)) * this->heatKernel[iSite];
      }

      VertexData<double> integrand;
      #pragma omp parallel for 
      for (size_t iSite = 0; iSite < options.nSites; iSite++){

        VertexData<double> d2minusPhi(mesh, 0); 
        for (Vertex v :  mesh.vertices()){
          d2minusPhi[v] = (logMapPerSite[iSite][v].norm2() - phiWeights[iSite]);
        }
        
        VertexData<double> rho = exp((phiWeights[iSite] - maxWeight)/(4 * options.shortTime)) * this->heatKernel[iSite];

        integrand = d2minusPhi * rho * options.measure * geom.vertexDualAreas;

        double sumTerm = 0;
        for (Vertex v : mesh.vertices()){
          sumTerm += rho[v] * options.measure[v] * geom.vertexDualAreas[v];
        }

        //add the gradient component 
        grad(iSite) = desiredMass - sumTerm;
      }

      //compute the first term 
      double weightSum = 0; 
      for (int i = 0; i < phiWeights.size(); i++) weightSum += phiWeights(i);
      double firstTerm = weightSum * desiredMass;

      //compute the second term
      double secondTerm = 0;
      for (Vertex v : mesh.vertices()) secondTerm += integrand[v];

      double obj = firstTerm + secondTerm;

      return obj;
  }
};

} // namespace surface
} // namespace geometrycentral
