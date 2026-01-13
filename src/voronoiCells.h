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
  size_t iterations = 300;                 // number of iterations to run for
  size_t descIter = 10;                   // number of steps of gradient descent for weights update 
  double stepSize = 0.8;                  // step size for steps towards cell centers
  bool useDelaunay = true;                // solve on an intrinsic Delaunay triangulation of the input
  bool computeDistributions = false;      // return the indicator functions for each cell (`result.siteDistributions`)
  double tCoef = 1.0;                     // diffusion time for vector heat method
  size_t nSubIterations = 1;              // number of iterations to use when computing Karcher means
  double eps = 1e-6;                      // stopping criterion on the sum of Karcher mean updates
  int seed = 42;                          // random seed for site intialization
  bool useLineSearch = false;             // whether to use line search for weights update
  VertexData<double> measure;             // the measure we're trying to quantize
  double shortTime;                       // short time parameter for heat diffusion
  double epsWeights = 1e-8;               // stopping criterion on the weight optimization gradient norm
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
void alignPointsOnIsolineFast(SurfaceMesh& mesh, VertexPositionGeometry& geom, alignOptions& options, polyscope::SurfaceMesh &psMesh);


//perform line search by computing objective 
VoronoiResult computeSitesWithFunction(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                      VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh);

//computing equal weights using LBFGS
std::vector<double> computePhiWeights(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, std::vector<double>& initWeights, std::vector<VectorHeatMethodSolver>& vSolvers, VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh);


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

      for (int i = 0; i < options.nSites; i++){
        this->logMapPerSite.emplace_back(VertexData<Vector2>(this->mesh, Vector2{0., 0.}));
      }
      
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

      std::vector<VertexData<double>> rhs;
      for (int i = 0; i < options.nSites; i++) {
        this->heatKernel.emplace_back(VertexData<double>(this->mesh, 0.0));
        rhs.push_back(computeRHS(std::vector<SurfacePoint>{sites[i]}));
      }


      #pragma omp parallel for 
      for (size_t iSite = 0; iSite < options.nSites; iSite++) {
        int tid = omp_get_thread_num(); // thread ID
        //Just do a scalar diffuse
        this->heatKernel[iSite] = vSolvers[tid].scalarDiffuse(rhs[iSite]);
      }
    }

    double operator()(const Eigen::VectorXd& phiWeights, Eigen::VectorXd& grad){ 

      geom.requireVertexDualAreas();
      
      //we use this for numerical stability
      double maxWeight = *std::max_element(phiWeights.begin(), phiWeights.end());
      
      //find the desired mass (per site)
      double desiredMass = 0;
      for (Vertex v : mesh.vertices()){
        desiredMass +=  geom.vertexDualAreas[v] * options.measure[v];
      }
      desiredMass /= options.nSites;
  
      // NEW OBJECTIVE (ENTROPIC OT)
      double obj = 0;
      double eps = 4*options.shortTime;
      grad.setZero();
      for (Vertex v : mesh.vertices()) {
        std::vector<double> argExp(options.nSites);
        // double innerProduct;
        for (int iSite = 0; iSite < options.nSites; iSite++) {
          argExp[iSite] = phiWeights[iSite] / eps;
        }

        // adjust the argExp for numerical stability
        double maxArgExp = *std::max_element(argExp.begin(), argExp.end());
        for (int iSite = 0; iSite < options.nSites; iSite++) {
          argExp[iSite] -= maxArgExp;
        }

        // double maxArgExp = 0; // todo later
        double sumExp = 0;
        for (int iSite = 0; iSite < options.nSites; iSite++) {
          sumExp += this->heatKernel[iSite][v] * exp(argExp[iSite]);
        }
        for (int iSite = 0; iSite < options.nSites; iSite++) {
          grad(iSite) += - (this->heatKernel[iSite][v] * exp(argExp[iSite]) / sumExp) * options.measure[v] * geom.vertexDualAreas[v];
        }
        obj += -eps * (maxArgExp + log(sumExp * desiredMass)) * options.measure[v] * geom.vertexDualAreas[v]; // adding maxArgExp here compensates for the subtraction before          
      }

      // Add second term
      for (int iSite = 0; iSite < options.nSites; iSite++) {
        obj += phiWeights[iSite] * desiredMass;
        grad(iSite) += desiredMass;
      }
      // Negate everything to make it a mimization problem
      grad = -grad;
      return -obj;
  }

  double checkGrad(Eigen::VectorXd& phi, double h) {
    int n = phi.size();
    Eigen::VectorXd grad(n), dummy(n), grad_FD(n);
    double f = (*this)(phi, grad); // analytical gradient
    double error = 0;
    for (int i = 0; i < n; i++) {
      phi(i) += h;
      double fr = (*this)(phi, dummy); // f(x+h)
      phi(i) -= 2*h;
      double fl = (*this)(phi, dummy); // f(x-h)
      phi(i) += h; // reset to initial value
      grad_FD(i) = (fr - fl) / (2*h);
    }
    return (grad - grad_FD).norm() / grad.norm();
  }
};

} // namespace surface
} // namespace geometrycentral
