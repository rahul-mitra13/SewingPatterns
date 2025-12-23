#include "voronoiCells.h"
#include "math.h"
#include "float.h"
#include "knitting_utils.h"
#include "helpers.h"

using namespace std;

namespace geometrycentral {
namespace surface {

namespace {
const bool VORONOI_PRINT = false;
}

// To print vectors easily
template<class T> std::ostream &operator<<(std::ostream &os, std::vector<T> v) {
  os << "["; if (v.size() > 0) os << v[0];
  for(int i = 1; i < v.size(); i++) os << ", " << v[i];
  os << "]";
  return os;
}

// The default trace options
const VoronoiOptions defaultVoronoiOptions{};

VoronoiResult computeGeodesicCentroidalVoronoiTessellationWithWeights(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh) {

  // Define one vector heat method solver per thread
  vector<VectorHeatMethodSolver> vSolvers;
  for (int i = 0; i < omp_get_max_threads(); i++)
    vSolvers.emplace_back(geom, options.tCoef);
  VectorHeatMethodSolver& vSolver = vSolvers[0]; // temporarily, for the serial code


  // This computes e^(φ_j/4t) δ(x_j) for every site j
  auto computeRHSWithWeights = [&](const std::vector<SurfacePoint>& points, std::vector<double>& weights, double shortTime) -> std::vector<VertexData<double>> {

    std::vector<VertexData<double>> rhs;
    double maxWeight = *std::max_element(weights.begin(), weights.end());

    for (size_t i = 0; i < points.size(); i++){
      VertexData<double> rhsi(mesh, 0);
      SurfacePoint facePoint = points[i].inSomeFace();
      Halfedge he = facePoint.face.halfedge();
      rhsi[he.vertex()] += exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.x;
      rhsi[he.next().vertex()] +=  exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.y;
      rhsi[he.next().next().vertex()] += exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.z;
      rhs.push_back(rhsi);
    }

    return rhs;
  };

  //This computes the objective of the proxy function we're trying to use for line search 
  auto proxyEnergy = [&](const std::vector<double>& cellMasses, int LloydIter) -> double {

    if (LloydIter == 0) return DBL_MAX;
    double minCellMass = *min_element(cellMasses.begin(), cellMasses.end());
    double maxCellMass = *max_element(cellMasses.begin(), cellMasses.end());
    return std::fabs(maxCellMass - minCellMass);
  };
  

  // Set points to start
  // TODO: this should be done outside this function, ideally
  //use a random seed 
  std::mt19937 rng(options.seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::uniform_int_distribution<int> pointTypeDist(0, 2); // 0: vertex, 1: edge, 2: face
  for (size_t i = 0; i < options.nSites; ++i) {
      int type = pointTypeDist(rng);
      type = 2;

      if (type == 0) {
        // Random vertex
        Vertex v = mesh.vertex(rng() % mesh.nVertices());
        options.initialSites.emplace_back(v);
      } else if (type == 1) {
        // Random edge interior
        Edge e = mesh.edge(rng() % mesh.nEdges());
        double t = uniform01(rng); // interpolation along the edge
        options.initialSites.emplace_back(e, t);
      } else if (type == 2) {
        // Random face interior
        Face f = mesh.face(rng() % mesh.nFaces());

        // Random barycentric coordinates
        double u = uniform01(rng);
        double v = uniform01(rng);
        if (u + v > 1.0) {
          u = 1.0 - u;
          v = 1.0 - v;
        }
        double w = 1.0 - u - v;
        Halfedge he = f.halfedge();
        SurfacePoint pt(f, Vector3{u, v, w}); // SurfacePoint using barycentric coords in a triangle
        options.initialSites.push_back(pt);
      }

      // H(type);
      // H(siteLocations.back().face);
  }
  std::vector<SurfacePoint> siteLocations = options.initialSites;

  //debug on the square 
  //f2592 and f6688 are next to each other 
  // siteLocations.push_back(SurfacePoint(mesh.face(2592), Vector3{0.25, 0.25, 0.5}));
  // siteLocations.push_back(SurfacePoint(mesh.face(6688), Vector3{1./3., 1./3., 1./3.}));

  size_t nSites = siteLocations.size();
  VoronoiResult result;
  result.steps.resize(nSites);
  result.stepSiteDistribution.resize(nSites);
  result.initialSites = siteLocations;

  // Handle case where there's no sites at all
  if (siteLocations.empty()) {
    result.siteLocations = siteLocations;
    return result;
  }

  geom.requireShapeLengthScale();
  double M = 0; // maximum measure
  for (Vertex v : mesh.vertices())
    M = max(M, measure[v]);
  double L = geom.shapeLengthScale; // characteristic length
  double lips = sqrt(options.nSites) * M * L / 2; // Lipschitz constant

  // For some reason, running a dummy scalarDiffuse / computeLogMap first
  // fixes thread issues lol. Many thanks to
  // https://github.com/nmwsharp/geometry-central/issues/108
  VertexData<double> dummyRHS(mesh, 42.0);
  SurfacePoint dummySite = mesh.vertex(0);
  for (int i = 0; i < omp_get_max_threads(); i++) {
    vSolvers[i].scalarDiffuse(dummyRHS);
    vSolvers[i].computeLogMap(dummySite);
  }  

  geom.requireVertexDualAreas();

  // double stepSize = 1e-6;
  std::vector<double> phiWeights(nSites, 0.), oldPhiWeights(nSites);
  std::vector<double> cellMasses(nSites); //masses on each cell
  double shortTime;
  geom.requireEdgeLengths();
  //compute the diffusion time 
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = options.tCoef * meanEdgeLength * meanEdgeLength;

  //compute mass per cell 
  double desiredMass = 0.; 
  for (Vertex v : mesh.vertices()){
    desiredMass += geom.vertexDualAreas[v] * measure[v];
  }
  desiredMass /= nSites;
  
  // == Iterations
  bool stop = false; // flag if stopping criterion on sumUpdateNorm is reached
  double sumUpdateNorm = 1; // just a starting point to determine a tolerance on weights grad norm
  for (size_t iIter = 0; iIter < options.iterations && !stop; iIter++) {//Lloyd iterations

    std::cout << "Logging info..." << std::endl;
    std::cout << "LLoyd iteration number " << iIter << std::endl;
    
    // UPDATE WEIGHTS WITH FIXED SITES (using Karcher mean)

    //double alpha = 1.0, beta = 0.9;
    // double alpha = 1e-3, beta = 0.; // standard gradient descent 

    // // Nesterov acceleration
    // vector<double> phiWeightsY(nSites);

    std::vector<VertexData<double>> rhs;
    VertexData<double> normD;
    VertexData<double> rho(mesh, 0.0); // just for sanity check
    VertexData<double> normRHS(mesh, 0); // sum on all sites

    // New L-BFGS stuff
    options.epsWeights = max(1e-8, 1e-4 * sumUpdateNorm); // we want the weights to be converged, but only relative to the sites
    phiWeights = computePhiWeights(mesh, geom, phiWeights, vSolvers, options, measure, psMesh);

    cout << endl;
    // double minCellMass = *min_element(cellMasses.begin(), cellMasses.end());
    // double maxCellMass = *max_element(cellMasses.begin(), cellMasses.end());
    // cout << "Cell mass bounds: " << minCellMass << ", " << maxCellMass << endl;

    // if (iIter == 1)
    //   break;

    // UPDATE SITES WITH FIXED WEIGHTS (using Karcher mean)
    double energy = 0;

    options.nSubIterations = 1;
    for (size_t iSubIter = 0; iSubIter < options.nSubIterations; iSubIter++) {

      // Compute the normsalizer distribution
      rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
      normRHS = VertexData<double>(mesh, 0); // sum on all sites
      for (int i = 0; i < siteLocations.size(); i++)
        normRHS += rhs[i];
      normD = vSolver.scalarDiffuse(normRHS);

      energy = 0;
      vector<double> updateNorm(nSites);

      int totalTime = 0;

      // We want to do the mesh registration serially
      vector<VertexData<double>> fracDKarcher;
      for (int i = 0; i < omp_get_max_threads(); i++)
        fracDKarcher.emplace_back(VertexData<double>(mesh, 0));

      // // This doesn't work because it's an intrisic geometry fml
      // vector<Vector3> sites;
      // for (auto &site: siteLocations) {
      //   sites.push_back(site.interpolate(geom.vertexPositions));
      // }
      // polyscope::registerPointCloud("sites", sites);

      #pragma omp parallel for reduction(+:energy)
      for (size_t iSite = 0; iSite < nSites; iSite++) {

        int tid = omp_get_thread_num(); // thread ID

        SurfacePoint site = siteLocations[iSite];
        // H(site.face);

        // === Compute the nearest distribution
        fracDKarcher[tid] = vSolvers[tid].scalarDiffuse(rhs[iSite]);


        // // polyscope::registerPointCloud("sites", siteLocations);
        // psMesh.addVertexScalarQuantity("rhs", rhs[iSite]);
        // psMesh.addVertexScalarQuantity("fracDKarcher", fracDKarcher[tid]);
        // psMesh.addVertexScalarQuantity("normD", normD);
        // polyscope::show();

        for (Vertex v : mesh.vertices()) {
          fracDKarcher[tid][v] /= normD[v];
        }


        //WEIGHT THE DISTRIBUTION BY THE CURL MEASURE (vertex dual areas below when finding weight)
        for (Vertex v : mesh.vertices()){ 
          fracDKarcher[tid][v] *= measure[v];
        }

        //visualize the cells
        //result.cellEvolution[iSite][iIter] = thisFracDKarcher;

        // === Compute the log map
        // auto t1 = chrono::high_resolution_clock::now();
        //VertexData<Vector2> logmap = vSolvers[tid].computeLogMap(site, LogMapStrategy::AffineAdaptive);
        VertexData<Vector2> logmap = vSolvers[tid].computeLogMap(site);
        // auto t2 = chrono::high_resolution_clock::now();
        // totalTime += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();

        // Evaluate energy and gradient contribution
        Vector2 updateSum{0, 0};
        double updateWSum = 0.0;

        for (Vertex v : mesh.vertices()) {

          double weight = fracDKarcher[tid][v] * geom.vertexDualAreas[v];
          Vector2 logVal = logmap[v];
          double dist = logVal.norm();

          // H(fracDKarcher[tid][v]);
          // H(geom.vertexDualAreas[v]);

          updateSum += weight * logVal;
          updateWSum += weight;

          energy += dist * dist * weight; // TODO: use a per-site array to avoid data race
        }

        Vector2 update = updateSum / updateWSum;
        updateNorm[iSite] = update.norm();

        // Take a step
        TraceGeodesicResult traceResult = traceGeodesic(geom, site, options.stepSize * update);

        site = traceResult.endPoint;
        //viz the path it's taking 
        result.steps[iSite].push_back(site);
        //viz the distribution at every step
        result.stepSiteDistribution[iSite].push_back(fracDKarcher[tid]);

        siteLocations[iSite] = site;
      }

      H(energy);

      options.initialSites = siteLocations;

      sumUpdateNorm = 0;
      for (int i = 0; i < nSites; i++) sumUpdateNorm += updateNorm[i];

      // Convergence criterion: max update should be a fraction of edge length
      double maxUpdateNorm = *max_element(updateNorm.begin(), updateNorm.end());
      double eps_rel = 1e-2; // feel free to tweak this
      cout << "Lloyd convergence: " << maxUpdateNorm << " / " << eps_rel * meanEdgeLength << std::endl; 

      if (maxUpdateNorm < eps_rel * meanEdgeLength) {
        stop = true;
        break;
      }
    }//end of Karcher mean step to find mean of distribution
    cout << endl;

    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
    std::cout << "-------------------------" << std::endl;
  }

  // geom.unrequireVertexDualAreas();

  result.siteLocations = siteLocations;

  if (options.computeDistributions) {
    result.hasDistributions = true;

    // Compute the normalizer distribution
    std::vector<VertexData<double>> rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime); // make it sharp just for visu
    VertexData<double> normRHS(mesh, 0); // sum on all sites
    for (int i = 0; i < siteLocations.size(); i++)
      normRHS += rhs[i];
    VertexData<double> normD = vSolver.scalarDiffuse(normRHS);

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> thisFracD = vSolver.scalarDiffuse(rhs[iSite]);
      for (Vertex v : mesh.vertices()) thisFracD[v] /= normD[v];

      //WEIGHT THE DISTRIBUTION BY THE CURL MEASURE
      for (Vertex v : mesh.vertices()){ 
        // H(thisFracD[v]);
        // thisFracD[v] *= measure[v]; // density (this is what we want, i think)
        // thisFracD[v] *= measure[v] * geom.vertexDualAreas[v]; // mass
      }

      result.siteDistributions.push_back(thisFracD);
    }
  }
  return result;
}

// Project surface point on the isoline of a given target time value.
// `point` is modified in place.
// Note that this is not exact as the geodesic might not be orthogonal to the time function isolines,
// but it should be a pretty good guess.
void projectOnIsoline(SurfacePoint& point, double target, SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, alignOptions& options, FaceData<Vector2>& timeFunctionGrad) {

  HeatMethodDistanceSolver heatSolver(geom);
  VertexData<double> dist = heatSolver.computeDistance(point);

  double bestDist = DBL_MAX;
  for (Edge e : mesh.edges()) {
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    double t1 = options.timeFunction[v1], t2 = options.timeFunction[v2];
    //should also ask edge to be vertical
    if (fmin(t1, t2) < target && target < fmax(t1, t2)) { // edge is a candidate as it crosses the target isoline
      double eDist = (dist[v1] + dist[v2]) / 2;
      if (eDist < bestDist) {
        bestDist = eDist;
        double t = (target - t1) / (t2 - t1); // exact location along the edge
        point = SurfacePoint(e, t);
      }
    }
  }
}


void alignPointsOnIsolineFast(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, alignOptions& options, polyscope::SurfaceMesh &psMesh) {

  // Compute the gradient of the time function in the tangent plane of faces
  FaceData<Vector2> timeFunctionGrad = computeTimeFunctionFaceGradIntrinsic(geom, options.timeFunction);

  // For each pair of sings, project both to the midpoint of their time values
  for (auto &[s1,s2] : options.pairedSites) {  
    double t1 = s1.interpolate(options.timeFunction), t2 = s2.interpolate(options.timeFunction);
    projectOnIsoline(s1, (t1+t2)/2, mesh, geom, options, timeFunctionGrad);
    projectOnIsoline(s2, (t1+t2)/2, mesh, geom, options, timeFunctionGrad);
    ensure(abs(s1.interpolate(options.timeFunction) - s2.interpolate(options.timeFunction)) < 1e-6); // sanity check
  }
}

//computing equal weights using LBFGS
std::vector<double> computePhiWeights(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, std::vector<double>& initWeights, vector<VectorHeatMethodSolver>& vSolvers, VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh){
  
  //set the measure in Voronoi options (TO DO: specify this elsewhere) 
  options.measure = measure;

  //compute shortTime parameter (TO DO: specify this elsewhere)
  geom.requireEdgeLengths();
  //compute the diffusion time 
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  options.shortTime = options.tCoef * meanEdgeLength * meanEdgeLength;

  // Set points to start
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  int n = siteLocations.size();

  // Set up parameters
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = options.epsWeights;
  param.epsilon_rel = 1e-8;
  param.max_iterations = 1000;
  param.max_linesearch = 10000;
 
  // Create solver and function object
  LBFGSpp::LBFGSSolver<double> solver(param);
  F_OT fun(geom, mesh, options, siteLocations);
  fun.requireHeatKernel(vSolvers);

  Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(initWeights.data(), initWeights.size());
  Eigen::VectorXd grad(n);

  // x will be overwritten to be the best point found
  double fx = 0;
  int niter = solver.minimize(fun, x, fx);
  std::cout << niter << " L-BFGS iterations" << std::endl;
  return vector<double>(x.data(), x.data()+x.size());
}

} // namespace surface
} // namespace geometrycentral