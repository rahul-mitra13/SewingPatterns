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
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  //use a random seed 
  std::mt19937 rng(options.seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::uniform_int_distribution<int> pointTypeDist(0, 2); // 0: vertex, 1: edge, 2: face
  for (size_t i = 0; i < options.nSites; ++i) {
      int type = pointTypeDist(rng);

      if (type == 0) {
        // Random vertex
        Vertex v = mesh.vertex(rng() % mesh.nVertices());
        siteLocations.emplace_back(v);
      } else if (type == 1) {
        // Random edge interior
        Edge e = mesh.edge(rng() % mesh.nEdges());
        double t = uniform01(rng); // interpolation along the edge
        siteLocations.emplace_back(e, t);
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
        siteLocations.push_back(pt);
      }
  }

  //debug on the square 
  //f2592 and f6688 are next to each other 
  // siteLocations.push_back(SurfacePoint(mesh.face(2592), Vector3{0.25, 0.25, 0.5}));
  // siteLocations.push_back(SurfacePoint(mesh.face(6688), Vector3{1./3., 1./3., 1./3.}));

  size_t nSites = siteLocations.size();
  VoronoiResult result;
  result.steps.resize(nSites);
  result.stepSiteDistribution.resize(nSites);
  result.initialSites = siteLocations;

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
  std::cout << "shortTime outside function = " << shortTime << std::endl;

  //compute mass per cell 
  double desiredMass = 0.; 
  for (Vertex v : mesh.vertices()){
    //std::cout << "dual area at vertex " << v << " " << geom.vertexDualAreas[v] << std::endl;
    desiredMass += geom.vertexDualAreas[v] * measure[v];
  }
  desiredMass /= nSites;
  std::cout << "desiredMass per cell = " << desiredMass << std::endl;

  // H(SUITESPARSE_USE_OPENMP);

  // == Iterations
  bool stop = false; // flag if stopping criterion on sumUpdateNorm is reached
  double sumUpdateNorm = 1; // just a starting point to determine a tolerance on weights grad norm
  for (size_t iIter = 0; iIter < options.iterations && !stop; iIter++) {//Lloyd iterations

    std::cout << "Logging info..." << std::endl;
    std::cout << "LLoyd iteration number " << iIter << std::endl;
    
    // UPDATE WEIGHTS WITH FIXED SITES (using Karcher mean)

    //double alpha = 1.0, beta = 0.9;
    double alpha = 1e-3, beta = 0.; // standard gradient descent 

    // Nesterov acceleration
    vector<double> phiWeightsY(nSites);

    std::vector<VertexData<double>> rhs;
    VertexData<double> normD;
    VertexData<double> rho(mesh, 0.0); // just for sanity check
    VertexData<double> normRHS(mesh, 0); // sum on all sites
    double gradNorm = 0;

    // // New L-BFGS stuff
    // Eigen::VectorXd phi = computePhiWeights(mesh, geom, options, measure, psMesh);
    // for (int iSite = 0; iSite < options.nSites; iSite++)
    //   phiWeights[iSite] = phi[iSite];

    //gradient descent to find the weights 
    for (size_t i = 0; i < options.descIter; i++){

      if (options.useLineSearch){
        //Line Search
        double oldEnergy;
        oldEnergy = proxyEnergy(cellMasses, iIter);
        double stepSize = alpha;
        for (int lineSearchIter = 0; lineSearchIter < 8; lineSearchIter++){
          //try taking a step
          // Compute the y_k of Nesterov's AGD
          for (int j = 0; j < nSites; j++) {
            phiWeightsY[j] = phiWeights[j] + beta * (phiWeights[j] - oldPhiWeights[j]);
            oldPhiWeights[j] = phiWeights[j]; // don't need old anymore, we can update
          }

          // Now compute gradient of y_k
          rhs = computeRHSWithWeights(siteLocations, phiWeightsY, shortTime);

          // Compute the normalizer distribution
          normRHS.fill(0);
          for (int i = 0; i < siteLocations.size(); i++)
            normRHS += rhs[i];
          normD = vSolver.scalarDiffuse(normRHS);

          // We want to do the mesh registration serially
          vector<VertexData<double>> fracDWeights;
          for (int i = 0; i < omp_get_max_threads(); i++)
            fracDWeights.emplace_back(VertexData<double>(mesh, 0));

          #pragma omp parallel for reduction(+:gradNorm)
          for (size_t j = 0; j < nSites; j++){

            int tid = omp_get_thread_num(); // thread ID

            SurfacePoint site = siteLocations[j];
            //scalarDiffuse takes something that's a mass and returns a density
            fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);
        
            //normalize and weigh the distibution by the curl measure 
            double updateWSum = 0.0;
            for (Vertex v : mesh.vertices()) {
              fracDWeights[tid][v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
              updateWSum += fracDWeights[tid][v];
            }

            gradNorm += (desiredMass-updateWSum) * (desiredMass-updateWSum);
            phiWeights[j] = phiWeightsY[j] + stepSize * (desiredMass - updateWSum);
            cellMasses[j] = updateWSum;
          }
          double newEnergy = proxyEnergy(cellMasses, false);
          // Accept step if good
          if (newEnergy < oldEnergy) {
            break;
          }
          // Otherwise decrease step size and repeat
          stepSize *= 0.5;

        }//end of line search 
      }
      else{
        // Compute the y_k of Nesterov's AGD
        for (int j = 0; j < nSites; j++) {
          phiWeightsY[j] = phiWeights[j] + beta * (phiWeights[j] - oldPhiWeights[j]);
          oldPhiWeights[j] = phiWeights[j]; // don't need old anymore, we can update
        }

        // Now compute gradient of y_k
        rhs = computeRHSWithWeights(siteLocations, phiWeightsY, shortTime);

        // Compute the normalizer distribution
        normRHS.fill(0);
        for (int i = 0; i < siteLocations.size(); i++)
          normRHS += rhs[i];
        normD = vSolver.scalarDiffuse(normRHS);

        // We want to do the mesh registration serially
        vector<VertexData<double>> fracDWeights;
        for (int i = 0; i < omp_get_max_threads(); i++)
          fracDWeights.emplace_back(VertexData<double>(mesh, 0));

        #pragma omp parallel for reduction(+:gradNorm)
        for (size_t j = 0; j < nSites; j++){

          int tid = omp_get_thread_num(); // thread ID

          SurfacePoint site = siteLocations[j];
          //scalarDiffuse takes something that's a mass and returns a density
          fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);
        
          //normalize and weigh the distibution by the curl measure 
          double updateWSum = 0.0;
          for (Vertex v : mesh.vertices()) {
            fracDWeights[tid][v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
            updateWSum += fracDWeights[tid][v];
          }

          gradNorm += (desiredMass-updateWSum) * (desiredMass-updateWSum);
          phiWeights[j] = phiWeightsY[j] + alpha * (desiredMass - updateWSum);
          cellMasses[j] = updateWSum;
        }
      }

      gradNorm = sqrt(gradNorm);
      std::cout << "gradNorm: " << gradNorm << "\t\r" << std::flush;

      psMesh.addVertexScalarQuantity("rho", rho);

      // if (gradNorm < 1e-2 * sumUpdateNorm) {
      if (gradNorm < 1e-9) {
        cout << endl << "Performed " << i+1 << " iterations for weights update.";
        break;
      }
    }//end of gradient descent

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


      #pragma omp parallel for reduction(+:energy)
      for (size_t iSite = 0; iSite < nSites; iSite++) {

        int tid = omp_get_thread_num(); // thread ID

        SurfacePoint site = siteLocations[iSite];

        // === Compute the nearest distribution
        fracDKarcher[tid] = vSolvers[tid].scalarDiffuse(rhs[iSite]);

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

      // H(totalTime);
      sumUpdateNorm = 0;
      for (int i = 0; i < nSites; i++) sumUpdateNorm += updateNorm[i];
      cout << "sumUpdateNorm = " << sumUpdateNorm << "\t\r" << flush;

      if (sumUpdateNorm < options.eps) {
        stop = true;
        break;
      }
    }//end of Karcher mean step to find mean of distribution
    cout << endl;

    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
    std::cout << "-------------------------" << std::endl;
  }

  geom.unrequireVertexDualAreas();

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

  SurfacePoint initPoint = point; // keep a copy
  Vector2 v = timeFunctionGrad[initPoint.face]; // search direction
  v = v.normalize();

  geom.requireShapeLengthScale(); // used to scale the search direction

  // Binary search in direction v to find the target time value
  double lower = -1, upper = +1.1; // is it enough? For some reason if mid=0, traceGeodesic fails
  while(upper - lower > 1e-6) {
    double mid = (lower+upper) / 2;
    TraceGeodesicResult traceResult = traceGeodesic(geom, initPoint, mid*v*geom.shapeLengthScale);
    point = traceResult.endPoint;
    double t = point.interpolate(options.timeFunction);
    if (t > target)
      upper = mid;
    else 
      lower = mid;
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
  }
}

//computing equal weights using LBFGS
Eigen::VectorXd computePhiWeights(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                        VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh){
  
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

  // Define one vector heat method solver per thread
  vector<VectorHeatMethodSolver> vSolvers;
  for (int i = 0; i < omp_get_max_threads(); i++)
    vSolvers.emplace_back(geom, options.tCoef);

  // For some reason, running a dummy scalarDiffuse / computeLogMap first
  // fixes thread issues lol. Many thanks to
  // https://github.com/nmwsharp/geometry-central/issues/108
  VertexData<double> dummyRHS(mesh, 42.0);
  SurfacePoint dummySite = mesh.vertex(0);
  for (int i = 0; i < omp_get_max_threads(); i++) {
    vSolvers[i].scalarDiffuse(dummyRHS);
    vSolvers[i].computeLogMap(dummySite);
  }
  
  
  // Set points to start
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  //use a random seed 
  std::mt19937 rng(options.seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::uniform_int_distribution<int> pointTypeDist(0, 2); // 0: vertex, 1: edge, 2: face
  for (size_t i = 0; i < options.nSites; ++i) {
      int type = pointTypeDist(rng);

      if (type == 0) {
        // Random vertex
        Vertex v = mesh.vertex(rng() % mesh.nVertices());
        siteLocations.emplace_back(v);
      } else if (type == 1) {
        // Random edge interior
        Edge e = mesh.edge(rng() % mesh.nEdges());
        double t = uniform01(rng); // interpolation along the edge
        siteLocations.emplace_back(e, t);
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
        siteLocations.push_back(pt);
      }
  }

  int n = siteLocations.size();

  // Set up parameters
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-5;
  param.max_iterations = 1000;
  param.max_linesearch = 10000;
 
  // Create solver and function object
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchMoreThuente> solver(param);
  F_OT fun(geom, mesh, options, siteLocations);
  fun.requireHeatKernel(vSolvers);
  fun.requireLogMap(vSolvers);

  // Initial guess
  // Eigen::VectorXd x(options.nSites);
  // for (int i = 0; i < options.nSites; i++){
  //   x(i) = uniform01(rng);
  // }
  Eigen::VectorXd x = Eigen::VectorXd::Ones(options.nSites);

  // Check gradient
  std::vector<double> eps {1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};
  for (auto h : eps)
    H(fun.checkGrad(x, h));
  polyscope::show();

  // x will be overwritten to be the best point found
  double fx;
  int niter = solver.minimize(fun, x, fx);
  x = -x; // max -> min
 
  std::cout << niter << " iterations" << std::endl;
  // std::cout << "x = \n" << x.transpose() << std::endl;
  // std::cout << "f(x) = " << fx << std::endl;
  
  return x;
}

} // namespace surface
} // namespace geometrycentral