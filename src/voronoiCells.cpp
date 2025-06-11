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
const VoronoiOptions defaultVoronoiOptions;

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

    //divide by factor of (4 * shortTime)
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

  // Set points to start
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  // if (siteLocations.empty()) {
  //   for (size_t i = 0; i < options.nSites; i++) {
  //     Face startF = mesh.face(randomIndex(mesh.nFaces()));
  //     double u = unitRand();
  //     Vector3 bCoord{u, 0.5 * (1.0 - u), 0.5 * (1.0 - u)};
  //     SurfacePoint fp{startF, bCoord};
  //     siteLocations.push_back(fp);
  //   }
  // }

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
  SurfacePoint dummySite = siteLocations[0];
  for (int i = 0; i < omp_get_max_threads(); i++) {
    vSolvers[i].scalarDiffuse(dummyRHS);
    vSolvers[i].computeLogMap(dummySite);
  }  

  size_t nSites = siteLocations.size();

  geom.requireVertexDualAreas();

  VoronoiResult result;
  result.steps.resize(nSites);

  //trying to find equal mass power cells 
  size_t descIter = 1000;
  // double stepSize = 1e-6;
  std::vector<double> phiWeights(nSites, 0.), oldPhiWeights(nSites);
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

    //double alpha = 2, beta = 0.2;
    double alpha = 2, beta = 0; // standard gradient descent 

    // Nesterov acceleration
    vector<double> phiWeightsY(nSites);

    std::vector<VertexData<double>> rhs;
    VertexData<double> normD;
    VertexData<double> rho(mesh, 0.0); // just for sanity check
    VertexData<double> normRHS(mesh, 0); // sum on all sites
    
    std::vector<double> cellMasses(nSites);

    //gradient descent to find the weights 
    for (size_t i = 0; i < descIter; i++){

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

      // double mean = 0;

      // We want to do the mesh registration serially
      vector<VertexData<double>> fracDWeights;
      for (int i = 0; i < omp_get_max_threads(); i++)
        fracDWeights.emplace_back(VertexData<double>(mesh, 0));

        

      double gradNorm = 0;

      #pragma omp parallel for reduction(+:gradNorm)
      for (size_t j = 0; j < nSites; j++){

        int tid = omp_get_thread_num(); // thread ID

        SurfacePoint site = siteLocations[j];
        fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);
        // cout << fracDWeights[tid] << endl;
        // fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);//scalarDiffuse takes something that's a mass and returns a density

        // Disabled because data race
        // for (Vertex v : mesh.vertices()) {
        //   rho[v] += thisFracDWeights[v] / normD[v];
        // }

        //normalize and weigh the distibution by the curl measure 
        //also integrate lol
        double updateWSum = 0.0;
        for (Vertex v : mesh.vertices()) {
          fracDWeights[tid][v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
          //thisFracDWeights[v] *= measure[v]/normD[v];
          updateWSum += fracDWeights[tid][v];
        }

        gradNorm += (desiredMass-updateWSum) * (desiredMass-updateWSum);

        // std::cout << "updateWSum = " << updateWSum << std::endl;
        //newWeight = phiWeights[j] + (((1/iIter) * shortTime) * (desiredMass - updateWSum));//do a time-decaying step size here

        phiWeights[j] = phiWeightsY[j] + alpha * (desiredMass - updateWSum);

        // phiWeights[j] += 1e-0 * (desiredMass - updateWSum);//don't do a time-decaying step size

        // // Nesterov (need to fix)
        // double y = phiWeights[j] + beta * (phiWeights[j] - oldPhiWeights[j]);
        // oldPhiWeights[j] = phiWeights[j]; // update old
        // phiWeights[j] += 

        // double z = beta * phiWeights[j] + 1e-0 * (desiredMass - updateWSum);
        // phiWeights[j] = phiWeights[j] + alpha * z;

        // newWeight = phiWeights[j] + ((shortTime) * (desiredMass - updateWSum));//don't do a time-decaying step size
        // mean += newWeight; 

        cellMasses[j] = updateWSum;
      }

      gradNorm = sqrt(gradNorm);
      // H(gradNorm);
      std::cout << "gradNorm: " << gradNorm << "\t\r" << std::flush;

      psMesh.addVertexScalarQuantity("rho", rho);

      if (gradNorm < 1e-2 * sumUpdateNorm) {
        cout << endl << "Performed " << i+1 << " iterations for weights update.";
        break;
      }

      //ensure that the average of all the weights is 0
      // mean /= nSites;
      // for (size_t k = 0; k < nSites; k++){
      //   newPhiWeights[k] -= mean;
      // }
      // phiWeights = newPhiWeights;
    }

    cout << endl;


    // std::cout << "weights after gradient descent: " << std::endl;
    // for (int i = 0; i < phiWeights.size(); i++){
    //   std::cout << "weight at site " << i << ": " << phiWeights[i] << std::endl;
    // }

    double minCellMass = *min_element(cellMasses.begin(), cellMasses.end());
    double maxCellMass = *max_element(cellMasses.begin(), cellMasses.end());
    cout << "Cell mass bounds: " << minCellMass << ", " << maxCellMass << endl;

    // UPDATE SITES WITH FIXED WEIGHTS (using Karcher mean)
    double energy = 0;
    // std::vector<SurfacePoint> newSiteLocations;

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

        //updateSum /= updateWSum;
        Vector2 update = updateSum / updateWSum;
        updateNorm[iSite] = update.norm();

        // Take a step
        TraceGeodesicResult traceResult = traceGeodesic(geom, site, options.stepSize * update);

        site = traceResult.endPoint;
        //viz the path it's taking 
        result.steps[iSite].push_back(site);

        siteLocations[iSite] = site;
      }

      H(energy);

      // H(totalTime);
      sumUpdateNorm = 0;
      for (int i = 0; i < nSites; i++) sumUpdateNorm += updateNorm[i];
      cout << "sumUpdateNorm = " << sumUpdateNorm << "\t\r" << flush;

      if (sumUpdateNorm < options.eps) {
        stop = true;
        break;
      }

      // newSiteLocations.push_back(site);
    }
    cout << endl;

    // siteLocations = newSiteLocations;
    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
    std::cout << "-------------------------" << std::endl;
  }

  geom.unrequireVertexDualAreas();

  result.siteLocations = siteLocations;

  if (options.computeDistributions) {
    result.hasDistributions = true;

    // Compute the normalizer distribution
    std::vector<VertexData<double>> rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
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
        thisFracD[v] *= measure[v] * geom.vertexDualAreas[v];
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

  // Binary search in direction v to find the target time value
  double lower = -1, upper = +1.1; // is it enough? For some reason if mid=0, traceGeodesic fails
  while(upper - lower > 1e-6) {
    double mid = (lower+upper) / 2;
    TraceGeodesicResult traceResult = traceGeodesic(geom, initPoint, mid*v);
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

VoronoiResult alignPointsOnIsoline(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                    alignOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh){
  
  
  //mostly the same code as above (with some edits :P)

  // Define one vector heat method solver per thread
  vector<VectorHeatMethodSolver> vSolvers;
  for (int i = 0; i < omp_get_max_threads(); i++)
    vSolvers.emplace_back(geom, options.tCoef);
  VectorHeatMethodSolver& vSolver = vSolvers[0]; // temporarily, for the serial code


  // This computes e^(φ_j/4t) δ(x_j) for every site j
  auto computeRHSWithWeights = [&](const std::vector<SurfacePoint>& points, std::vector<double>& weights, double shortTime) -> std::vector<VertexData<double>> {

    std::vector<VertexData<double>> rhs;
    double maxWeight = *std::max_element(weights.begin(), weights.end());

    //divide by factor of (4 * shortTime)
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

  // Set points to start
  std::vector<SurfacePoint> siteLocations;
  //make a map from positive site location index to time function values on it's oppositely signed pair 
  std::map<int, double> mapToPairedIsoVal;
  //if using the positive curl signal, set the initial sets to the first entry in the pair 
  if (options.usingPosCurl){
    for (int i = 0; i < options.pairedSites.size(); i++){
      auto p = options.pairedSites[i];
      siteLocations.push_back(p.first);
      mapToPairedIsoVal[i] = p.second.interpolate(options.timeFunction);
    }
  }

  // For some reason, running a dummy scalarDiffuse / computeLogMap first
  // fixes thread issues lol. Many thanks to
  // https://github.com/nmwsharp/geometry-central/issues/108
  VertexData<double> dummyRHS(mesh, 42.0);
  SurfacePoint dummySite = siteLocations[0];
  for (int i = 0; i < omp_get_max_threads(); i++) {
    vSolvers[i].scalarDiffuse(dummyRHS);
    vSolvers[i].computeLogMap(dummySite);
  }  

  size_t nSites = siteLocations.size();

  geom.requireVertexDualAreas();

  VoronoiResult result;
  result.steps.resize(nSites);

  //trying to find equal mass power cells 
  size_t descIter = 1000;
  // double stepSize = 1e-6;
  std::vector<double> phiWeights(nSites, 0.), oldPhiWeights(nSites);
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

  //step size for aligning points 
  double lambda = options.lambda;

  // == Iterations
  for (size_t iIter = 0; iIter < options.iterations; iIter++) {//Lloyd iterations

    std::cout << "Logging info..." << std::endl;
    std::cout << "LLoyd iteration number " << iIter << std::endl;
    
    // UPDATE WEIGHTS WITH FIXED SITES (using Karcher mean)

    double alpha = 2, beta = 0.2;
    // double alpha = 2, beta = 0; // standard gradient descent 

    // Nesterov acceleration
    vector<double> phiWeightsY(nSites);

    std::vector<VertexData<double>> rhs;
    VertexData<double> normD;
    VertexData<double> rho(mesh, 0.0); // just for sanity check
    VertexData<double> normRHS(mesh, 0); // sum on all sites
    
    std::vector<double> cellMasses(nSites);

    //gradient descent to find the weights 
    for (size_t i = 0; i < descIter; i++){

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

      // double mean = 0;

      // We want to do the mesh registration serially
      vector<VertexData<double>> fracDWeights;
      for (int i = 0; i < omp_get_max_threads(); i++)
        fracDWeights.emplace_back(VertexData<double>(mesh, 0));

        

      double gradNorm = 0;

      #pragma omp parallel for reduction(+:gradNorm)
      for (size_t j = 0; j < nSites; j++){

        int tid = omp_get_thread_num(); // thread ID

        SurfacePoint site = siteLocations[j];
        fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);
        // cout << fracDWeights[tid] << endl;
        // fracDWeights[tid] = vSolvers[tid].scalarDiffuse(rhs[j]);//scalarDiffuse takes something that's a mass and returns a density

        // Disabled because data race
        // for (Vertex v : mesh.vertices()) {
        //   rho[v] += thisFracDWeights[v] / normD[v];
        // }

        //normalize and weigh the distibution by the curl measure 
        //also integrate lol
        double updateWSum = 0.0;
        for (Vertex v : mesh.vertices()) {
          fracDWeights[tid][v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
          //thisFracDWeights[v] *= measure[v]/normD[v];
          updateWSum += fracDWeights[tid][v];
        }

        gradNorm += (desiredMass-updateWSum) * (desiredMass-updateWSum);

        // std::cout << "updateWSum = " << updateWSum << std::endl;
        //newWeight = phiWeights[j] + (((1/iIter) * shortTime) * (desiredMass - updateWSum));//do a time-decaying step size here

        phiWeights[j] = phiWeightsY[j] + alpha * (desiredMass - updateWSum);

        // phiWeights[j] += 1e-0 * (desiredMass - updateWSum);//don't do a time-decaying step size

        // // Nesterov (need to fix)
        // double y = phiWeights[j] + beta * (phiWeights[j] - oldPhiWeights[j]);
        // oldPhiWeights[j] = phiWeights[j]; // update old
        // phiWeights[j] += 

        // double z = beta * phiWeights[j] + 1e-0 * (desiredMass - updateWSum);
        // phiWeights[j] = phiWeights[j] + alpha * z;

        // newWeight = phiWeights[j] + ((shortTime) * (desiredMass - updateWSum));//don't do a time-decaying step size
        // mean += newWeight; 

        cellMasses[j] = updateWSum;
      }

      gradNorm = sqrt(gradNorm);
      // H(gradNorm);
      std::cout << "gradNorm: " << gradNorm << "\t\r" << std::flush;


      if (gradNorm < 1e-3)
        break;

      psMesh.addVertexScalarQuantity("rho", rho);
      
      //ensure that the average of all the weights is 0
      // mean /= nSites;
      // for (size_t k = 0; k < nSites; k++){
      //   newPhiWeights[k] -= mean;
      // }
      // phiWeights = newPhiWeights;
    }

    cout << endl;


    // std::cout << "weights after gradient descent: " << std::endl;
    // for (int i = 0; i < phiWeights.size(); i++){
    //   std::cout << "weight at site " << i << ": " << phiWeights[i] << std::endl;
    // }

    double minCellMass = *min_element(cellMasses.begin(), cellMasses.end());
    double maxCellMass = *max_element(cellMasses.begin(), cellMasses.end());
    cout << "Cell mass bounds: " << minCellMass << ", " << maxCellMass << endl;

    // UPDATE SITES WITH FIXED WEIGHTS (using Karcher mean)
    double energy = 0;
    // std::vector<SurfacePoint> newSiteLocations;

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


      #pragma omp parallel for
      for (size_t iSite = 0; iSite < nSites; iSite++) {
        
        Vector2 alignmentTerm;//alignment term to move singularities near the same isoline
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

          // energy += dist * dist * weight; // TODO: use a per-site array to avoid data race
        }

        //align points on the same isoline of the time function 
        if (options.usingPosCurl){
          double posTimeFuncVal = siteLocations[iSite].interpolate(options.timeFunction);
          double negTimeFuncVal = mapToPairedIsoVal[iSite];

          // approximate gradient of time function pulled back to tangent plane
          double hi = options.timeFunction[siteLocations[iSite].inSomeFace().face.halfedge().vertex()];
          double hj = options.timeFunction[siteLocations[iSite].inSomeFace().face.halfedge().next().vertex()];
          double hk = options.timeFunction[siteLocations[iSite].inSomeFace().face.halfedge().next().next().vertex()];
        
          Vector2 xi = logmap[siteLocations[iSite].inSomeFace().face.halfedge().vertex()];
          Vector2 xj = logmap[siteLocations[iSite].inSomeFace().face.halfedge().next().vertex()];
          Vector2 xk = logmap[siteLocations[iSite].inSomeFace().face.halfedge().next().next().vertex()]; 
        
          double area = fabs(cross(xj-xi,xk-xi));
        
          Vector2 gradh = ((hj-hi) * (xi - xk).rotate90()) / area + ((hk - hi) * (xj-xi).rotate90()) / area;

          // std::cout << "lambda = " << options.lambda << std::endl;
          // std::cout << "pos time func val = " << posTimeFuncVal << std::endl;
          // std::cout << "neg time func val = " << negTimeFuncVal << std::endl;
          // std::cout << "posTimeFuncVal - negTimeFuncVal = "  << (posTimeFuncVal - negTimeFuncVal) << std::endl;
          // std::cout << "gradh = " << gradh << std::endl;
        
          alignmentTerm = -(posTimeFuncVal - negTimeFuncVal) * gradh;
        }

        cout << "lambda = " << options.lambda << endl;
        cout << "alignmentTerm = " << alignmentTerm << endl;

        //add the alignment term to the update
        Vector2 update = (1 - options.lambda) * (updateSum / updateWSum) +  (options.lambda) * alignmentTerm;

        updateNorm[iSite] = update.norm();

        // Take a step
        TraceGeodesicResult traceResult = traceGeodesic(geom, site, options.stepSize * update);

        site = traceResult.endPoint;
        //viz the path it's taking 
        result.steps[iSite].push_back(site);

        siteLocations[iSite] = site;
      }
      // H(totalTime);
      double sumUpdateNorm = 0;
      for (int i = 0; i < nSites; i++) sumUpdateNorm += updateNorm[i];
      cout << "sumUpdateNorm = " << sumUpdateNorm << "\t\r" << flush;

      // newSiteLocations.push_back(site);
    }
    cout << endl;

    // siteLocations = newSiteLocations;
    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
    std::cout << "-------------------------" << std::endl;
  }

  geom.unrequireVertexDualAreas();

  result.siteLocations = siteLocations;

  if (options.computeDistributions) {
    result.hasDistributions = true;

    // Compute the normalizer distribution
    std::vector<VertexData<double>> rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
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
        thisFracD[v] *= measure[v] * geom.vertexDualAreas[v];
      }

      result.siteDistributions.push_back(thisFracD);
    }
  }
  return result;
}

} // namespace surface
} // namespace geometrycentral