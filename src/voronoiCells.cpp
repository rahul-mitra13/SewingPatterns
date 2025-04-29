#include "voronoiCells.h"
#include "math.h"
#include "float.h"

namespace geometrycentral {
namespace surface {

namespace {
const bool VORONOI_PRINT = false;
}

// The default trace options
const VoronoiOptions defaultVoronoiOptions;

VoronoiResult computeGeodesicCentroidalVoronoiTessellationWithWeights(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options, VertexData<double>& measure) {

  if (options.useDelaunay) {
    std::unique_ptr<SignpostIntrinsicTriangulation> intTri(new SignpostIntrinsicTriangulation(mesh, geom));
    intTri->flipToDelaunay();

    // Translate initial sites to intrinsic triangulation if any were given
    for (size_t iS = 0; iS < options.initialSites.size(); iS++) {
      options.initialSites[iS] = intTri->equivalentPointOnIntrinsic(options.initialSites[iS]);
    }

    // Get solutions on intrinsic triangulation
    options.useDelaunay = false;
    VoronoiResult result = computeGeodesicCentroidalVoronoiTessellationWithWeights(*intTri->intrinsicMesh, *intTri, options, measure);

    // Translate solutions back to original triangulation
    for (size_t iS = 0; iS < result.siteLocations.size(); iS++) {
      result.siteLocations[iS] = intTri->equivalentPointOnInput(result.siteLocations[iS]);
    }

    return result;
  }

  VectorHeatMethodSolver vSolver(geom, options.tCoef);

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

  auto computeRHSWithWeights = [&](const std::vector<SurfacePoint>& points, std::vector<double>& weights, double shortTime) -> VertexData<double> {
    VertexData<double> rhs(mesh, 0);

    //find the max of the weights 
    double maxWeight = DBL_MIN;
    for (size_t i = 0; i < weights.size(); i++){
      if (weights[i] > maxWeight) maxWeight = weights[i];
    }
    // //update the weights subtracting off the max weights so that the exponent doesn't blow up
    // for (size_t i = 0; i < weights.size(); i++){
    //   weights[i] -= maxWeight;
    // }

    //divide by factor of (4 * shortTime)
    for (size_t i = 0; i < points.size(); i++){
      SurfacePoint facePoint = points[i].inSomeFace();
      Halfedge he = facePoint.face.halfedge();
      rhs[he.vertex()] += exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.x;
      rhs[he.next().vertex()] +=  exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.y;
      rhs[he.next().next().vertex()] += exp((weights[i] - maxWeight)/(4 * shortTime)) * facePoint.faceCoords.z;
    }

    return rhs;
  };

  // Set points to start
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  if (siteLocations.empty()) {
    for (size_t i = 0; i < options.nSites; i++) {
      Face startF = mesh.face(randomIndex(mesh.nFaces()));
      double u = unitRand();
      Vector3 bCoord{u, 0.5 * (1.0 - u), 0.5 * (1.0 - u)};
      SurfacePoint fp{startF, bCoord};
      siteLocations.push_back(fp);
    }
  }

  size_t nSites = siteLocations.size();

  geom.requireVertexDualAreas();

  VoronoiResult result;
  result.steps.resize(nSites);

  //trying to find equal mass power cells 
  size_t descIter = 1;
  double stepSize = 1e-6;
  std::vector<double> phiWeights(nSites, 0.);
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

  //for visualization of the sites
  result.cellEvolution.resize(nSites);
  for (int i = 0; i < nSites; i++){
    result.cellEvolution[i].resize(options.iterations);
  }

  // == Iterations
  for (size_t iIter = 0; iIter < options.iterations; iIter++) {//Lloyd iterations

    std::cout << "Logging info..." << std::endl;
    std::cout << "LLoyd iteration number " << iIter << std::endl;
    
    // Compute the normalizer distribution
    VertexData<double> normRHS = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
    VertexData<double> normD = vSolver.scalarDiffuse(normRHS);
    
    //gradient descent to find the weights 
    for (size_t i = 0; i < descIter; i++){
      double mean = 0.;
      std::vector<double> newPhiWeights(nSites, 0.);
      for (size_t j = 0; j < nSites; j++){//don't we want to move this loop out of the gradient descent loop? 

        double newWeight = 0.; 
        SurfacePoint site = siteLocations[j];
        VertexData<double> unitRHSWeights = computeRHSWithWeights({site}, phiWeights, shortTime);
        VertexData<double> thisFracDWeights = vSolver.scalarDiffuse(unitRHSWeights);//scalarDiffuse takes something that's a mass and returns a density

        //std::cout << "At site " << j << " weight: " << (phiWeights[j] / (4 * shortTime)) << std::endl;
        //visualize thisFracDWeights 
        result.cellEvolution[j][iIter] = thisFracDWeights;

        double updateWSum = 0.0;
        //normalize and weigh the distibution by the curl measure 
        //also integrate lol
        for (Vertex v : mesh.vertices()){
          thisFracDWeights[v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
          //thisFracDWeights[v] *= measure[v]/normD[v];
          updateWSum += thisFracDWeights[v];
        }

        std::cout << "updateWSum = " << updateWSum << std::endl;
        //newWeight = phiWeights[j] + (((1/iIter) * shortTime) * (desiredMass - updateWSum));//do a time-decaying step size here
        newWeight = phiWeights[j] + ((shortTime) * (desiredMass - updateWSum));//don't do a time-decaying step size
        mean += newWeight; 
        //update the new phi weights
        newPhiWeights[j] = newWeight;
      }
      //ensure that the average of all the weights is 0
      // mean /= nSites;
      // for (size_t k = 0; k < nSites; k++){
      //   newPhiWeights[k] -= mean;
      // }
      phiWeights = newPhiWeights;
    }


    std::cout << "weights after gradient descent: " << std::endl;
    for (int i = 0; i < phiWeights.size(); i++){
      std::cout << "weight at site " << i << ": " << phiWeights[i] << std::endl;
    }


    normRHS = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
    normD = vSolver.scalarDiffuse(normRHS);

    double energy = 0;
    std::vector<SurfacePoint> newSiteLocations;

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> unitRHSKarcher = computeRHSWithWeights({site}, phiWeights, shortTime);
      VertexData<double> thisFracDKarcher = vSolver.scalarDiffuse(unitRHSKarcher);

      for (Vertex v : mesh.vertices()) {
        thisFracDKarcher[v] /= normD[v];
      }

      //WEIGHT THE DISTRIBUTION BY THE CURL MEASURE (vertex dual areas below when finding weight)
      for (Vertex v : mesh.vertices()){ 
        thisFracDKarcher[v] *= measure[v];
      }

      //visualize the cells
      //result.cellEvolution[iSite][iIter] = thisFracDKarcher;


      for (size_t iSubIter = 0; iSubIter < options.nSubIterations; iSubIter++) {

        // === Compute the log map
        VertexData<Vector2> logmap = vSolver.computeLogMap(site);

        // Evaluate energy and gradient contribution
        Vector2 updateSum{0, 0};
        double updateWSum = 0.0;

        for (Vertex v : mesh.vertices()) {

          double weight = thisFracDKarcher[v] * geom.vertexDualAreas[v];
          Vector2 logVal = logmap[v];
          double dist = logVal.norm();

          updateSum += weight * logVal;
          updateWSum += weight;

          energy += dist * dist * weight;
        }
        //updateSum /= updateWSum;
        Vector2 update = updateSum / updateWSum;

        // Take a step
        TraceGeodesicResult traceResult = traceGeodesic(geom, site, options.stepSize * update);
        site = traceResult.endPoint;
        //viz the path it's taking 
        result.steps[iSite].push_back(site);
      }

      newSiteLocations.push_back(site);
    }

    siteLocations = newSiteLocations;
    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
    std::cout << "-------------------------" << std::endl;
  }

  geom.unrequireVertexDualAreas();

  result.siteLocations = siteLocations;

  if (options.computeDistributions) {
    result.hasDistributions = true;

    // Compute the normalizer distribution
    VertexData<double> normRHS = computeRHS(siteLocations);
    VertexData<double> normD = vSolver.scalarDiffuse(normRHS);

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> unitRHS = computeRHS({site});
      VertexData<double> thisFracD = vSolver.scalarDiffuse(unitRHS);
      for (Vertex v : mesh.vertices()) thisFracD[v] /= normD[v];

      //WEIGHT THE DISTRIBUTION BY THE CURL MEASURE
      for (Vertex v : mesh.vertices()){ 
        thisFracD[v] *= measure[v];
      }

      result.siteDistributions.push_back(thisFracD);
    }
  }
  return result;
}
} // namespace surface
} // namespace geometrycentral
