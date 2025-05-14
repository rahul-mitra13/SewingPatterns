#include "voronoiCells.h"
#include "math.h"
#include "float.h"
#include "helpers.h"

namespace geometrycentral {
namespace surface {

namespace {
const bool VORONOI_PRINT = false;
}

// The default trace options
const VoronoiOptions defaultVoronoiOptions;

VoronoiResult computeGeodesicCentroidalVoronoiTessellationWithWeights(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options, VertexData<double>& measure, polyscope::SurfaceMesh &psMesh) {

  if (options.useDelaunay) {
    std::unique_ptr<SignpostIntrinsicTriangulation> intTri(new SignpostIntrinsicTriangulation(mesh, geom));
    intTri->flipToDelaunay();

    // Translate initial sites to intrinsic triangulation if any were given
    for (size_t iS = 0; iS < options.initialSites.size(); iS++) {
      options.initialSites[iS] = intTri->equivalentPointOnIntrinsic(options.initialSites[iS]);
    }

    // Get solutions on intrinsic triangulation
    options.useDelaunay = false;
    VoronoiResult result = computeGeodesicCentroidalVoronoiTessellationWithWeights(*intTri->intrinsicMesh, *intTri, options, measure, psMesh);

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
  size_t descIter = 100;
  // double stepSize = 1e-6;
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

  

  // == Iterations
  for (size_t iIter = 0; iIter < options.iterations; iIter++) {//Lloyd iterations

    std::cout << "Logging info..." << std::endl;
    std::cout << "LLoyd iteration number " << iIter << std::endl;
    
    std::vector<VertexData<double>> rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
    VertexData<double> normD;
    VertexData<double> rho(mesh, 0.0); // just for sanity check
    VertexData<double> normRHS(mesh, 0); // sum on all sites
    
    std::vector<double> cellMasses(nSites);

    //gradient descent to find the weights 
    for (size_t i = 0; i < descIter; i++){

      rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);

      // Compute the normalizer distribution
      normRHS.fill(0);
      for (int i = 0; i < siteLocations.size(); i++)
        normRHS += rhs[i];
      normD = vSolver.scalarDiffuse(normRHS);

      // double mean = 0;

      double gradNorm = 0;

      for (size_t j = 0; j < nSites; j++){

        SurfacePoint site = siteLocations[j];
        VertexData<double> thisFracDWeights = vSolver.scalarDiffuse(rhs[j]);//scalarDiffuse takes something that's a mass and returns a density

        for (Vertex v : mesh.vertices()) {
          rho[v] += thisFracDWeights[v] / normD[v];
        }

        //normalize and weigh the distibution by the curl measure 
        //also integrate lol
        double updateWSum = 0.0;
        for (Vertex v : mesh.vertices()){
          thisFracDWeights[v] *= (geom.vertexDualAreas[v] * measure[v]) / normD[v];//multiplying by area (density -> mass)
          //thisFracDWeights[v] *= measure[v]/normD[v];
          updateWSum += thisFracDWeights[v];
        }

        gradNorm += (desiredMass-updateWSum) * (desiredMass-updateWSum);

        // std::cout << "updateWSum = " << updateWSum << std::endl;
        //newWeight = phiWeights[j] + (((1/iIter) * shortTime) * (desiredMass - updateWSum));//do a time-decaying step size here
        phiWeights[j] += 1e-0 * (desiredMass - updateWSum);//don't do a time-decaying step size
        // newWeight = phiWeights[j] + ((shortTime) * (desiredMass - updateWSum));//don't do a time-decaying step size
        // mean += newWeight; 

        cellMasses[j] = updateWSum;
      }

      gradNorm = sqrt(gradNorm);
      // H(gradNorm);
      std::cout << "gradNorm: " << gradNorm << "\t\r" << std::flush;


      if (gradNorm < 1e-6)
        break;

      psMesh.addVertexScalarQuantity("rho", rho);
      
      //ensure that the average of all the weights is 0
      // mean /= nSites;
      // for (size_t k = 0; k < nSites; k++){
      //   newPhiWeights[k] -= mean;
      // }
      // phiWeights = newPhiWeights;
    }


    // std::cout << "weights after gradient descent: " << std::endl;
    // for (int i = 0; i < phiWeights.size(); i++){
    //   std::cout << "weight at site " << i << ": " << phiWeights[i] << std::endl;
    // }

    std::cout << "Cell masses:" << std::endl;
    for (int i = 0; i < nSites; i++)
      std::cout << cellMasses[i] << std::endl;

    // Compute the normalizer distribution
    rhs = computeRHSWithWeights(siteLocations, phiWeights, shortTime);
    normRHS = VertexData<double>(mesh, 0); // sum on all sites
    for (int i = 0; i < siteLocations.size(); i++)
      normRHS += rhs[i];
    normD = vSolver.scalarDiffuse(normRHS);


    double energy = 0;
    std::vector<SurfacePoint> newSiteLocations;

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> thisFracDKarcher = vSolver.scalarDiffuse(rhs[iSite]);

      for (Vertex v : mesh.vertices()) {
        thisFracDKarcher[v] /= normD[v];
      }

      //WEIGHT THE DISTRIBUTION BY THE CURL MEASURE (vertex dual areas below when finding weight)
      for (Vertex v : mesh.vertices()){ 
        thisFracDKarcher[v] *= measure[v];
      }

      //visualize the cells
      //result.cellEvolution[iSite][iIter] = thisFracDKarcher;

      // TODO: Swap loops
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
