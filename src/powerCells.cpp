#include "powerCells.h"

std::vector<std::pair<SurfacePoint, SurfacePoint>> computeCourseSingularities(powerCellOptions& options){

    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    auto saddleLoops = options.saddleLoops;
    double period = options.period;


    // Create the Heat Method solver
    HeatMethodDistanceSolver heatSolver(*gluedGeometry);
    //first compute the measure in the course direction 
    VertexData<double> courseMeasure = computeCourseCurlMeasure(options); 
    //all distance from the heat vertices
    VertexData<double> allDist(gluedMesh, DBL_MAX);//need to handle this in a better way

    //if we want to mask the saddle 
    if (options.maskSaddle && saddleLoops.size() != 0){
        std::vector<Vertex> heatSourceVerts; 
        for (int i = 0; i < saddleLoops.size(); i++){
            std::vector<double> path = saddleLoops[i];
            for (int j = 0; j < path.size(); j++){
                if (std::fabs(path[j]) > 0){
                    Edge e = gluedGeometry->mesh.edge(j);
                    heatSourceVerts.push_back(e.halfedge().tailVertex());
                    heatSourceVerts.push_back(e.halfedge().tipVertex());
                }
            }
        }
        //also run the diffusion from boundary vertices
        for (Vertex v : gluedMesh.vertices()){
            if (v.isBoundary()) heatSourceVerts.push_back(v);
        }
        options.heatSourceVerts = heatSourceVerts;
        //mask the saddle 
        allDist = heatSolver.computeDistance(heatSourceVerts);
        VertexData<double> courseWeighting(globalMesh);
        double maxVal = std::numeric_limits<double>::min();
        double maxSourceVal = std::numeric_limits<double>::min();
        //for all the source vertices, find the max value 
        for (Vertex v : heatSourceVerts){
            maxSourceVal = std::max(maxSourceVal, allDist[v]);
        }
        //find the max val over all distances
        for (Vertex v : gluedMesh.vertices()){
            maxVal = std::max(maxVal, allDist[v]);
        }
        //shift down all the source vertex values
        for (Vertex v : heatSourceVerts){
            allDist[v] -= maxSourceVal;
        }
        //clip all values to 0 
        for (Vertex v : gluedMesh.vertices()){
            allDist[v] = std::max(allDist[v], 0.0);
        }
        //have a hard mask in the course direction
        for (Vertex v : globalMesh.vertices()){
            courseWeighting[v] = (allDist[v] > 1.5 * period);
            courseMeasure[v] = courseWeighting[v] * courseMeasure[v];
        }
    }//GENERAL TO-DO: \sigma should integrate to 0 around saddle loops

    psMesh -> addVertexScalarQuantity("Curl measure from PC file", courseMeasure);
    polyscope::show();

    //find the decompositions per face 
    FaceData<int> faceComponentsGlued =  componentsCutByLoops(*gluedGeometry, saddleLoops);

    // find k = #components
    int maxID = 0;
    for (Face f : gluedMesh.faces()) {
        maxID = std::max(maxID, faceComponentsGlued[f]);
    }
    int k = maxID + 1;
    //prepare buckets and assign vertex buckets
    std::vector<std::vector<Vertex>> buckets(k);
    //assign vertices (multi-membership)
    for (Vertex v : gluedMesh.vertices()) {
        std::vector<char> seen(k, 0);
        for (Face f : v.adjacentFaces()) {
            int cid = faceComponentsGlued[f];
            if (!seen[cid]) {
                buckets[cid].push_back(v);
                seen[cid] = 1;
            }
        }
    }

    std::vector<std::pair<SurfacePoint, SurfacePoint>> allMatchedPairs;


    //solve the problem per bucket of vertices
    for (int i = 0; i < buckets.size(); i++){
        VertexData<double> bucketCurl(gluedMesh, 0.);
        auto bucket = buckets[i];
        //update the measure for vertices in this bucket
        for (Vertex v : bucket){
            bucketCurl[v] = courseMeasure[v];
        }
        auto matchedPairs = computeBucketSingularities(options, bucketCurl, allDist);
        std::vector<Vector3> positiveCenters;
        for (int i = 0; i < matchedPairs.size(); i++){
            auto sp = matchedPairs[i].first;
            positiveCenters.push_back(sp.interpolate(globalGeometry->vertexPositions));
        }
        std::vector<Vector3> negativeCenters;
        for (int i = 0; i < matchedPairs.size(); i++){
            auto sp = matchedPairs[i].second;
            negativeCenters.push_back(sp.interpolate(globalGeometry->vertexPositions));
        }
        allMatchedPairs.insert(allMatchedPairs.end(), matchedPairs.begin(), matchedPairs.end());
        psMesh->addVertexScalarQuantity("Curl for bucket " + std::to_string(i), bucketCurl);
        polyscope::registerPointCloud("(+) Aligned Pts for bucket " + std::to_string(i), positiveCenters);
        polyscope::registerPointCloud("(-) Aligned Pts for bucket " + std::to_string(i), negativeCenters);

        polyscope::show();
    }

    return allMatchedPairs;
}

//compute the course curl measure in the glued setting 
VertexData<double> computeCourseCurlMeasure(powerCellOptions& options){

    SurfaceMesh& globalMesh = options.globalGeometry->mesh;
    SurfaceMesh& gluedMesh = options.gluedGeometry->mesh;
    VertexData<double> curlGlobal(globalMesh);
    VertexData<double> curlGlued(gluedMesh);
    options.globalGeometry->requireFaceAreas();
    auto vertexMap = options.vertexMap;
    auto field = options.normalizedTFGrad;
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        double area = 0.0; // area of the 1-ring of faces
        for (Halfedge he : options.gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = options.globalGeometry->vertexPositions[hjk.tipVertex()] - options.globalGeometry->vertexPositions[hjk.tailVertex()];
            field[he.face()] = field[he.face()].normalize(); //always normalize the field
            area += options.globalGeometry->faceArea(he.face());
            sum += dot(hjkVec, field[he.face()]);
        }

        if (vi.isBoundary()) curlGlobal[vi] = 0.;//set curl to 0 on boundary vertices
        else curlGlobal[vi] = sum / area;
    }

    //convert curl to glued setting 
    for (Vertex v : globalMesh.vertices()){
        curlGlued[vertexMap[v.getIndex()]] = curlGlobal[v];
    }
   
    return curlGlued;
}

std::vector<std::pair<SurfacePoint, SurfacePoint>> computeBucketSingularities(powerCellOptions& options, VertexData<double>& curlMeasure, VertexData<double>& allDist){
    
    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    double period = options.period;
    VertexData<double> timeFunction = options.timeFunction;
    
    // Cap curl measure to avoid high concentration of singularities
    for (Vertex v : gluedMesh.vertices()) {
        curlMeasure[v] = fmin(curlMeasure[v], period / (3*globalGeometry->vertexDualAreas[v]));
        curlMeasure[v] = fmax(curlMeasure[v], -period / (3*globalGeometry->vertexDualAreas[v]));
    }
                          
    //now divide the measure into positive and negative 
    VertexData<double> posMeasure(gluedMesh, 0.0);
    VertexData<double> negMeasure(gluedMesh, 0.0);
    double totalPosMeasure = 0;
    double totalNegMeasure = 0;
    for (Vertex v : gluedMesh.vertices()){
        if (curlMeasure[v] > 0){
            posMeasure[v] = curlMeasure[v];
            totalPosMeasure += gluedGeometry->vertexDualAreas[v] * posMeasure[v];
        }
        else{
            negMeasure[v] = std::fabs(curlMeasure[v]);
            totalNegMeasure += gluedGeometry->vertexDualAreas[v] * negMeasure[v];
        }
    }
    double avgTotalMeasure = (totalPosMeasure + totalNegMeasure) / 2;
    
    //positive sites
    VoronoiOptions posOptions = defaultVoronoiOptions;
    posOptions.nSites = std::round(avgTotalMeasure / period);
    posOptions.useDelaunay = false;
    posOptions.computeDistributions = true;
    posOptions.seed = 42;
    VoronoiResult posVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(globalMesh, *globalGeometry, posOptions, posMeasure, *psMesh);
    std::vector<VertexData<double>> posSiteDistributions = posVoronoiCenters.siteDistributions;
    for (size_t i = 0; i < posSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){ 
            mass += posSiteDistributions[i][v] * posMeasure[v] * globalGeometry->vertexDualAreas[v];
        }
        std::cout << "positive mass at site " << i << " = " << mass << std::endl;
    }
    

    //negative sites 
    VoronoiOptions negOptions = defaultVoronoiOptions;
    negOptions.nSites = posOptions.nSites;
    negOptions.useDelaunay = false;
    negOptions.computeDistributions = true;
    negOptions.seed = 42;
    VoronoiResult negVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(globalMesh, *globalGeometry, negOptions, negMeasure, *psMesh);
    std::vector<VertexData<double>> negSiteDistributions = negVoronoiCenters.siteDistributions;
    for (size_t i = 0; i < negSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){ 
            mass += negSiteDistributions[i][v] * negMeasure[v] * globalGeometry->vertexDualAreas[v];
        }
        std::cout << "negative mass at site " << i << " = " << mass << std::endl;
    }
    
    //render surface points and pair up time functions
    std::vector<Vector3> positiveCenters;
    std::vector<std::pair<SurfacePoint, double>> posSiteTimeFunctions;
    std::vector<std::pair<SurfacePoint, double>> negSiteTimeFunctions;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
       auto sp = posVoronoiCenters.siteLocations[i];
       positiveCenters.push_back(sp.interpolate(globalGeometry->vertexPositions));
       posSiteTimeFunctions.emplace_back(std::make_pair(sp, sp.interpolate(timeFunction)));
    }
    std::vector<Vector3> negativeCenters;
    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
       auto sp = negVoronoiCenters.siteLocations[i];
       negativeCenters.push_back(sp.interpolate(globalGeometry->vertexPositions));
       negSiteTimeFunctions.emplace_back(std::make_pair(sp, sp.interpolate(timeFunction)));
    }

    // sort by time function value ascending
    auto cmp = [](const std::pair<SurfacePoint,double>& a,
              const std::pair<SurfacePoint,double>& b) {
        return a.second < b.second;
    };

    std::sort(posSiteTimeFunctions.begin(), posSiteTimeFunctions.end(), cmp);
    std::sort(negSiteTimeFunctions.begin(), negSiteTimeFunctions.end(), cmp);

    // ensure sizes match
    int n = std::min(posSiteTimeFunctions.size(), negSiteTimeFunctions.size());

    // pair positive and negative sites by sorted time
    std::vector<std::pair<SurfacePoint, SurfacePoint>> matchedPairs;
    matchedPairs.reserve(n);

    for (int i = 0; i < n; i++) {
        matchedPairs.emplace_back(
        posSiteTimeFunctions[i].first,
        negSiteTimeFunctions[i].first
        );
    }

    //specify options for alignment
    alignOptions alignmentOptions;
    alignmentOptions.timeFunction = timeFunction;
    std::vector<std::pair<SurfacePoint, SurfacePoint>> matchedSurfacePoints;
    
    //align the points on the same isoline
    alignmentOptions.pairedSites = matchedPairs;
    alignPointsOnIsolineFast(globalMesh, *globalGeometry, alignmentOptions, *psMesh);
    matchedPairs = alignmentOptions.pairedSites;

    // Remove pairs whose endpoints are within cutoff distance of any saddle heat source
    double cutoff = 1.5 * period;

    std::vector<std::pair<SurfacePoint, SurfacePoint>> filteredPairs;
    filteredPairs.reserve(matchedPairs.size());

    for (auto& pr : matchedPairs) {
        const SurfacePoint& spPos = pr.first;
        const SurfacePoint& spNeg = pr.second;

        double dPos = spPos.interpolate(allDist);
        double dNeg = spNeg.interpolate(allDist);

        // Keep only if BOTH ends are far enough from saddle set
        if (dPos >= cutoff && dNeg >= cutoff) {
            filteredPairs.push_back(pr);
        }
    }

    matchedPairs = std::move(filteredPairs);
    return matchedPairs;
}