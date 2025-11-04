#include "powerCells.h"
void computeCourseSingularities(powerCellOptions& options){

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
        //mask the saddle 
        VertexData<double> allDist = heatSolver.computeDistance(heatSourceVerts);
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

    //solve the proble per bucket of vertices
    for (auto bucket : buckets){
        VertexData<double> bucketCurl(gluedMesh, 0.);
        //update the measure for vertices in this bucket
        for (Vertex v : bucket){
            bucketCurl[v] = courseMeasure[v];
        }
        computeBucketSingularities(options, bucketCurl);
    }
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

void computeBucketSingularities(powerCellOptions& options, VertexData<double>& curlMeasure){
    
    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    double period = options.period;
    std::cout << "period = " << period << std::endl;
    psMesh->addVertexScalarQuantity("Curl measure from PC file 0 ", curlMeasure);

    // Cap curl measure to avoid high concentration of singularities
    for (Vertex v : gluedMesh.vertices()) {
        curlMeasure[v] = fmin(curlMeasure[v], period / (3*globalGeometry->vertexDualAreas[v]));
        curlMeasure[v] = fmax(curlMeasure[v], -period / (3*globalGeometry->vertexDualAreas[v]));
    }

    psMesh->addVertexScalarQuantity("Curl measure from PC file 1 ", curlMeasure);
                              
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
    psMesh->addVertexScalarQuantity("Pos measure from PC file 2 ", posMeasure);
    psMesh->addVertexScalarQuantity("Neg measure from PC file 3 ", negMeasure);

    //positive sites
    VoronoiOptions posOptions = defaultVoronoiOptions;
    posOptions.nSites = std::round(avgTotalMeasure / period);
    posOptions.useDelaunay = false;
    posOptions.computeDistributions = true;
    posOptions.seed = 42;
    VoronoiResult posVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(globalMesh, *globalGeometry, posOptions, posMeasure, *psMesh);
    std::vector<std::vector<VertexData<double>>> posStepSiteDistribution = posVoronoiCenters.stepSiteDistribution;
    std::vector<std::vector<SurfacePoint>> posSteps = posVoronoiCenters.steps;
    std::vector<SurfacePoint> posInitialSites = posVoronoiCenters.initialSites;

    //negative sites 
    VoronoiOptions negOptions = defaultVoronoiOptions;
    negOptions.nSites = posOptions.nSites;
    negOptions.useDelaunay = false;
    negOptions.computeDistributions = true;
    negOptions.seed = 42;
    VoronoiResult negVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(globalMesh, *globalGeometry, negOptions, negMeasure, *psMesh);
    std::vector<std::vector<VertexData<double>>> negStepSiteDistribution = negVoronoiCenters.stepSiteDistribution;
    std::vector<std::vector<SurfacePoint>> negSteps = negVoronoiCenters.steps;
    std::vector<SurfacePoint> negInitialSites = negVoronoiCenters.initialSites;

    std::vector<Vector3> positiveCenters;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
       positiveCenters.push_back(posVoronoiCenters.siteLocations[i].interpolate(globalGeometry->vertexPositions));
    }
    polyscope::registerPointCloud("Voronoi sites unaligned from PC file(+)", positiveCenters)->setEnabled(true);
    
    std::vector<Vector3> negativeCenters;
    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
       negativeCenters.push_back(negVoronoiCenters.siteLocations[i].interpolate(globalGeometry->vertexPositions));
    }
    polyscope::registerPointCloud("Voronoi sites unaligned from PC file (-)", negativeCenters)->setEnabled(false);
    //polyscope::show();
}