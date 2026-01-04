#include "powerCells.h"

std::vector<std::pair<SurfacePoint, SurfacePoint>> computeCourseSingularities(powerCellOptions& options){

    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    auto saddleLoops = options.saddleLoops;
    double period = options.period;
    //cut edges
    std::unordered_set<Edge> cutEdges;


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
                    cutEdges.insert(e);
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
    //polyscope::show();

    //find the decompositions per face 
    FaceData<int> faceComponentsGlued =  componentsCutByLoops(*gluedGeometry, saddleLoops);

    //--------------------------------------------------------------------------//
    //Build component meshes 
    std::unordered_map<int, std::vector<Face>> facesByComponent;
    for (Face f : gluedGeometry->mesh.faces()) {
        int cid = faceComponentsGlued[f];
        facesByComponent[cid].push_back(f);
    }
    std::vector<ComponentMesh> componentMeshes;

    for (auto& [cid, faces] : facesByComponent) {
        componentMeshes.push_back(
            extractComponent(
                globalMesh,
                *globalGeometry,
                faces
            )
        );
    }
    // for (size_t i = 0; i < componentMeshes.size(); ++i) {
    //     auto& cm = componentMeshes[i];
    //     std::cout << "Does component " << i << " have boundary? " << cm.mesh->hasBoundary() << std::endl;
    //     std::cout << "Is component " << i << " manifold? " << cm.mesh->isManifold() << std::endl;
    //     std::cout << "Number of boundary loops in component " << i << " " << cm.mesh->nBoundaryLoops() << std::endl;
    //     cm.geom->requireDECOperators();
    //     // --- Register the surface mesh ---
    //     auto* psMesh = polyscope::registerSurfaceMesh(
    //         "Component " + std::to_string(i),
    //         cm.geom->vertexPositions,
    //         cm.mesh->getFaceVertexList()
    //     );

    //     // --- Collect boundary vertices ---
    //     std::vector<Vector3> boundaryPositions;
    //     boundaryPositions.reserve(cm.mesh->nVertices());

    //     for (Vertex v : cm.mesh->vertices()) {
    //         if (v.isBoundary()) {
    //             boundaryPositions.push_back(cm.geom->vertexPositions[v]);
    //         }
    //     }

    //     // --- Register boundary vertices as a point cloud ---
    //     if (!boundaryPositions.empty()) {
    //         auto* psBoundary = polyscope::registerPointCloud(
    //             "Component " + std::to_string(i) + " boundaries",
    //             boundaryPositions
    //         );
    //     }
    // }
    // polyscope::show();

    //---------------------------------------------------------------------------//

    // find k = #components
    int maxID = 0;
    for (Face f : gluedMesh.faces()) {
        maxID = std::max(maxID, faceComponentsGlued[f]);
    }
    int k = maxID + 1;
    bool cylindricalParameterization = false;
    if (k == 1){
        cylindricalParameterization = true;//first test cylindrical parameterization on cylinders only
    }
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
        auto matchedPairs = computeBucketSingularities(options, bucketCurl, allDist, cylindricalParameterization);
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
        polyscope::registerPointCloud("(+) Aligned Pts for bucket " + std::to_string(i), positiveCenters)->setEnabled(false);
        polyscope::registerPointCloud("(-) Aligned Pts for bucket " + std::to_string(i), negativeCenters)->setEnabled(false);
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

std::vector<std::pair<SurfacePoint, SurfacePoint>> computeBucketSingularities(powerCellOptions& options, VertexData<double>& curlMeasure, VertexData<double>& allDist,
                                                                            bool cylindricalParameterization){
    
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

    // pair positive and negative sites by sorted time
    std::vector<std::pair<SurfacePoint, SurfacePoint>> matchedPairs;

    if (options.usePathCost){
        //match pairs via an OT problem where the cost is ||\Delta h||^2 + ||Path Cost||^2
        matchedPairs = timeFuncAndDistOTMatching(options, posVoronoiCenters.siteLocations, negVoronoiCenters.siteLocations);
    }
    else{
        // sort by time function value ascending
        auto cmp = [](const std::pair<SurfacePoint,double>& a,
              const std::pair<SurfacePoint,double>& b) {
            return a.second < b.second;
        };

        std::sort(posSiteTimeFunctions.begin(), posSiteTimeFunctions.end(), cmp);
        std::sort(negSiteTimeFunctions.begin(), negSiteTimeFunctions.end(), cmp);

        // ensure sizes match
        int n = std::min(posSiteTimeFunctions.size(), negSiteTimeFunctions.size());

        matchedPairs.reserve(n);

        for (int i = 0; i < n; i++) {
            matchedPairs.emplace_back(
            posSiteTimeFunctions[i].first,
            negSiteTimeFunctions[i].first
            );
        }
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

    auto sameFace = [&](const SurfacePoint& spPos, const SurfacePoint& spNeg){
        Face fPos1 = spPos.edge.halfedge().face();
        Face fPos2 = spPos.edge.halfedge().twin().face();

        Face fNeg1 = spNeg.edge.halfedge().face();
        Face fNeg2 = spNeg.edge.halfedge().twin().face();

        
        if (fPos1 == fNeg1) return true;
        if (fPos1 == fNeg2) return true;
        if (fPos2 == fNeg1) return true;
        if (fPos2 == fNeg2) return true;
        return false;
    };

    //do some filtering of the matched pairs
    for (auto& pr : matchedPairs) {
        const SurfacePoint& spPos = pr.first;
        const SurfacePoint& spNeg = pr.second;

        double dPos = spPos.interpolate(allDist);
        double dNeg = spNeg.interpolate(allDist);
        // Keep only if BOTH ends are far enough from saddle set and if they're not from the same face! 
        if (dPos >= cutoff && dNeg >= cutoff && !sameFace(spPos, spNeg)) {
            filteredPairs.push_back(pr);
        }
    }

    matchedPairs = std::move(filteredPairs);
    return matchedPairs;
}

//extract components without preserving boundaries between cylindrical components
//in particular, nBoundaryLoops for each cylindrical component = 0
ComponentMesh extractComponent(SurfaceMesh& srcMesh, VertexPositionGeometry& srcGeom, const std::vector<Face>& faces){
    // --- assign new indices to used vertices ---
    std::unordered_map<Vertex, size_t> oldToIndex;
    std::vector<Vertex> indexToOld;

    for (Face f : faces) {
        for (Vertex v : f.adjacentVertices()) {
            if (oldToIndex.count(v) == 0) {
                size_t idx = indexToOld.size();
                oldToIndex[v] = idx;
                indexToOld.push_back(v);
            }
        }
    }

    // --- build polygon list ---
    std::vector<std::vector<size_t>> polygons;
    polygons.reserve(faces.size());

    for (Face f : faces) {
        std::vector<size_t> poly;
        for (Vertex v : f.adjacentVertices()) {
            poly.push_back(oldToIndex[v]);
        }
        polygons.push_back(poly);
    }

    // --- create new mesh ---
    auto newMesh = std::make_unique<SurfaceMesh>(polygons);

    // --- transfer vertex positions ---
    VertexData<Vector3> newPositions(*newMesh);

    for (size_t i = 0; i < indexToOld.size(); ++i) {
        Vertex oldV = indexToOld[i];
        Vertex newV = newMesh->vertex(i);
        newPositions[newV] = srcGeom.vertexPositions[oldV];
    }

    auto newGeom = std::make_unique<VertexPositionGeometry>(*newMesh, newPositions);

    // --- vertex build maps ---
    ComponentMesh result;
    result.mesh = std::move(newMesh);
    result.geom = std::move(newGeom);

    result.newToOldV.resize(indexToOld.size());
    for (size_t i = 0; i < indexToOld.size(); ++i) {
        Vertex newV = result.mesh->vertex(i);
        Vertex oldV = indexToOld[i];
        result.oldToNewV[oldV] = newV;
        result.newToOldV[i] = oldV;
    }

    // --- build edge maps ---
    result.newToOldE.resize(result.mesh->nEdges());

    for (Edge newE : result.mesh->edges()) {
        // endpoints in new mesh
        Vertex newA = newE.firstVertex();
        Vertex newB = newE.secondVertex();

        // corresponding old vertices
        Vertex oldA = result.newToOldV[newA.getIndex()];
        Vertex oldB = result.newToOldV[newB.getIndex()];

        // check if old mesh had this edge
        Edge oldE;
        for (Halfedge he : oldA.outgoingHalfedges()) {
            if (he.tipVertex() == oldB) {
                oldE = he.edge();
                break;
            }
        }
        if (oldE != Edge()) {
            result.oldToNewE[oldE] = newE;
            result.newToOldE[newE.getIndex()] = oldE;
        } else {
            // This can happen on boundaries or cut edges
            result.newToOldE[newE.getIndex()] = Edge();
        }
    }

    return result;
}


std::vector<std::pair<SurfacePoint, SurfacePoint>> timeFuncAndDistOTMatching(powerCellOptions& options, std::vector<SurfacePoint>& posVoronoiCenters, 
                                                                            std::vector<SurfacePoint>& negVoronoiCenters){

    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    double period = options.period;
    VertexData<double> timeFunction = options.timeFunction;
    HalfedgeData<double> heWeights = options.heWeights;

    // --------------------------------------------
    // Build site list (+1 = pos, -1 = neg)
    // --------------------------------------------
    std::vector<std::pair<SurfacePoint, int>> sites;
    std::vector<std::pair<SurfacePoint, SurfacePoint>> toReturn;

    for (const auto& sp : posVoronoiCenters)
        sites.emplace_back(sp, +1);

    for (const auto& sp : negVoronoiCenters)
        sites.emplace_back(sp, -1);

    int numSites = (int)sites.size();

    // --------------------------------------------
    // Cost matrix + raw path cost storage
    // --------------------------------------------
    Eigen::MatrixXd C(numSites, numSites);
    C.setConstant(1e7); // default high cost

    Eigen::MatrixXd rawPathCost(numSites, numSites);
    rawPathCost.setZero();

    double minPathCost = std::numeric_limits<double>::infinity();
    double maxPathCost = -std::numeric_limits<double>::infinity();

    // --------------------------------------------
    // First pass: compute raw costs
    // --------------------------------------------
    for (int i = 0; i < numSites; ++i) {
        int singValI = sites[i].second;

        for (int j = 0; j < numSites; ++j) {
            int singValJ = sites[j].second;

            // Only allow + → -
            if (singValI < 0 && singValJ > 0) continue;
            if (singValI == singValJ) continue;

            // --- time function term ---
            double timeFuncI = sites[i].first.interpolate(timeFunction);
            double timeFuncJ = sites[j].first.interpolate(timeFunction);
            double deltaH = std::abs(timeFuncI - timeFuncJ);

            // --- path cost term ---
            Vertex vStart = sites[i].first.nearestVertex();
            Vertex vEnd   = sites[j].first.nearestVertex();

            double pathCost =
                computePathCost(*globalGeometry, heWeights, vStart, vEnd);

            rawPathCost(i, j) = pathCost;

            minPathCost = std::min(minPathCost, pathCost);
            maxPathCost = std::max(maxPathCost, pathCost);
        }
    }

    // --------------------------------------------
    // Second pass: normalize + fill C
    // --------------------------------------------
    double denom = maxPathCost - minPathCost;
    if (denom < 1e-12) denom = 1.0; // safety

    for (int i = 0; i < numSites; ++i) {
        int singValI = sites[i].second;

        for (int j = 0; j < numSites; ++j) {
            int singValJ = sites[j].second;

            if (singValI < 0 && singValJ > 0) continue;
            if (singValI == singValJ) continue;

            // normalized path cost in [0,1]
            double normPathCost =
                (rawPathCost(i, j) - minPathCost) / denom;

            // time function term (already in reasonable scale)
            double timeFuncI = sites[i].first.interpolate(timeFunction);
            double timeFuncJ = sites[j].first.interpolate(timeFunction);
            double deltaH = std::abs(timeFuncI - timeFuncJ);

            // combine (weights optional)
            double alpha = 1.0; // path importance
            double beta  = 1.0; // isoline importance

            C(i, j) = alpha * normPathCost + beta * deltaH;
        }
    }

    //solve the bipartite matching model
    using namespace std;
    try {
        
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //add a matrix of Gurobi variables that signify the matching
        std::vector<std::vector<GRBVar>> T;
    
        for (int i = 0; i < numSites; i++){
            std::vector<GRBVar> currRow;
            for (int j = 0; j < numSites; j++){
                GRBVar t = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
                currRow.push_back(t);
            }
            T.push_back(currRow);
        }

        //first constraint (T\mathbb{1} = \mathbb{1})
        for (int i = 0; i < numSites; i++){
            GRBLinExpr lhs = 0;
            for (int j = 0; j < numSites; j++){
                lhs += T[i][j];
            }
            model.addConstr(lhs == 1.0);
        }

        //second constraint (\mathbb{1}^T T = \mathbb{1})
        for (int i = 0; i < numSites; i++){
            GRBLinExpr lhs = 0;
            for (int j = 0; j < numSites; j++){
                lhs += T[j][i];
            }
            model.addConstr(lhs == 1.0);
        }

        //add the objective
        GRBLinExpr matProd = 0;
        for (int i = 0; i < numSites; i++){
            for (int j = 0; j < numSites; j++){
                matProd += C(i, j) * T[i][j];
            }
        }

        GRBLinExpr obj = matProd;

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();
        
        for (int i = 0; i < numSites; i++){
            for (int j = 0; j < numSites; j++){
                if (T[i][j].get(GRB_DoubleAttr_X) > 0 && (sites[i].second > 0 && sites[j].second < 0)){
                    std::cout << "match surface point " << sites[i].first << " to surface point " << sites[j].first << std::endl;
                    toReturn.push_back(std::make_pair(sites[i].first, sites[j].first));//always return pairs in (posSite, negSite) order
                }
            }
        }

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return toReturn;

}