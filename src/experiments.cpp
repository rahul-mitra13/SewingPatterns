#include "experiments.h"

//--------------------------Implementation of strategy 1---------------------------------------//
std::tuple<HalfedgeData<double>, EdgeData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap, 
                                                                    polyscope::SurfaceMesh& psMesh, globalBoundaryConditions& boundaryConditions, double period){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    HalfedgeData<double> oldGluedSigmaTilde(gluedMesh);
    HalfedgeData<double> newGluedSigmaTilde(gluedMesh);
    HalfedgeData<double> globalSigmaTilde(globalMesh);
    //tracking objective values
    double currObj, oldObj;
    //vertex singularities for polyscope viz 
    VertexData<double> vertexSingularities(globalMesh, 0.0);
    //edge singularities for polyscope viz 
    EdgeData<double> edgeSingularities(globalMesh, 0.0);
    //curl per vertex 
    VertexData<double> vertexCurl(globalMesh);
    //curl per edge
    EdgeData<double> edgeCurl(globalMesh);
    //singular vertex indices 
    std::vector<std::pair<int, int>> singularVertices;
    //singular edges for gurobi optimization
    std::vector<std::pair<int, int>> singularEdges;
    //iso values we've placed singularities
    std::vector<double> usedIsoVals;
    
    //striping information 
    //global data
    CornerData<double> stripeValuesSigmaCourse;
    //global data
    FaceData<int> stripeIndicesSigmaCourse;
    std::vector<Vector3> positionsCourse;
    std::vector<std::array<int, 2>> edgesCourse;

    int maxPairs = 4;
    int numPairs = 0;

    double eps = 1e-3;

    //compute curl per vertex of gradient field (in the global setting)
    vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, globalFaceGradients, gluedOneRingMap);
    //compute curl per edge of the gradient field (in the global setting)
    edgeCurl = computeEdgeCurl(globalGeometry, gluedGeometry, globalFaceGradients, globalToGluedEdgeMap);

    //visualize the vertex curl 
    //this quantity never changes through the optimization 
    psMesh.addVertexScalarQuantity("vertex curl", vertexCurl);
    //visualize the edge curl 
    psMesh.addEdgeScalarQuantity("edge curl after placing " + std::to_string(numPairs) + " singularity pairs", edgeCurl);

    //find max/min curl vertex over entire mesh
    //std::pair<int, int> singVertexPair = findVertexSingularityPair(globalGeometry, gluedGeometry, vertexCurl,
    //                                    gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, 0.0, numPairs, true);
    
    //find max/min curl edge over entire mesh 
    std::pair<int, int> singEdgePair = findEdgeSingularityPair(globalGeometry, gluedGeometry, edgeCurl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, 0.0, numPairs, true);
    
    //vertexSingularities[globalMesh.vertex(singVertexPair.first)] = 1.0;
    //vertexSingularities[globalMesh.vertex(singVertexPair.second)] = -1.0;
    edgeSingularities[globalMesh.edge(singEdgePair.first)] = 1.0;
    edgeSingularities[globalMesh.edge(singEdgePair.second)] = -1.0;
    //singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.first], 1.0));
    //singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.second], -1.0));
    

    //make the model
    //we will solve it in the glued mesh setting 
    Model model;
    std::vector<std::array<double, 3>> faceGradients(gluedMesh.nFaces());
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints; 
    for (Face f : globalMesh.faces()){
        faceGradients[f.getIndex()] = std::array{globalFaceGradients[f][0], globalFaceGradients[f][1], globalFaceGradients[f][2]};
    }
    model.setFaceGradients(faceGradients);
    model.setPeriod(period);
    model.setBdyEdges(boundaryConditions.courseBdyEdges);

    //add bdy-bdy path constraints for non-collapse
    for (int i = 0; i < boundaryConditions.bdyBdyPathConstraints.size(); i++){
        //visualizing bdy-bdy edge constraints
        // EdgeData<double> bdyBdyPath(globalMesh);
        // for (Edge e : globalMesh.edges()){
        //     bdyBdyPath[e] = boundaryConditions.bdyBdyPathConstraints[i][globalToGluedEdgeMap[e.getIndex()]];
        // }
        // psMesh.addEdgeScalarQuantity("bdy bdy path " + std::to_string(i), bdyBdyPath);
        edgePathConstraints.push_back(std::make_pair(boundaryConditions.bdyBdyPathConstraints[i], 1.0));
    }
    //set non-collapse constraint
    model.setEdgePathConstraints(edgePathConstraints);

    //solve the model without any singularities 
    std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;

    //find first pair of singular edges  (vertex-curl approach)
    // std::pair<int, int> singularEdgeIndices = findSingularEdges(globalGeometry, gluedGeometry, singVertexPair, gluedOneRingMap, model);
    // singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.first], 1));
    // singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.second], -1));

    //first pair of singular edges (edge-curl approach)
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.first], -1));
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.second], 1));
    model.setSingularEdges(singularEdges);
    numPairs++;
    //psMesh.addVertexScalarQuantity("vertex singularity after inserting " + std::to_string(numPairs) + " singularity pairs", vertexSingularities);
    psMesh.addEdgeScalarQuantity("edge singularities after inserting " + std::to_string(numPairs) + " singularity pairs", edgeSingularities);
    //solve the model with 1 pair of singularities
    std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);

    //FINDING SINGULAR EDGES DIRECTLY USING EDGE CURL (HAVE TO WRITE THIS UP)
    while(true){
        //gradient of sigmaTilde in the global setting 
        FaceData<Vector3> gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, oldGluedSigmaTilde);
        psMesh.addFaceVectorQuantity("grad sigmaTilde after placing " + std::to_string(numPairs) + " singularity pairs", gradSigmaTilde);
        //use the previous iterates in place of the gradient of the time function when re-optimizing
        // std::vector<std::array<double, 3>> newGrads(gluedMesh.nFaces());
        // for (Face f : globalMesh.faces()){
        //     newGrads[f.getIndex()] = std::array{gradSigmaTilde[globalMesh.face(f.getIndex())][0], gradSigmaTilde[globalMesh.face(f.getIndex())][1], 
        //                             gradSigmaTilde[globalMesh.face(f.getIndex())][2]};
        // }
        // model.setFaceGradients(newGrads);
        //recompute edge curl using gradSigmaTilde
        edgeCurl = computeEdgeCurl(globalGeometry, gluedGeometry, gradSigmaTilde, globalToGluedEdgeMap);
        psMesh.addEdgeScalarQuantity("edge curl after placing " + std::to_string(numPairs) + " singularity pairs", edgeCurl);
        double isoVal = findIsoValWithMaxAvgEdgeCurl(globalGeometry, gluedGeometry, edgeCurl, gluedTimeFunction,
                                                    usedIsoVals, globalToGluedVertexMap);
        if (std::fabs(isoVal - (-1.0) < 1e-15)){
            std::cout << "Breaking cause we can't find any more sensible isovalues " << std::endl;
            break;//our function couldn't find any more isovalues
        } 
        std::cout << "next isoval = " << isoVal << std::endl;
        singEdgePair = findEdgeSingularityPair(globalGeometry, gluedGeometry, edgeCurl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, isoVal, numPairs, false);
        if (std::find(singularEdges.begin(), singularEdges.end(), std::make_pair(globalToGluedEdgeMap[singEdgePair.first], -1.0)) != singularEdges.end()
            || std::find(singularEdges.begin(), singularEdges.end(), std::make_pair(globalToGluedEdgeMap[singEdgePair.second], 1.0)) != singularEdges.end()){
            //don't use edges you've used before
            usedIsoVals.push_back(isoVal);
            continue;
        }
        usedIsoVals.push_back(isoVal);
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.first], -1));
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.second], 1));
        model.setSingularEdges(singularEdges);
        //solve the optimization problem 
        //sigmaTilde is in the glued mesh setting 
        std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
        if (currObj > oldObj){//current objective become worse or staying the same 
            std::cout << "Breaking cause objective no longer improving..." << std::endl;
            std::cout << "oldObj " << oldObj << std::endl;
            std::cout << "currObj " << currObj << std::endl;
            break;
        }
        //increment, update and visualize
        edgeSingularities[globalMesh.edge(singEdgePair.first)] = 1.0;
        edgeSingularities[globalMesh.edge(singEdgePair.second)] = -1.0;
        numPairs++;
        psMesh.addEdgeScalarQuantity("edge singularities after inserting " + std::to_string(numPairs) + " singularity pairs", edgeSingularities);
        oldObj = currObj;
        oldGluedSigmaTilde = newGluedSigmaTilde;
        std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
        courseStripes -> setRadius(0.001);
        courseStripes -> setEnabled(false);
    }
    
    return std::tie(oldGluedSigmaTilde, edgeSingularities);

}

//compute the per-face gradient of a 1-form (in global setting)
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    FaceData<Vector3> gradients(globalMesh);

    globalGeometry.requireFaceNormals();
    for (Face f : globalMesh.faces()){
        Eigen::MatrixXd faceSystem(3, 3);
        Eigen::VectorXd faceSigmas(3);
        Vector3 hijVec = globalGeometry.vertexPositions[f.halfedge().tipVertex()] - globalGeometry.vertexPositions[f.halfedge().tailVertex()];
        Vector3 hjkVec = globalGeometry.vertexPositions[f.halfedge().next().tipVertex()] - 
                            globalGeometry.vertexPositions[f.halfedge().next().tailVertex()];
        faceSigmas(0) = sigmaTilde[gluedMesh.face(f.getIndex()).halfedge()];
        faceSigmas(1) = sigmaTilde[gluedMesh.face(f.getIndex()).halfedge().next()];
        faceSigmas(2) = 0.0;

        faceSystem << hijVec[0], hijVec[1], hijVec[2], 
                        hjkVec[0], hjkVec[1], hjkVec[2], 
                        globalGeometry.faceNormals[f][0], globalGeometry.faceNormals[f][1], globalGeometry.faceNormals[f][2];
        
        Eigen::VectorXd soln = faceSystem.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(faceSigmas);
        gradients[f] = Vector3{soln(0), soln(1), soln(2)};
        //normalize the computed gradient 
        gradients[f] = gradients[f].normalize();
    }
    return gradients;
}

//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes (in global setting)
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> curl(globalMesh);
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        for (Halfedge he : gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()];
            sum += dot(hjkVec, field[he.face()]);
        }
        curl[vi] = sum;
    }


    return curl;
}

//compute curl per edge in the global setting 
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                FaceData<Vector3>& globalFaceGradients, std::map<int, int>& globalToGluedEdgeMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    EdgeData<double> globalCurl(globalMesh);
    EdgeData<double> gluedCurl(gluedMesh, 0.0);

    HalfedgeData<double> globalOmegaTilde(globalMesh);
    HalfedgeData<double> gluedOmegaTilde(gluedMesh);

    gluedGeometry.requireEdgeLengths();

    for (Face f : globalMesh.faces()){
        Vector3 grad = globalFaceGradients[f];
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = globalGeometry.vertexPositions[he.tipVertex()] - globalGeometry.vertexPositions[he.tailVertex()];
            //omega tilde in the global setting
            globalOmegaTilde[he] = dot(heVec, grad);
        }
    }

    for (Face f : gluedMesh.faces()){
        gluedOmegaTilde[f.halfedge()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge()];
        gluedOmegaTilde[f.halfedge().next()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge().next()];
        gluedOmegaTilde[f.halfedge().next().next()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge().next().next()];
    }

    for (Edge e : gluedMesh.edges()){
        if (e.isBoundary()) continue;
        gluedCurl[e] = (gluedOmegaTilde[e.halfedge()] + gluedOmegaTilde[e.halfedge().twin()]) / gluedGeometry.edgeLengths[e];
    }

    globalCurl = convertGluedToGlobalEdgeFunction(globalGeometry, gluedGeometry, gluedCurl, globalToGluedEdgeMap);

    return globalCurl;

}

//find max/min curl vertex for a given isoline (on global setting)
std::pair<int, int> findVertexSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllVertices){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);
    
    if (useAllVertices){
        int maxVertex, minVertex = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        //find max curl over entire mesh 
        //find vertex with max  curl
        for (Vertex v : globalMesh.vertices()){
            if (curl[v] > maxCurl){
                maxCurl = curl[v];
                maxVertex = v.getIndex();
            }
        }
        double isoVal = globalTimeFunction[globalMesh.vertex(maxVertex)];
        usedIsoVals.push_back(isoVal);
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        //polyscope::registerCurveNetwork("vertex isoline: " + std::to_string(isoVal) + " for pair " + std::to_string(numPairs) , iV, iE);
        //find min curl on the same isoline as the max vertex 
        //find vertex with min curl
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            for (Vertex v : currFace.adjacentVertices()){
                if (curl[v] < minCurl){
                    minCurl = curl[v];
                    minVertex = v.getIndex();
                }
            }
        }
        return std::make_pair(maxVertex, minVertex);
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        int maxVertex, minVertex = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        //polyscope::registerCurveNetwork("vertex isoline: " + std::to_string(isoVal) + " for pair " + std::to_string(numPairs) , iV, iE);
        for (int i = 0; i < f.size(); i++){
            //if (std::find(usedVertices.begin(), usedVertices.end(), v.getIndex()) != usedVertices.end()) continue;//don't use vertices we've already used 
            Face currFace = globalMesh.face(f[i]);
            for (Vertex v : currFace.adjacentVertices()){
                if (curl[v] > maxCurl){
                    maxCurl = curl[v];
                    maxVertex = v.getIndex();
                }
                if (curl[v] < minCurl){
                    minCurl = curl[v];
                    minVertex = v.getIndex();
                }
            }
        }
        return std::make_pair(maxVertex, minVertex);
    }
}

//find max/min curl edge for a given isoline (on global setting)
std::pair<int, int> findEdgeSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllEdges){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);


    if (useAllEdges){
        int maxEdge, minEdge = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        //find max curl over entire mesh 
        //find edge with max  curl
        for (Edge e : globalMesh.edges()){
            if (curl[e] > maxCurl){
                maxCurl = curl[e];
                maxEdge = e.getIndex();
            }
        }
        //average isovalue of the two vertices on that edge
        double isoVal = 0.5 * (globalTimeFunction[globalMesh.edge(maxEdge).halfedge().tailVertex()] +
                                     globalTimeFunction[globalMesh.edge(maxEdge).halfedge().tipVertex()]);
        
        usedIsoVals.push_back(isoVal);
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        isoline -> setRadius(0.001);
        isoline -> setEnabled(false);
        //find min curl on the same isoline as the max edge 
        //find edge with min curl
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            for (Edge e : currFace.adjacentEdges()){
                //check if current isoline crosses this edge
                if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                 || (isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                    if (curl[e] < minCurl){
                        minCurl = curl[e];
                        minEdge = e.getIndex();
                    }
                }
            }
        }
        return std::make_pair(maxEdge, minEdge);
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        isoline -> setEnabled(false);
        isoline -> setRadius(0.001);
        int maxEdge, minEdge = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            for (Edge e : currFace.adjacentEdges()){
                //check if current isoline crosses this edge
                if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    ||(isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                    if (curl[e] > maxCurl){
                        maxCurl = curl[e];
                        maxEdge = e.getIndex();
                    }
                }
                //check if current isoline crosses this edge
                if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    ||(isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                    if (curl[e] < minCurl){
                        minCurl = curl[e];
                        minEdge = e.getIndex();
                    }
                }
            }
        }
        return std::make_pair(maxEdge, minEdge);
    }
}

//solve the optimization problem for strategy 1 (in the glued setting)
std::tuple<HalfedgeData<double>, double> computeStrategy1_oneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);
    double objectiveVal;
    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Solving strategy 1 opimtization " << std::endl;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::vector<double>> waleBdyPathConstraints = gbModel.getWaleBdyPathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<std::array<double, 3>> gradients = gbModel.getFaceGradients();

    //build the gradient operator 
    Eigen::MatrixXd V(gluedMesh.nVertices(), 3);
    Eigen::MatrixXi F(gluedMesh.nFaces(), 3);
    std::tie(V, F) = getVertexPositionsandFaceLists(globalGeometry);
    // Compute the global gradient operator: #F*3 by #V
    Eigen::SparseMatrix<double> grad;
    igl::grad(V,F,grad);
    Eigen::SparseMatrix<double, Eigen::RowMajor> G(grad);
    //require the face areas
    gluedGeometry.requireFaceAreas();

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        //defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //add boundary integral value for wale direction stripes
        std::vector<GRBVar> waleBdyIntegerConstraints;
        for (size_t i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            waleBdyIntegerConstraints.push_back(k_i);
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        //compute nP over every face 
        //add the second constraint while we're here
        //second constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            nP[f.getIndex()] = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0.0);
        }

        //third constraint 
        //ensure that values are opposite sign across halfedges
        //if it's not a singular halfedge 
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            //invert the sign of the constraint
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] == -1.0 * p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //fourth constraint 
        //add constraints in the wale direction
        for (int i = 0; i < waleBdyPathConstraints.size(); i++){
            std::vector<double> path = waleBdyPathConstraints[i];
            GRBLinExpr pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(j).halfedge().getIndex()] = std::fabs(path[j]);
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = std::fabs(path[j]);
                }
            }
            for (int k = 0; k < gluedMesh.nHalfedges(); k++){
                pathIntegral += hePath[k] * sigma[k];
            }
            model.addConstr(pathIntegral == period * waleBdyIntegerConstraints[i]);
        }

        //compute a piecewise linear function over the vertices of the mesh 
        std::vector<GRBLinExpr> u(gluedMesh.nVertices());
        std::vector<std::vector<GRBLinExpr>> gradU(gluedMesh.nFaces(), std::vector<GRBLinExpr>(3));
        for (Face f : gluedMesh.faces()){

            u[f.halfedge().vertex().getIndex()] = 0.0;
            u[f.halfedge().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()]) - (nP[f.getIndex()]/3.0);
            u[f.halfedge().next().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()] + sigma[f.halfedge().next().getIndex()]) - ((2.0 * nP[f.getIndex()])/3.0);

            GRBLinExpr currGradU = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][2] = currGradU;
        }

        //set up the difference for L2 norms 
        //GRBQuadExpr gradDiff_L2 = 0;
        //set up the difference for L1 norms
        GRBLinExpr gradDiff_L1 = 0;
        std::vector<std::array<double, 3>> grads;
        //store the per face difference 
        std::vector<GRBQuadExpr> difference(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            //objective is 2-norm squared 
            // GRBQuadExpr diffX = (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]) * (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]);
            // GRBQuadExpr diffY = (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]) * (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]);
            // GRBQuadExpr diffZ = (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]) * (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]);

            //objective is 1-norm
            GRBLinExpr diffX = (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]);
            GRBVar absDiffX = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "abs_DiffX");
            GRBLinExpr diffY = (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]);
            GRBVar absDiffY = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "abs_DiffY");
            GRBLinExpr diffZ = (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]);
            GRBVar absDiffZ = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "abs_DiffZ");
            model.addConstr(absDiffX >= diffX);
            model.addConstr(absDiffX >= -diffX);
            model.addConstr(absDiffY >= diffY);
            model.addConstr(absDiffY >= -diffY);
            model.addConstr(absDiffZ >= diffZ);
            model.addConstr(absDiffZ >= -diffZ);

            double currArea = gluedGeometry.faceAreas[f];
            gradDiff_L1 += currArea * (absDiffX + absDiffY + absDiffZ);
            difference[f.getIndex()] = currArea * (diffX + diffY + diffZ);
        }

        GRBQuadExpr obj = gradDiff_L1;

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        std::cout << "Objective after placing " << singularEdges.size()/2 << " singularity pairs: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }     
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    std::cout << "--------------------------------" << std::endl;
    return std::tie(gluedOneForm, objectiveVal);
}

//trying to solve for a 1-form that only optimizes for equally spaced stripes 
//Matteo's idea: The objective term is || ||\delta \sigma||^2 - 1||^2 
//Have a bdy-bdy path integral for non-collapse behavior
std::tuple<HalfedgeData<double>, double> computeEquallySpacedOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                                std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);
    double objectiveVal;
    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Solving equally spaced stripes 1-form " << std::endl;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::vector<double>> waleBdyPathConstraints = gbModel.getWaleBdyPathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<std::array<double, 3>> gradients = gbModel.getFaceGradients();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();

    //build the gradient operator 
    Eigen::MatrixXd V(gluedMesh.nVertices(), 3);
    Eigen::MatrixXi F(gluedMesh.nFaces(), 3);
    std::tie(V, F) = getVertexPositionsandFaceLists(globalGeometry);
    // Compute the global gradient operator: #F*3 by #V
    Eigen::SparseMatrix<double> grad;
    igl::grad(V,F,grad);
    Eigen::SparseMatrix<double, Eigen::RowMajor> G(grad);
    //require the face areas
    gluedGeometry.requireFaceAreas();

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        //defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //add boundary integral value for wale direction stripes
        std::vector<GRBVar> waleBdyIntegerConstraints;
        for (size_t i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            waleBdyIntegerConstraints.push_back(k_i);
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        //compute nP over every face 
        //add the second constraint while we're here
        //second constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        //std::vector<GRBVar> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            nP[f.getIndex()] = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0.0);
        }

        //third constraint 
        //ensure that values are opposite sign across halfedges
        //if it's not a singular halfedge 
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            //invert the sign of the constraint
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] == -1.0 * p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //fourth constraint 
        //add constraints in the wale direction
        for (int i = 0; i < waleBdyPathConstraints.size(); i++){
            std::vector<double> path = waleBdyPathConstraints[i];
            GRBLinExpr pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(j).halfedge().getIndex()] = std::fabs(path[j]);
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = std::fabs(path[j]);
                }
            }
            for (int k = 0; k < gluedMesh.nHalfedges(); k++){
                pathIntegral += hePath[k] * sigma[k];
            }
            model.addConstr(pathIntegral == period * waleBdyIntegerConstraints[i]);
        }

        //fifth constraint: add bdy-bdy path constraint 
        for (int i = 0; i < edgePathConstraints.size(); i++){
            std::vector<double> path = edgePathConstraints[i].first;
            GRBLinExpr pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(i).halfedge().getIndex()] = path[j];
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(i).halfedge().twin().getIndex()] = path[j];
                }
            }
            for (int k = 0; k < gluedMesh.nHalfedges(); k++){
                pathIntegral += hePath[k] * sigma[k];
            }
            model.addConstr(pathIntegral == period * edgePathConstraints[i].second);
        }

        //compute a piecewise linear function over the vertices of the mesh 
        std::vector<GRBLinExpr> u(gluedMesh.nVertices());
        std::vector<std::vector<GRBLinExpr>> gradU(gluedMesh.nFaces(), std::vector<GRBLinExpr>(3));
        for (Face f : gluedMesh.faces()){

            u[f.halfedge().vertex().getIndex()] = 0.0;
            u[f.halfedge().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()]) - (nP[f.getIndex()]/3.0);
            u[f.halfedge().next().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()] + sigma[f.halfedge().next().getIndex()]) - ((2.0 * nP[f.getIndex()])/3.0);

            GRBLinExpr currGradU = 0.0;
            //GRBVar currGradU = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][2] = currGradU;
        }

        GRBQuadExpr obj;
        std::vector<std::array<double, 3>> grads;
        //store the per face difference 
        std::vector<GRBQuadExpr> difference(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            //objective is (||\delta \sigma||^2 - 1)
            //need to sort out how to actually implement the above objective
            GRBQuadExpr Xterm = (gradU[f.getIndex()][0] * gradU[f.getIndex()][0]);
            GRBQuadExpr Yterm = (gradU[f.getIndex()][1] * gradU[f.getIndex()][1]);
            GRBQuadExpr Zterm = (gradU[f.getIndex()][2] * gradU[f.getIndex()][2]);
            // GRBVar Xterm = (gradU[f.getIndex()][0] * gradU[f.getIndex()][0]);
            // GRBVar Yterm = (gradU[f.getIndex()][1] * gradU[f.getIndex()][1]);
            // GRBVar Zterm = (gradU[f.getIndex()][2] * gradU[f.getIndex()][2]);
            double currArea = gluedGeometry.faceAreas[f];
            obj +=  currArea * (Xterm + Yterm + Zterm - 1);
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        std::cout << "Objective after placing " << singularEdges.size()/2 << " singularity pairs: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }     
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    std::cout << "--------------------------------" << std::endl;
    return std::tie(gluedOneForm, objectiveVal);

}

//find the isoval with max average deviation from \frac{\nabla h}{||h||} in the vertex setting
double findIsoValWithMaxVertexCurlDeviation(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& gradSigmaTilde, 
                                  FaceData<Vector3>& globalFaceGradients, VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh; 

    for (int i = 0; i < usedIsoVals.size(); i++){
        std::cout << "used iso vals = " << usedIsoVals[i] << std::endl;
    }

    double start = 0.1;
    double stepSize = 0.1;
    double end = 1.0;
    double curr = start;
    double maxAvgDeviation = -DBL_MAX;
    double maxAvgDeviationIsoVal = 0.0;
    double eps = 1e-6;
    bool skipFlag = false;
    globalGeometry.requireFaceAreas();
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);

    while (curr < end){
        for (int i = 0; i < usedIsoVals.size(); i++){
            //don't isoVals you've used before
            if (std::fabs(usedIsoVals[i] - curr) < eps){
                skipFlag = true;
            }
        }
        if (skipFlag){
            curr += stepSize;
            continue;
        }
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, curr);

        double currAvgDeviation = 0.0;
        double currDeviationSum = 0.0;
        for (int i = 0; i < f.size(); i++){
            currDeviationSum += globalGeometry.faceAreas[globalMesh.face(f[i])] * norm(gradSigmaTilde[globalMesh.face(f[i])] - globalFaceGradients[globalMesh.face(f[i])]); 
        }
        currAvgDeviation = currDeviationSum / f.size();
        if (currAvgDeviation > maxAvgDeviation){
            maxAvgDeviation = currAvgDeviation;
            maxAvgDeviationIsoVal = curr;
        }
        skipFlag = false;
        curr += stepSize;
    }
    return maxAvgDeviationIsoVal;

}

//find the isoval with max average curl in the edge setting
double findIsoValWithMaxAvgEdgeCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, EdgeData<double>& curl, 
                                            VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, 
                                            std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh; 

    for (int i = 0; i < usedIsoVals.size(); i++){
    //    std::cout << "used iso vals = " << usedIsoVals[i] << std::endl;
    }

    double stepSize = 0.01;
    double end = 1.0;
    double curr = stepSize;
    double maxDeviation = -DBL_MAX;
    double maxDeviationIsoVal = -1.0;
    double eps = 1e-8;
    bool skipFlag = false;
    globalGeometry.requireFaceAreas();
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);
    double currAvgDeviation = 0.0;
    double currDeviationSum = 0.0;

    while (curr < end){
        for (int i = 0; i < usedIsoVals.size(); i++){
            //don't isoVals you've used before
            if (std::fabs(usedIsoVals[i] - curr) < eps){
            //    std::cout << "skipping isoval " << usedIsoVals[i] << std::endl;
                skipFlag = true;
            }
        }
        if (skipFlag){
            curr += stepSize;
            skipFlag = false;
            continue;
        }
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, curr);
        //reset values
        currAvgDeviation = 0.0;
        currDeviationSum = 0.0;
        int numEdges = 0;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            for (Edge e : currFace.adjacentEdges()){
                //check if current isoline crosses this edge
                if ((curr > globalTimeFunction[e.halfedge().tailVertex()] && curr < globalTimeFunction[e.halfedge().tipVertex()])
                    ||(curr > globalTimeFunction[e.halfedge().tipVertex()] && curr < globalTimeFunction[e.halfedge().tailVertex()])){
                    currDeviationSum += std::fabs(curl[e]);
                }
            }
        }
        //std::cout << "Sum for isoval " << curr << " is " << currDeviationSum << std::endl;
        if (currDeviationSum > maxDeviation){
            maxDeviation = currDeviationSum;
            maxDeviationIsoVal = curr;
        }

        curr += stepSize;
    }

    // std::cout << "maxDeviation " << maxDeviation << std::endl;
    // std::cout << "maxDeviationIsoVal " << maxDeviationIsoVal << std::endl;
    return maxDeviationIsoVal;
}

//-------------------------Experiment 2-----------------------------//
//optimizing min \sum_{{ij} \in E}||\omegaTilde_{ij} - \sigmaTild_{ij}||^2 + \lambda \sum{e \in E}||\sigma_{ij} + \sigma_{ji}||^2
//visualizing where ||\sigma_{ij} + \sigma_{ji}||^2 is highest 
void vizEdgeDifference(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                        FaceData<Vector3>& globalFaceGradients, polyscope::SurfaceMesh& psMesh, 
                        globalBoundaryConditions& boundaryConditions, double period, double lambda,
                        std::map<int,int>& globalToGluedEdgeMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    globalGeometry.requireVertexPositions();
    std::vector<double> gluedOmegaTilde(gluedMesh.nHalfedges());
    std::vector<double> globalOmegaTilde(globalMesh.nHalfedges());
    for (Face f : globalMesh.faces()){
        Vector3 grad = globalFaceGradients[f];
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = globalGeometry.vertexPositions[he.tipVertex()] - globalGeometry.vertexPositions[he.tailVertex()];
            //omega tilde in the global setting
            globalOmegaTilde[he.getIndex()] = dot(heVec, grad);
        }
    }

    for (Face f : gluedMesh.faces()){
        gluedOmegaTilde[f.halfedge().getIndex()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge().getIndex()];
        gluedOmegaTilde[f.halfedge().next().getIndex()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge().next().getIndex()];
        gluedOmegaTilde[f.halfedge().next().next().getIndex()] = globalOmegaTilde[globalMesh.face(f.getIndex()).halfedge().next().next().getIndex()];
    }
    
    EdgeData<double> gluedEdgeDiff(gluedMesh);
    EdgeData<double> globalEdgeDiff(globalMesh);

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        //defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : boundaryConditions.courseBdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }
    
        //compute nP over every face 
        //add the second constraint while we're here
        //second constraint - (d1*sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            nP[f.getIndex()] = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0.0);
        }

        //set up the difference
        GRBQuadExpr obj = 0;

        for (Halfedge he : gluedMesh.halfedges()){
            obj += (gluedOmegaTilde[he.getIndex()] - sigma[he.getIndex()]) * (gluedOmegaTilde[he.getIndex()] - sigma[he.getIndex()])
                    + lambda * ((sigma[he.getIndex()] + sigma[he.twin().getIndex()]) * (sigma[he.getIndex()] + sigma[he.twin().getIndex()]));
        }
        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 

        for (Edge e : gluedMesh.edges()){
            gluedEdgeDiff[e] = sigma[e.halfedge().getIndex()].get(GRB_DoubleAttr_X) + sigma[e.halfedge().twin().getIndex()].get(GRB_DoubleAttr_X);
        }

    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    globalEdgeDiff = convertGluedToGlobalEdgeFunction(globalGeometry, gluedGeometry, gluedEdgeDiff, globalToGluedEdgeMap);

    psMesh.addEdgeScalarQuantity("edge difference", globalEdgeDiff);

}

//------------------------Experiment 3----------------------------//
//Trying to find a harmonic 1-form in the halfedge optimization setting 
//use the thus far inserted singular edges (and edge-based curl of zero elsewhere), 
//enforce face-based curl of zero, and an integral of 1 along some path from either boundary to optimize 
//for a harmonic (half-edge) one-form subject to the singularity placements; do this at each 
//iteration and then normalize to find edge-based curl
std::tuple<HalfedgeData<double>, EdgeData<double>> harmonic1FormImpl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, 
                                                                    polyscope::SurfaceMesh& psMesh, globalBoundaryConditions& boundaryConditions, double period){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    //show the level set of the time function 
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);

    //curl per edge
    EdgeData<double> edgeCurl(globalMesh);
    //curl per face 
    FaceData<double> faceCurl(globalMesh);
    //striping 1-form quantities
    HalfedgeData<double> oldGluedSigmaTilde(gluedMesh);
    HalfedgeData<double> newGluedSigmaTilde(gluedMesh);
    //edge singularities 
    EdgeData<double> edgeSingularities(globalMesh);
    FaceData<double> faceSingularities(globalMesh);
    //iso values we've placed singularities
    std::vector<double> usedIsoVals;
    //singular edges in the gurobi optimization
    std::vector<std::pair<int, int>> singularEdges;
    //face indices in the gurobi optimization
    std::vector<int> faceIndices(gluedMesh.nFaces(), 0);
    //objective values
    double oldObj, currObj;
    //number of singularity pairs  
    int numPairs = 0;
    //max number of singularity pairs to insert 
    int maxPairs = 6;
    //gurobi model we will be solving 
    Model model;

    //striping information 
    //global data
    CornerData<double> stripeValuesSigmaCourse;
    //global data
    FaceData<int> stripeIndicesSigmaCourse;
    std::vector<Vector3> positionsCourse;
    std::vector<std::array<int, 2>> edgesCourse;
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints; 
    //require edge lengths 
    globalGeometry.requireEdgeLengths();
    for (int i = 0; i < boundaryConditions.bdyBdyPathConstraints.size(); i++){
        //visualizing bdy-bdy edge constraints
        EdgeData<double> bdyBdyPath(globalMesh);
        double pathLength = 0.;
        for (Edge e : globalMesh.edges()){
            bdyBdyPath[e] = boundaryConditions.bdyBdyPathConstraints[i][globalToGluedEdgeMap[e.getIndex()]];
            if (std::fabs(bdyBdyPath[e]) > 0 && !e.isBoundary()){
                pathLength += globalGeometry.edgeLengths[e];
            } 
        }
        psMesh.addEdgeScalarQuantity("bdy bdy path " + std::to_string(i), bdyBdyPath);
        edgePathConstraints.push_back(std::make_pair(boundaryConditions.bdyBdyPathConstraints[i], (pathLength/period)));
    }
    model.setPeriod(period);
    model.setBdyEdges(boundaryConditions.courseBdyEdges);
    model.setEdgePathConstraints(edgePathConstraints);
    model.setFaceIndices(faceIndices);

    //-------------PLACING SINGULARITY ON EDGES-----------------//
    //compute curl per edge of the gradient field (in the global setting)
    edgeCurl = computeEdgeCurl(globalGeometry, gluedGeometry,
                                globalFaceGradients, 
                                globalToGluedEdgeMap);

    psMesh.addEdgeScalarQuantity("edge curl after placing " + std::to_string(numPairs) + " singularity pairs", edgeCurl);
    faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, gluedGeometry, globalFaceGradients, globalToGluedEdgeMap);
    psMesh.addFaceScalarQuantity("face curl after placing " + std::to_string(numPairs) + " singularity pairs", faceCurl);
    //find max/min curl face over entire mesh 
    std::pair<int, int> singFacePair = findFaceSingularityPair(globalGeometry, gluedGeometry, faceCurl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, 0.0, numPairs, true);
    faceSingularities[singFacePair.first] = 1.0;
    faceSingularities[singFacePair.second] = -1.0;
    faceIndices[singFacePair.first] = 1.0;
    faceIndices[singFacePair.second] = -1.0;
    int negEdge = findSingularEdgeFromSingularFace(globalGeometry, gluedGeometry, singFacePair.first, globalFaceGradients[singFacePair.first]);
    int posEdge = findSingularEdgeFromSingularFace(globalGeometry, gluedGeometry, singFacePair.second, globalFaceGradients[singFacePair.second]);
    edgeSingularities[globalMesh.edge(posEdge)] = 1.0;
    edgeSingularities[globalMesh.edge(negEdge)] = -1.0;
    //first pair of singular edges (edge-curl approach)
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[negEdge], -1));
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[posEdge], 1));
    std::cout << "negative edge " << negEdge << std::endl;
    std::cout << "positive edge " << posEdge << std::endl;
    
    //solve the model without any singularities 
    std::tie(newGluedSigmaTilde, currObj) = computeHarmonic1Form(globalGeometry, gluedGeometry, model, globalToGluedVertexMap, edgeMappingsPairs, psMesh);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
    std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);
    model.setSingularEdges(singularEdges);
    numPairs++;
    psMesh.addEdgeScalarQuantity("edge singularities after inserting " + std::to_string(numPairs) + " singularity pairs", edgeSingularities);
    //solve the model with 1 pair of singularities
    std::tie(newGluedSigmaTilde, currObj) = computeHarmonic1Form(globalGeometry, gluedGeometry, model, globalToGluedVertexMap, edgeMappingsPairs, psMesh);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);

    
    //hard-coding stopping conditions for now
    while(numPairs < maxPairs){
        //gradient of sigmaTilde in the global setting 
        FaceData<Vector3> gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, oldGluedSigmaTilde);
        psMesh.addFaceVectorQuantity("grad sigmaTilde after placing " + std::to_string(numPairs) + " singularity pairs", gradSigmaTilde);
        //recompute edge curl using gradSigmaTilde
        edgeCurl = computeEdgeCurl(globalGeometry, gluedGeometry, gradSigmaTilde, globalToGluedEdgeMap);
        psMesh.addEdgeScalarQuantity("edge curl after placing " + std::to_string(numPairs) + " singularity pairs", edgeCurl);
        faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, gluedGeometry, gradSigmaTilde, globalToGluedEdgeMap);
        psMesh.addFaceScalarQuantity("face curl after placing " + std::to_string(numPairs) + " singularity pairs", faceCurl);
        double isoVal = findIsoValWithMaxAvgFaceCurl(globalGeometry, gluedGeometry, faceCurl, gluedTimeFunction,
                                                    usedIsoVals, globalToGluedVertexMap);
        if (std::fabs(isoVal - (-1.0) < 1e-15)){
            std::cout << "Breaking cause we can't find any more sensible isovalues " << std::endl;
            break;//our function couldn't find any more isovalues
        } 
        std::cout << "next isoval = " << isoVal << std::endl;
        // singEdgePair = findEdgeSingularityPair(globalGeometry, gluedGeometry, edgeCurl,
        //                                 gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, isoVal, numPairs, false);
        //find max/min curl face 
        singFacePair = findFaceSingularityPair(globalGeometry, gluedGeometry, faceCurl,
                                            gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, isoVal, numPairs, false);
        if (faceIndices[singFacePair.first] != 0
            || faceIndices[singFacePair.second] != 0){
            //don't use faces you've used before
            usedIsoVals.push_back(isoVal);
            continue;
        }
        faceSingularities[singFacePair.first] = 1.0;
        faceSingularities[singFacePair.second] = -1.0;
        faceIndices[singFacePair.first] = 1.0;
        faceIndices[singFacePair.second] = -1.0;
        int negEdge = findSingularEdgeFromSingularFace(globalGeometry, gluedGeometry, singFacePair.first, globalFaceGradients[singFacePair.first]);
        int posEdge = findSingularEdgeFromSingularFace(globalGeometry, gluedGeometry, singFacePair.second, globalFaceGradients[singFacePair.second]);
        //don't select edges we've seen before 
        if ((std::find(singularEdges.begin(), singularEdges.end(), std::make_pair(globalToGluedEdgeMap[posEdge], 1)) != singularEdges.end())
            || std::find(singularEdges.begin(), singularEdges.end(), std::make_pair(globalToGluedEdgeMap[negEdge], -1)) != singularEdges.end()){
            //don't use edges you've used before
            usedIsoVals.push_back(isoVal);
            continue;
        }
        edgeSingularities[globalMesh.edge(posEdge)] = 1.0;
        edgeSingularities[globalMesh.edge(negEdge)] = -1.0;
        //first pair of singular edges (edge-curl approach)
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[negEdge], -1));
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[posEdge], 1));
        usedIsoVals.push_back(isoVal);
        // singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.first], -1));
        // singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singEdgePair.second], 1));
        model.setSingularEdges(singularEdges);
        //solve the optimization problem 
        //sigmaTilde is in the glued mesh setting 
        std::tie(newGluedSigmaTilde, currObj) = computeHarmonic1Form(globalGeometry, gluedGeometry, model, globalToGluedVertexMap, edgeMappingsPairs, psMesh);
        // if (currObj > oldObj){//current objective become worse or staying the same 
        //     std::cout << "Breaking cause objective no longer improving..." << std::endl;
        //     std::cout << "oldObj " << oldObj << std::endl;
        //     std::cout << "currObj " << currObj << std::endl;
        //     break;
        // }
        //increment, update and visualize
        //edgeSingularities[globalMesh.edge(singEdgePair.first)] = 1.0;
        //edgeSingularities[globalMesh.edge(singEdgePair.second)] = -1.0;
        numPairs++;
        psMesh.addFaceScalarQuantity("face singularities after inserting " + std::to_string(numPairs) + " singularity pairs", faceSingularities);
        psMesh.addEdgeScalarQuantity("edge singularities after inserting " + std::to_string(numPairs) + " singularity pairs", edgeSingularities);
        oldObj = currObj;
        oldGluedSigmaTilde = newGluedSigmaTilde;
        std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
        courseStripes -> setRadius(0.001);
        courseStripes -> setEnabled(false);
    }

    return std::tie(oldGluedSigmaTilde, edgeSingularities);
    
    //----------------------------------------------------//

    //------------PLACING SINGULARITIES AT FACES------------//
    // faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, gluedGeometry, globalFaceGradients, globalToGluedEdgeMap);
    // psMesh.addFaceScalarQuantity("face curl after placing " + std::to_string(numPairs) + " singularity pairs", faceCurl);

    // //solve the model without any singularities 
    // std::tie(newGluedSigmaTilde, currObj) = computeHarmonic1Form(globalGeometry, gluedGeometry, model, globalToGluedVertexMap, psMesh);
    // oldObj = currObj;
    // oldGluedSigmaTilde = newGluedSigmaTilde;
    // std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
    // std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
    // auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
    //                                                 " singularity pairs", positionsCourse, edgesCourse);
    // courseStripes -> setRadius(0.001);
    // courseStripes -> setEnabled(false);
    // FaceData<Vector3> gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, oldGluedSigmaTilde);
    // psMesh.addFaceVectorQuantity("grad Sigma after placing " + std::to_string(numPairs) + " singularity pairs", gradSigmaTilde);
    
    // //add the first pair of singularities
    // //find max/min curl face over entire mesh 
    // std::pair<int, int> singFacePair = findFaceSingularityPair(globalGeometry, gluedGeometry, faceCurl,
    //                                     gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, 0.0, numPairs, true);
    // faceSingularities[singFacePair.first] = 1.0;
    // faceSingularities[singFacePair.second] = -1.0;
    // faceIndices[singFacePair.first] = 1.0;
    // faceIndices[singFacePair.second] = -1.0;
    // numPairs++;
    // model.setFaceIndices(faceIndices);
    // //solve the model with one pair of singularities 
    // std::tie(newGluedSigmaTilde, currObj) = computeHarmonic1Form(globalGeometry, gluedGeometry, model, globalToGluedVertexMap, psMesh);
    // oldObj = currObj;
    // oldGluedSigmaTilde = newGluedSigmaTilde;
    // std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
    // std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
    // courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
    //                                                 " singularity pairs", positionsCourse, edgesCourse);
    // courseStripes -> setRadius(0.001);
    // courseStripes -> setEnabled(false);
    // gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, oldGluedSigmaTilde);
    // psMesh.addFaceVectorQuantity("grad Sigma after placing " + std::to_string(numPairs) + " singularity pairs", gradSigmaTilde);
    // faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, gluedGeometry, gradSigmaTilde, globalToGluedEdgeMap);
    // psMesh.addFaceScalarQuantity("face curl after placing " + std::to_string(numPairs) + " singularity pairs", faceCurl);
    // return std::tie(oldGluedSigmaTilde, edgeSingularities);
}

//solve the optimization problem for the harmonic 1-form
std::tuple<HalfedgeData<double>, double> computeHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& globalToGluedVertexMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, polyscope::SurfaceMesh& psMesh){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);
    FaceData<Vector3> normalizedDelSigma(gluedMesh);
    EdgeData<double> edgeAveragedDelSigma(globalMesh);
    FaceData<double> d1edgeAveragedDelSigma(globalMesh);
    double objectiveVal;
    int numSingularFaces = 0;

    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<int> faceIndices = gbModel.getFaceIndices();
    

    //build the gradient operator 
    Eigen::MatrixXd V(gluedMesh.nVertices(), 3);
    Eigen::MatrixXi F(gluedMesh.nFaces(), 3);
    std::tie(V, F) = getVertexPositionsandFaceLists(globalGeometry);
    // Compute the global gradient operator: #F*3 by #V
    Eigen::SparseMatrix<double> grad;
    igl::grad(V,F,grad);
    Eigen::SparseMatrix<double, Eigen::RowMajor> G(grad);


    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();

    // Eigen::SparseMatrix<double, Eigen::RowMajor> d0 = gluedGeometry.d0;
    // Eigen::SparseMatrix<double, Eigen::RowMajor> d1 = gluedGeometry.d1;
    // Eigen::SparseMatrix<double, Eigen::RowMajor> hodge1 = gluedGeometry.hodge1;
    // Eigen::SparseMatrix<double, Eigen::RowMajor> hodge1Inverse = gluedGeometry.hodge1Inverse;
    Eigen::SparseMatrix<double, Eigen::RowMajor> hodge0Inverse = gluedGeometry.hodge0Inverse;
    Eigen::SparseMatrix<double, Eigen::RowMajor> hodge2 = gluedGeometry.hodge2;
    //Eigen::SparseMatrix<double, Eigen::RowMajor> term1 = (hodge1Inverse * d1.transpose() * hodge2 * d1);
    //Eigen::SparseMatrix<double, Eigen::RowMajor> term2 = (d0 * hodge0Inverse * d0.transpose() * hodge1);
    //Eigen::SparseMatrix<double, Eigen::RowMajor> OneLaplacian = term1 + term2;
    // std::cout << "Number of rows in the one laplacian: " << OneLaplacian.rows() << std::endl;
    // std::cout << "Number of cols in the one laplacian: " << OneLaplacian.cols() << std::endl;
    // std::cout << "Number of edges in the mesh: " << gluedMesh.nEdges() << std::endl;

    //implement d1 to have dimensions |F| x |HE| 
    Eigen::SparseMatrix<double, Eigen::RowMajor> d1(gluedMesh.nFaces(), gluedMesh.nHalfedges());
    for (Face f : gluedMesh.faces()){
        for (Halfedge he : f.adjacentHalfedges()){
            d1.coeffRef(f.getIndex(), he.getIndex()) = 1.0;
        }
    }
    //implement d0 to have dimension |HE| x |V|
    Eigen::SparseMatrix<double, Eigen::RowMajor> d0(gluedMesh.nHalfedges(), gluedMesh.nVertices());
    for (Halfedge he : gluedMesh.halfedges()){
        int sourceI = he.tailVertex().getIndex();
        int targetI = he.tipVertex().getIndex();
        d0.coeffRef(he.getIndex(), sourceI) = -1.0;
        d0.coeffRef(he.getIndex(), targetI) = 1.0;
    }
    
    // //implement the hogde1 as |HE| x |HE| 
    Eigen::SparseMatrix<double, Eigen::RowMajor> hodge1(gluedMesh.nHalfedges(), gluedMesh.nHalfedges());
    Eigen::SparseMatrix<double, Eigen::RowMajor> hodge1Inverse(gluedMesh.nHalfedges(), gluedMesh.nHalfedges());
    Eigen::VectorXd hodge1V(gluedMesh.nHalfedges());
    for (Edge e : gluedMesh.edges()) {
      double ratio = gluedGeometry.edgeCotanWeights[e];
      hodge1V[e.halfedge().getIndex()] = ratio;
      hodge1V[e.halfedge().twin().getIndex()] = ratio;
    }
    hodge1 = hodge1V.asDiagonal();
    hodge1Inverse = hodge1V.asDiagonal().inverse();

    Eigen::SparseMatrix<double, Eigen::RowMajor> term1 = (hodge1Inverse * d1.transpose() * hodge2 * d1);
    Eigen::SparseMatrix<double, Eigen::RowMajor> term2 = (d0 * hodge0Inverse * d0.transpose() * hodge1);
    Eigen::SparseMatrix<double, Eigen::RowMajor> OneLaplacian = term1 + term2;

    // Eigen::LLT<Eigen::MatrixXd> lltOfOneLaplacian(OneLaplacian); // compute the Cholesky decomposition of A
    // if(lltOfOneLaplacian.info() == Eigen::NumericalIssue)
    // {
    //     throw std::runtime_error("One Laplacian is possibly non semi-positive definitie matrix!");
    // }    

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 0.5);
        //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_SolutionLimit, 2);

        //sigma defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //add integer variable for bdy-bdy path constraints 
        std::vector<GRBVar> bdyBdyPathIntegers;
        for (size_t i = 0; i < edgePathConstraints.size(); i++){
            GRBVar k = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            bdyBdyPathIntegers.push_back(k);
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        //compute nP over every face 
        //add the second constraint while we're here
        //second constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (int r = 0; r < d1.outerSize(); ++r) {
            GRBLinExpr lhs = 0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d1, r); it; ++it ) {
                lhs += it.value() * sigma[it.col()];
            }
            nP[r] = lhs;
            if (faceIndices[r] > std::fabs(1e-8)) numSingularFaces++;
            model.addConstr(lhs == faceIndices[r], "Integral Constraint");
        }

        //third constraint 
        //specify singular halfedges and also specify that form values 
        //are equal across non-singular halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            //invert the sign of the constraint
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] == -1.0 * p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //fourth constraint: add bdy-bdy path constraint 
        for (int i = 0; i < edgePathConstraints.size(); i++){
            std::vector<double> path = edgePathConstraints[i].first;
            GRBLinExpr pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(j).halfedge().getIndex()] = std::fabs(path[j]);
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = std::fabs(path[j]);
                }
            }
            for (int k = 0; k < gluedMesh.nHalfedges(); k++){
                pathIntegral += hePath[k] * sigma[k];
            }
            if (std::fabs(edgePathConstraints[i].second) > 0)
                model.addConstr(pathIntegral == period * edgePathConstraints[i].second);
            else
                model.addConstr(pathIntegral == period * bdyBdyPathIntegers[i]);
        }

        //compute a piecewise linear function over the vertices of the mesh 
        std::vector<GRBLinExpr> u(gluedMesh.nVertices());
        std::vector<std::vector<GRBLinExpr>> gradU(gluedMesh.nFaces(), std::vector<GRBLinExpr>(3));
        for (Face f : gluedMesh.faces()){

            u[f.halfedge().vertex().getIndex()] = 0.0;
            u[f.halfedge().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()]) - (nP[f.getIndex()]/3.0);
            u[f.halfedge().next().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()] + sigma[f.halfedge().next().getIndex()]) - ((2.0 * nP[f.getIndex()])/3.0);

            GRBLinExpr currGradU = 0.0;
            //GRBVar currGradU = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * F.rows()); it; ++it){
                currGradU += it.value() * u[globalToGluedVertexMap[it.col()]];
            }
            gradU[f.getIndex()][2] = currGradU;
        }


        //set up the objective term
        GRBQuadExpr obj = 0;        
        //setting the objective to be min cot_e||\sigma||^2
        for (Halfedge he : gluedMesh.halfedges()){
            obj +=  gluedGeometry.edgeCotanWeights[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
            // obj += gluedGeometry.edgeLengths[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
        }

        //setting the objective to be min ||\nabla \sigma||^2
        //weighted by the face areas
        // for (Face f : gluedMesh.faces()){
        //     obj +=  (gradU[f.getIndex()][0] * gradU[f.getIndex()][0] + gradU[f.getIndex()][1] * gradU[f.getIndex()][1]
        //             + gradU[f.getIndex()][2] * gradU[f.getIndex()][2]);
        // }

        //setting the objective to be min ||\sigma^T \Delta^1 \sigma||
        // std::vector<GRBLinExpr> delSigma(gluedMesh.nHalfedges()); 
        // for (int r = 0; r < OneLaplacian.outerSize(); ++r) {
        //     GRBLinExpr lhs = 0;
        //     for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(OneLaplacian, r); it; ++it ) {
        //         lhs += it.value() * sigma[it.col()];
        //     }
        //     delSigma[r] = lhs;
        // }
        // for (Halfedge he : gluedMesh.halfedges()){
        //     obj += sigma[he.getIndex()] * delSigma[he.getIndex()];
        // }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }

        model.write("myModel.lp");

        //print out value of edge path constraints
        for (int i = 0; i < edgePathConstraints.size(); i++){
            std::vector<double> path = edgePathConstraints[i].first;
            double pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(j).halfedge().getIndex()] = std::fabs(path[j]);
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = std::fabs(path[j]);
                }
            }
            for (int k = 0; k < gluedMesh.nHalfedges(); k++){
                pathIntegral += hePath[k] * sigma[k].get(GRB_DoubleAttr_X);
            }
            std::cout << "value of path integral " << i << " is: " << pathIntegral << std::endl;
        } 

        //store normalized \del sigma as a vector per face
        for (Face f : gluedMesh.faces()){
            normalizedDelSigma[f] = Vector3{gradU[f.getIndex()][0].getValue(), gradU[f.getIndex()][1].getValue(), gradU[f.getIndex()][2].getValue()};
            normalizedDelSigma[f] = normalizedDelSigma[f].normalize();
        }

        //compute edge-averaged del sigma 
        //create a map from the mapped edges
        std::map<int, int> edgeMap;
        for (std::pair<int, int> pair : edgeMappingsPairs){
            edgeMap.insert({pair.first, pair.second});
        }

        //edges which we've handles already
        std::map<int, bool> seenEdges;
        for (Edge e : globalMesh.edges()){
            seenEdges.insert({e.getIndex(), false});
        }

        for (Edge e : globalMesh.edges()){
            if (seenEdges[e.getIndex()]) continue;
            if (e.halfedge().twin().isInterior()){//found an interior halfedge
                Vector3 eVector = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
                edgeAveragedDelSigma[e] = 0.5 * dot((normalizedDelSigma[e.halfedge().face()] + normalizedDelSigma[e.halfedge().twin().face()]),
                                    eVector);
                seenEdges[e.getIndex()] = true;
            }
            else{//found a boundary halfedge
                if (edgeMap.find(e.getIndex()) != edgeMap.end()){//found a stitched together edge
                    //take the average direction vector? 
                    Vector3 e1 = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
                    //just pick the original edge as the "canonical" direction in the global mesh
                    edgeAveragedDelSigma[e] = 0.5 * dot((normalizedDelSigma[e.halfedge().face()] + normalizedDelSigma[e.halfedge().twin().face()]), e1);
                    edgeAveragedDelSigma[globalMesh.edge(edgeMap.at(e.getIndex()))] = 0.5 * dot((normalizedDelSigma[e.halfedge().face()] + normalizedDelSigma[e.halfedge().twin().face()]), e1);
                    seenEdges[e.getIndex()] = true;
                    seenEdges[edgeMap.at(e.getIndex())] = true;
                }
                else{//found a boundary edge that's not stitched to anything
                    Vector3 eVector = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
                    edgeAveragedDelSigma[e] = dot(normalizedDelSigma[e.halfedge().face()], eVector);
                    seenEdges[e.getIndex()] = true;
                }
            }
        }

        //compute curl use edge-averaged del sigma
        for (Face f : globalMesh.faces()){
            int signIJ = f.halfedge().orientation() ? 1 : -1;
            int signJK = f.halfedge().next().orientation() ? 1 : -1;
            int signKI = f.halfedge().next().next().orientation() ? 1 : -1;
            double valIJ = signIJ * edgeAveragedDelSigma[f.halfedge().edge()];
            double valJK = signJK * edgeAveragedDelSigma[f.halfedge().next().edge()];
            double valKI = signKI * edgeAveragedDelSigma[f.halfedge().next().next().edge()];
            d1edgeAveragedDelSigma[f] = valIJ + valJK + valKI;
        }
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    psMesh.addFaceScalarQuantity("d1 * edge averaged del sigma after placing " + std::to_string(singularEdges.size() / 2) + " singular edges", d1edgeAveragedDelSigma);
    //psMesh.addFaceVectorQuantity("normalized del sigma after placing " + std::to_string(numSingularFaces) + " singular edges", normalizedDelSigma);
    return std::tie(gluedOneForm, objectiveVal);
}

//compute face curl by averaging edge curl over the edges in a face 
FaceData<double> computeAverageEdgeCurlonFaces(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                FaceData<Vector3>& globalFaceGradients, std::map<int, int>& globalToGluedEdgeMap){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    FaceData<double> faceCurl(globalMesh);
    EdgeData<double> edgeCurl = computeEdgeCurl(globalGeometry, gluedGeometry, globalFaceGradients, globalToGluedEdgeMap);

    for (Face f : globalMesh.faces()){
        double sum = 0.0;
        for (Edge e : f.adjacentEdges()){
            sum += edgeCurl[e];
        }
        faceCurl[f] = sum/3.;
    }
    return faceCurl;   
}

//find max/min curl face for a given isoline of the TIME FUNCTION 
//would probably want to change tracing the level sets of the time function with level sets of the harmonic 1-form
std::pair<int, int> findFaceSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, int numPairs, bool useAllFaces){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);


    if (useAllFaces){//find max/min curl over the pairs of all isolines over the mesh
        int maxFace, minFace = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        //find max curl over entire mesh 
        //find face with max  curl
        for (Face f : gluedMesh.faces()){
            if (curl[f] > maxCurl){
                maxCurl = curl[f];
                maxFace = f.getIndex();
            }
        }
        //average isovalue on that face
        double isoVal = (gluedTimeFunction[gluedMesh.face(maxFace).halfedge().tailVertex()] + gluedTimeFunction[gluedMesh.face(maxFace).halfedge().next().tailVertex()]
                            + gluedTimeFunction[gluedMesh.face(maxFace).halfedge().next().next().tailVertex()])/3.0;
        usedIsoVals.push_back(isoVal);
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        isoline -> setRadius(0.001);
        isoline -> setEnabled(false);
        //find min curl face on the same isoline as the max edge 
        //find face with min curl
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            if (curl[currFace] < minCurl){
                minCurl = curl[currFace];
                minFace = currFace.getIndex();
            }
        }
        return std::make_pair(maxFace, minFace);
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        isoline -> setEnabled(false);
        isoline -> setRadius(0.001);
        int maxFace, minFace = -1;
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            if (curl[currFace] > maxCurl){
                maxCurl = curl[currFace];
                maxFace = currFace.getIndex();
            }
            if (curl[currFace] < minCurl){
                minCurl = curl[currFace];
                minFace = currFace.getIndex();
            }
        }
        return std::make_pair(maxFace, minFace);
    }
}

int findSingularEdgeFromSingularFace(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, int singFaceIndex, Vector3 globalFaceGradient){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    double max = -DBL_MAX;
    Edge maxEdge;
    
    for (Edge e : globalMesh.face(singFaceIndex).adjacentEdges()){
        Halfedge he1 = e.halfedge();
        Halfedge he2 = e.halfedge().twin();
        double p1 = dot(globalGeometry.vertexPositions[he1.tipVertex()] - globalGeometry.vertexPositions[he1.tailVertex()], globalFaceGradient);
        double p2 = dot(globalGeometry.vertexPositions[he2.tipVertex()] - globalGeometry.vertexPositions[he2.tailVertex()], globalFaceGradient);
        if (p1 > max){
            max = p1;
            maxEdge = he1.edge();
        }
        else if (p2 > max){
            max = p2;
            maxEdge = he2.edge();
        }
    }
    std::cout << "max edge = " << maxEdge << std::endl;
    return maxEdge.getIndex();
}

//find the isoval with max average curl in the face setting
double findIsoValWithMaxAvgFaceCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<double>& curl, 
                                            VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, 
                                            std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh; 

    for (int i = 0; i < usedIsoVals.size(); i++){
        std::cout << "used iso vals = " << usedIsoVals[i] << std::endl;
    }

    double stepSize = 0.05;
    double end = 1.0;
    double curr = stepSize;
    double maxDeviation = -DBL_MAX;
    double maxDeviationIsoVal = -1.0;
    double eps = 1e-8;
    bool skipFlag = false;
    globalGeometry.requireFaceAreas();
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);
    double currAvgDeviation = 0.0;
    double currDeviationSum = 0.0;

    while (curr < end){
        for (int i = 0; i < usedIsoVals.size(); i++){
            //don't isoVals you've used before
            if (std::fabs(usedIsoVals[i] - curr) < eps){
            //    std::cout << "skipping isoval " << usedIsoVals[i] << std::endl;
                skipFlag = true;
            }
        }
        if (skipFlag){
            curr += stepSize;
            skipFlag = false;
            continue;
        }
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, curr);
        //reset values
        currAvgDeviation = 0.0;
        currDeviationSum = 0.0;
        int numEdges = 0;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            currDeviationSum += std::fabs(curl[currFace]);
        }

        if (currDeviationSum > maxDeviation){
            maxDeviation = currDeviationSum;
            maxDeviationIsoVal = curr;
        }
        curr += stepSize;
    }

    // std::cout << "maxDeviation " << maxDeviation << std::endl;
    // std::cout << "maxDeviationIsoVal " << maxDeviationIsoVal << std::endl;
    return maxDeviationIsoVal;
}

