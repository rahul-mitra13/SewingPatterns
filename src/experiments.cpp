#include "experiments.h"

//--------------------------Implementation of strategy 1---------------------------------------//
std::tuple<HalfedgeData<double>, VertexData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::map<int, int>& globalToGluedVertexMap, std::map<int, int>& globalToGluedEdgeMap,
                                                                    polyscope::SurfaceMesh& psMesh, globalBoundaryConditions& boundaryConditions, double period){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    HalfedgeData<double> oldGluedSigmaTilde(gluedMesh);
    HalfedgeData<double> newGluedSigmaTilde(gluedMesh);
    HalfedgeData<double> globalSigmaTilde(globalMesh);
    double currObj;
    VertexData<double> vertexSingularities(globalMesh, 0.0);
    //curl per vertex 
    VertexData<double> vertexCurl(globalMesh);
    //curl per edge
    EdgeData<double> edgeCurl(globalMesh);
    //singular vertex indices 
    std::vector<std::pair<int, int>> singularVertices;
    //singular edges
    std::vector<std::pair<int, int>> singularEdges;
    //iso values we've placed singularities
    std::vector<double> usedIsoVals;
    //tracking objective values 
    double oldObj;
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
    edgeCurl = computeEdgeCurl(globalGeometry, globalFaceGradients);



    //visualize the vertex curl 
    //this quantity never changes through the optimization 
    psMesh.addVertexScalarQuantity("vertex curl", vertexCurl);
    //visualize the edge curl 
    psMesh.addEdgeScalarQuantity("edge curl", edgeCurl);

    //find max/min vertex over entire mesh
    std::pair<int, int> singVertexPair = findVertexSingularityPair(globalGeometry, gluedGeometry, vertexCurl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, 0.0, numPairs, true);
    
    vertexSingularities[globalMesh.vertex(singVertexPair.first)] = 1.0;
    vertexSingularities[globalMesh.vertex(singVertexPair.second)] = -1.0;
    singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.first], 1.0));
    singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.second], -1.0));
    

    //make the model
    //we will solve it in the glued mesh setting 
    Model model;
    std::vector<std::array<double, 3>> faceGradients(gluedMesh.nFaces());
    for (Face f : globalMesh.faces()){
        faceGradients[f.getIndex()] = std::array{globalFaceGradients[f][0], globalFaceGradients[f][1], globalFaceGradients[f][2]};
    }
    model.setFaceGradients(faceGradients);
    model.setPeriod(period);
    model.setBdyEdges(boundaryConditions.courseBdyEdges);

    //solve the model without any singularities 
    std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;

    //find first pair of singular edges  
    std::pair<int, int> singularEdgeIndices = findSingularEdges(globalGeometry, gluedGeometry, singVertexPair, gluedOneRingMap, model);
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.first], 1));
    singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.second], -1));
    model.setSingularEdges(singularEdges);
    numPairs++;
    psMesh.addVertexScalarQuantity("vertex singularity after inserting " + std::to_string(numPairs) + " singularity pairs", vertexSingularities);
    //solve the model with 1 pair of singularities
    std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
    oldObj = currObj;
    oldGluedSigmaTilde = newGluedSigmaTilde;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
    courseStripes -> setRadius(0.001);

    while (true){
    
        //gradient of sigmaTilde in the global setting 
        FaceData<Vector3> gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, oldGluedSigmaTilde);
        psMesh.addFaceVectorQuantity("grad sigmaTilde after placing " + std::to_string(numPairs) + " singularity pairs", gradSigmaTilde);
        double isoVal = findIsoValWithMaxDeviation(globalGeometry, gluedGeometry, gradSigmaTilde, globalFaceGradients, gluedTimeFunction, 
                        usedIsoVals, globalToGluedVertexMap);
        // usedIsoVals.push_back(isoVal);
        singVertexPair = findVertexSingularityPair(globalGeometry, gluedGeometry, vertexCurl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, usedIsoVals, isoVal, numPairs, false);
        
        if (std::find(singularVertices.begin(), singularVertices.end(), std::make_pair(globalToGluedVertexMap[singVertexPair.first], 1.0)) != singularVertices.end()
            || std::find(singularVertices.begin(), singularVertices.end(), std::make_pair(globalToGluedVertexMap[singVertexPair.second], -1.0)) != singularVertices.end()){
            usedIsoVals.push_back(isoVal);
            continue;
        }
        usedIsoVals.push_back(isoVal);
        singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.first], 1.0));
        singularVertices.push_back(std::make_pair(globalToGluedVertexMap[singVertexPair.second], -1.0));
        singularEdgeIndices = findSingularEdges(globalGeometry, gluedGeometry, singVertexPair, gluedOneRingMap, model);
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.first], 1));
        singularEdges.push_back(std::make_pair(globalToGluedEdgeMap[singularEdgeIndices.second], -1));
        model.setSingularEdges(singularEdges);
        //solve the optimization problem
        //sigmaTilde is in the glued setting 
        std::tie(newGluedSigmaTilde, currObj) = computeStrategy1_oneForm(globalGeometry, gluedGeometry, model, globalToGluedVertexMap);
        if (currObj > oldObj){//current objective become worse or staying the same 
            std::cout << "Breaking cause objective no longer improving..." << std::endl;
            std::cout << "oldObj " << oldObj << std::endl;
            std::cout << "currObj " << currObj << std::endl;
            break;
        }
        
        vertexSingularities[globalMesh.vertex(singVertexPair.first)] = 1.0;
        vertexSingularities[globalMesh.vertex(singVertexPair.second)] = -1.0;

        //increment, update and visualize
        numPairs++;
        psMesh.addVertexScalarQuantity("vertex singularity after inserting " + std::to_string(numPairs) + " singularity pairs", vertexSingularities);
        oldObj = currObj;
        oldGluedSigmaTilde = newGluedSigmaTilde;
        std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, oldGluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
        courseStripes -> setRadius(0.001);
    }

    return std::tie(oldGluedSigmaTilde, vertexSingularities);

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
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    EdgeData<double> curl(globalMesh);

    for (Edge e : globalMesh.edges()){
        Vector3 heVec = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
        Vector3 heTwinVec = globalGeometry.vertexPositions[e.halfedge().tailVertex()] - globalGeometry.vertexPositions[e.halfedge().tipVertex()];
        //consider dividing by cotan weights, edge length, 1/3 * adjacent face areas etc. 
        curl[e] = dot(globalFaceGradients[e.halfedge().face()], heVec) + dot(globalFaceGradients[e.halfedge().twin().face()], heTwinVec);
    }

    return curl;

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
std::pair<int, int> findEdgeSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            std::vector<double>& usedIsoVals, double isoVal, bool useAllEdges){

    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
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
    std::vector<std::pair<int, int>> singularFaceIndices = gbModel.getSingularFaces();
    std::vector<int> faceIndices = gbModel.getFaceIndices();
    double period = gbModel.getPeriod();
    bool hasIntegrabilityConstraint = gbModel.getIntegrabilityConstraint();
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
                    hePath[gluedMesh.edge(i).halfedge().getIndex()] = path[j];
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(i).halfedge().twin().getIndex()] = path[j];
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

        //set up the difference
        GRBQuadExpr gradDiff = 0;
        std::vector<std::array<double, 3>> grads;
        //store the per face difference 
        std::vector<GRBQuadExpr> difference(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBQuadExpr diffX = (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]) * (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]);
            GRBQuadExpr diffY = (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]) * (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]);
            GRBQuadExpr diffZ = (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]) * (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]);
            double currArea = gluedGeometry.faceAreas[f];
            gradDiff += currArea * (diffX + diffY + diffZ);
            difference[f.getIndex()] = currArea * (diffX + diffY + diffZ);
        }

        GRBQuadExpr obj = gradDiff;

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

//find the isoval with max average deviation from \frac{\nabla h}{||h||}
double findIsoValWithMaxDeviation(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& gradSigmaTilde, 
                                  FaceData<Vector3>& globalFaceGradients, VertexData<double>& gluedTimeFunction, std::vector<double>& usedIsoVals, std::map<int, int>& globalToGluedVertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh; 

    //first normalize the global face gradients 
    for (Face f : globalMesh.faces()){
        globalFaceGradients[f] = globalFaceGradients[f].normalize();
    }
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

//-------------------------Experiment 2-----------------------------//
//optimizing min \sum_{{ij} \in E}||\omegaTilde_{ij} - \sigmaTild_{ij}||^2 + \lambda \sum{e \in E}||\sigma_{ij} + \sigma_{ji}||^2
//visualizing where ||\sigma_{ij} + \sigma_{ji}||^2 is highest 
void vizEdgeDifference(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                        FaceData<Vector3>& globalFaceGradients, polyscope::SurfaceMesh& psMesh, 
                        globalBoundaryConditions& boundaryConditions, double period, double lambda){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    globalGeometry.requireVertexPositions();
    std::vector<double> omegaTilde(globalMesh.nHalfedges());
    for (Face f : globalMesh.faces()){
        Vector3 grad = globalFaceGradients[f];
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = globalGeometry.vertexPositions[he.tipVertex()] - globalGeometry.vertexPositions[he.tailVertex()];
            omegaTilde[he.getIndex()] = dot(heVec, grad);
        }
    }
    EdgeData<double> edgeDiff(globalMesh);

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
            obj += (omegaTilde[he.getIndex()] - sigma[he.getIndex()]) * (omegaTilde[he.getIndex()] - sigma[he.getIndex()])
                    + lambda * ((sigma[he.getIndex()] + sigma[he.twin().getIndex()]) * (sigma[he.getIndex()] + sigma[he.twin().getIndex()]));
        }
        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 

        //put the computed one-form into and edge quantity
        for (Face f : gluedMesh.faces()){
            edgeDiff[globalMesh.face(f.getIndex()).halfedge().edge()] = sigma[f.halfedge().getIndex()].get(GRB_DoubleAttr_X) + sigma[f.halfedge().twin().getIndex()].get(GRB_DoubleAttr_X);
            edgeDiff[globalMesh.face(f.getIndex()).halfedge().next().edge()] = sigma[f.halfedge().next().getIndex()].get(GRB_DoubleAttr_X) + sigma[f.halfedge().next().twin().getIndex()].get(GRB_DoubleAttr_X);
            edgeDiff[globalMesh.face(f.getIndex()).halfedge().next().next().edge()] = sigma[f.halfedge().next().next().getIndex()].get(GRB_DoubleAttr_X) + sigma[f.halfedge().next().next().twin().getIndex()].get(GRB_DoubleAttr_X);
        }
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    psMesh.addEdgeScalarQuantity("edge difference", edgeDiff);

}