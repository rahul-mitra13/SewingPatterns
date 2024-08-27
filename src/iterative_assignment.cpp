#include "iterative_assignment.h"


EdgeData<double> omega;
Model gbModel;
std::vector<int> usedFaceIndices;
Eigen::SparseMatrix<double, Eigen::RowMajor> dOne; 
Eigen::SparseMatrix<double, Eigen::RowMajor> dZero;
std::map<int, int> globalToGluedVertexMap;

//----------------First strategy------------------------//
//This should reall happen in the intrinsic setting
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, VertexData<double>& timeFunction, EdgeData<double>& omega, 
                                            double period, std::map<int, int>& vertexMap, std::vector<std::pair<int, int>>& edgeMappingsPairs, globalBoundaryConditions& globalBdyConditions){

    SurfaceMesh& mesh = geometry.mesh; 
    FaceData<int> singPositions(mesh);
    //build differential operators in this stitched mesh setting
    geometry.requireDECOperators();
    dZero = geometry.d0;
    dOne = geometry.d1;
    //convert dOne to column major for easy updates of columns 
    Eigen::SparseMatrix<double, Eigen::ColMajor> d_oneColMajor;
    d_oneColMajor = dOne;
    //invert the signs to account for "stitched together" panels in d0
    for (std::pair<int, int> p : edgeMappingsPairs){
        int iE2 = p.second;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(dZero, iE2); it; ++it){
            it.valueRef() = -it.value();
        }
    }
    //invert the signs to account for "stitched together" panels
    for (std::pair<int, int> p : edgeMappingsPairs){
        int iE2 = p.second;
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(d_oneColMajor, iE2); it; ++it){
            it.valueRef() = -it.value();
        }
    }
    //specify boundary-boundary edge path constraint 
    //i.e., non-collapse constraint
    Vertex v0;
    Vertex v1; 
    int ctr = 0;
    int gluedStartVertexI = globalBdyConditions.courseStartBoundaryVertices[0];
    int gluedEndVertexI = globalBdyConditions.courseEndBoundaryVertices[0];
    for (std::pair<int, int> p : vertexMap){
        if (p.second == gluedStartVertexI){
            v0 = mesh.vertex(p.first);
            ctr++;
        }
        if (p.second == gluedEndVertexI){
            v1 = mesh.vertex(p.first);
            ctr++;
        }
        if (ctr == 2) break;
    }
    //set up the gurobi model
    gbModel.setPeriod(period);
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<double> weights; 
    std::tie(vertices, edges, weights) = getVerticesAndEdgesInShortestEdgePath(geometry, v0, v1);
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints; 
    edgePathConstraints.push_back(std::make_pair(weights, 1.0));
    gbModel.setEdgePathConstraints(edgePathConstraints);
    psMesh.addEdgeScalarQuantity("bdy-bdy path", weights);

    //convert it back to row major format 
    dOne = d_oneColMajor;
    int numPairs = 0;
    int maxSingularityPairs = 2;
    //convert omega to an Eigen::VectorXd
    Eigen::VectorXd omegaEig(mesh.nEdges());
    for (Edge e : mesh.edges()){
        omegaEig[e.getIndex()] = omega[e];
    }
    while(numPairs < maxSingularityPairs){
        showd1Omega(geometry, psMesh, omegaEig, numPairs);
        std::pair<int, int> singPair; 
        singPair = findSingularityPair(geometry, timeFunction, omegaEig);
        usedFaceIndices.push_back(singPair.first);
        usedFaceIndices.push_back(singPair.second);
        singPositions[mesh.face(singPair.first)] = 1.0;
        singPositions[mesh.face(singPair.second)] = -1.0;
        std::vector<int> singularFaces(mesh.nFaces(), 0);
        singularFaces[singPair.first] = 1.0;
        singularFaces[singPair.second] = -1.0;
        gbModel.setFaceIndices(singularFaces);
        gbModel.setBdyEdges(globalBdyConditions.courseBdyEdges);
        omegaEig = computeIterativeOneForm(geometry, psMesh, gbModel, edgeMappingsPairs);
        numPairs++;
    }
    return singPositions;
}

//render d1*omega
void showd1Omega(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Eigen::VectorXd& omega, int numPairs){
    SurfaceMesh& mesh = geometry.mesh;
    Eigen::VectorXd d1_omega = dOne * omega;
    psMesh.addFaceScalarQuantity("d1_times_omega_" + std::to_string(numPairs), d1_omega);
}

//find a pair of singularities given certain \omega values
//these singularities will be on the same isolevel set of the time function
std::pair<int, int> findSingularityPair(VertexPositionGeometry& geometry, VertexData<double>& timeFunction,  Eigen::VectorXd& omega){

    SurfaceMesh& mesh = geometry.mesh;
    Eigen::VectorXd d1_omega = dOne * omega;

    std::pair<int, int> toReturn;
    
    //find face with max absolute curl
    double maxAbsCurl = -DBL_MAX;
    int maxFace = -1;

    for (int i = 0; i < d1_omega.size(); i++){
        if (std::find(usedFaceIndices.begin(), usedFaceIndices.end(), i) != usedFaceIndices.end()) continue;//don't use faces that have already been set
        if (std::fabs(d1_omega(i)) > maxAbsCurl){
            maxAbsCurl = d1_omega(i);
            maxFace = i;
        }
    }

    Face maxCurlFace = mesh.face(maxFace);
    double isoVal = (timeFunction[maxCurlFace.halfedge().vertex()] +  timeFunction[maxCurlFace.halfedge().next().vertex()] 
                    + timeFunction[maxCurlFace.halfedge().next().next().vertex()]) / 3.; 

    Eigen::MatrixXd iV;
    Eigen::MatrixXd iE;
    std::vector<int> f;
    std::tie(iV, iE, f) = getTimeFunctionIsoLine(geometry, timeFunction, isoVal);
    polyscope::registerCurveNetwork("Isoline", iV, iE);

    //find singularity of similar curl but opposite sign
    double minDiff = DBL_MAX;
    int minFace = -1;
    for (int i = 0; i < f.size(); i++){
        if (std::find(usedFaceIndices.begin(), usedFaceIndices.end(), f[i]) != usedFaceIndices.end()) continue;//don't use faces that have already been set
        if ((d1_omega(f[i]) * d1_omega(maxFace)) < 0){//max abs curl and current face differ in sign
            double diffSignCurl = std::fabs(d1_omega(f[i]));
            if (std::fabs(maxAbsCurl - diffSignCurl) < minDiff){
                minDiff = std::fabs(maxAbsCurl - diffSignCurl);
                minFace = f[i];
            }
        }
    }
    if (d1_omega(maxFace) > 0 && d1_omega(minFace) < 0)
        return std::make_pair(maxFace, minFace);
    else
        return std::make_pair(minFace, maxFace);

}

//get the vertices (of a curve network), edges (of a curve network) and faces that a particular isovalue of the time function passes through 
//generate isolines for the time function given a specific isoVal
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, std::vector<int>> getTimeFunctionIsoLine(VertexPositionGeometry& geometry, VertexData<double>& timeFunction, double isoVal){

    SurfaceMesh& mesh = geometry.mesh;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::tie(V, F) = getVertexPositionsandFaceLists(geometry);
    Eigen::MatrixXd S(V.size(), 1);
    for (Vertex v : mesh.vertices()){
        S(v.getIndex(), 0) = timeFunction[v];
    }
    Eigen::MatrixXd vals(1, 1);
    vals(0, 0) = isoVal;
    Eigen::MatrixXd iV;
    Eigen::MatrixXd iE;
    Eigen::VectorXd I;
    igl::isolines(V, F, S, vals, iV, iE, I);

    std::vector<int> passes;
    for(int f = 0;f<F.rows();f++)
    {
        if(( 
        S(F(f,0))<isoVal ||
        S(F(f,1))<isoVal ||
        S(F(f,2))<isoVal)
        &&
        (
        S(F(f,0))>isoVal ||
        S(F(f,1))>isoVal ||
        S(F(f,2))>isoVal))
        {
            passes.push_back(f);
        }
    }

    return std::tie(iV, iE, passes);
}

Eigen::VectorXd computeIterativeOneForm(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, Model& gbModel, std::vector<std::pair<int, int>>& edgeMappingsPairs){
    
    SurfaceMesh& mesh = geometry.mesh; 
    double period = gbModel.getPeriod();
    std::vector<int> faceIndices = gbModel.getFaceIndices();
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();
    Eigen::VectorXd toReturn(mesh.nEdges());
    EdgeData<double> oneForm(mesh);

    std::cout << "size of edge path constraints: " << edgePathConstraints.size() << std::endl;

    //solve the model
    using namespace std;
    try {
        //NOTE: APPARENTLY += is not the fastest
        //should look up the gurobi way of doing this

        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //add variable 1-form variable sigma (per edge)
        vector<GRBVar> sigma;
        
        for (size_t i = 0; i < mesh.nEdges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variabls
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdyEdgeIndex : bdyEdges){
            model.addConstr(sigma[bdyEdgeIndex] == 0.0, "Boundary Constraint");
        }

        //second constraint - (d1*sigma) == +/- 1 based on curl
        //d1*sigma == 0 on every other face
        for (int r = 0; r < dOne.outerSize(); ++r ){
            GRBLinExpr lhs = 0;
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(dOne, r); it; ++it ) {
                    lhs += it.value() * sigma[it.col()];
                }
            model.addConstr(lhs == period * faceIndices[r]);
        }

        //third constraint - sigma across stiched edges should be equal 
        for (int i = 0; i < edgeMappingsPairs.size(); i++){
            int iEdge1 = edgeMappingsPairs[i].first;
            int iEdge2 = edgeMappingsPairs[i].second;
            model.addConstr(sigma[iEdge2] == sigma[iEdge1], "Stitched Edge Constraint");
        }

        //fourth constraint - boundary-boundary integral is 1 
        for (int i = 0; i < edgePathConstraints.size(); i++){
            GRBLinExpr pathIntegral = 0;
            for (int j = 0; j < mesh.nEdges(); j++){
                pathIntegral += edgePathConstraints[i].first[j] * sigma[j];
            }
            model.addConstr(pathIntegral == edgePathConstraints[i].second);
        }

        //trying to match energies
        GRBQuadExpr diff1_sum = 0;

        for (int i = 0; i < mesh.nEdges(); i++){
            //diff1_sum += ((sigma[i] - d0_f_avg_gb[i]) * (sigma[i] - d0_f_avg_gb[i]));
            
            //trying the dirichlet stuff
            //diff1_sum +=  geometry.edgeCotanWeight(mesh.edge(i)) * ((sigma[i] * sigma[i]));
            diff1_sum += ((sigma[i] * sigma[i]));
        }
        
        //GRBQuadExpr obj = diff1_sum + diff2_sum;

        //trying the dirichlet stuff 
        GRBQuadExpr obj = diff1_sum;

        model.setObjective(obj);
        model.optimize();

        //put the computed one-form into an edge vector
        for (Edge e : mesh.edges()){
            oneForm[e] = sigma[e.getIndex()].get(GRB_DoubleAttr_X);
        }

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        cout << "-----------------------------------------------" << endl;
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    FaceData<Vector3> gradients = computeOneFormGrad(geometry, oneForm);
    psMesh.addFaceVectorQuantity("Iterative gradients", gradients); 
    EdgeData<double> matchingEnergy = computeMatchingOneForm(geometry, psMesh, 0, gradients, edgeMappingsPairs);
    
    for (Edge e : mesh.edges()){
        toReturn(e.getIndex()) = matchingEnergy[e];
    }
    return toReturn;
}

//compute the gradient (per face) of a 1-form
FaceData<Vector3> computeOneFormGrad(VertexPositionGeometry& geometry, EdgeData<double>& oneForm){
    
    SurfaceMesh& mesh = geometry.mesh;
    //solving per face gradient
    //FACE COUNT DOES NOT INCLUDE BOUNDARY FACES
    Eigen::MatrixXd gradient(mesh.nFaces(), 3);
    
    //re-name to sigma so I can just re-use some code. Ew 
    Eigen::VectorXd sigma(mesh.nEdges());
    for (Edge e : mesh.edges()){
        sigma(e.getIndex()) = oneForm[e];
    }

    FaceData<Vector3> face_gradients(mesh);

    Eigen::MatrixXd faceSystem(3, 3);
    Eigen::VectorXd faceSigmas(3);
    for (Face f : mesh.faces()){
        Halfedge he = f.halfedge();
      	Vector3 face_normal = geometry.faceNormal(f);
      	int j = 0;
      	do{	
      	 	Edge edgeObject = he.edge();
      	    int edgeIndex = edgeObject.getIndex();
      	    if (j == 2){//last row of RHS
      	     	faceSigmas(j) = 0;
      	    }
      	    else{
      	    	faceSigmas(j) = sigma(edgeIndex);
      	    }
      	    for (int i = 0; i < 3; i++){
      	    	Vector3 he_vector = geometry.vertexPositions[edgeObject.halfedge().tipVertex()] 
                                        - geometry.vertexPositions[edgeObject.halfedge().tailVertex()];
      	    	if (j == 2){//last row of LHS
      	    		faceSystem(j, i) = face_normal[i];
      	    	}
      	    	else{
      	    		faceSystem(j, i) = he_vector[i];
      	    	}
      	    }
      	    j++;
      	    he = he.next();
      	}while(he != f.halfedge());
       	Eigen::VectorXd soln = faceSystem.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(faceSigmas);
       	gradient(f.getIndex(), 0) = soln(0);
       	gradient(f.getIndex(), 1) = soln(1);
       	gradient(f.getIndex(), 2) = soln(2);    
        Vector3 face_grad{soln(0), soln(1), soln(2)};
        face_gradients[f] = face_grad;
    }

    return face_gradients;
}
//-------------------End of first strategy---------------------//


//------------------Second strategy----------------------------//
FaceData<int> getGreedySingularityPositions(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, Model& model, 
                                                VertexData<double>& gluedTimeFunction, std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, std::vector<bool>& orientations){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh; 
    SurfaceMesh& globalMesh = globalGeometry.mesh;

    FaceData<int> singularityPositions(globalMesh);
    gluedGeometry.requireDECOperators();
    dOne = gluedGeometry.d1;
    globalToGluedVertexMap = vertexMap;
    //unset the integrability constraint in the model 
    gbModel = model;
    gbModel.setIntegrabilityConstraint(false);
    int maxSingularityPairs = 0; 
    int numPairs = 0; 
    std::vector<std::pair<int, int>> singularFaces;
    //vizualize the non-integrability
    EdgeData<double> sigmaTilde = computeOneForm(globalGeometry, gluedGeometry, gbModel, vertexMap, edgeMap, psMesh);
    EdgeData<double> sigmaTildeGlobal = convertGluedToGlobalEdgeFunction(globalGeometry, gluedGeometry, sigmaTilde, edgeMap);
    Eigen::Map<Eigen::VectorXd> sigmaTildeEig(sigmaTilde.raw().data(), gluedMesh.nEdges());
    Eigen::VectorXd d1Sigma = dOne * sigmaTildeEig; 
    psMesh.addFaceScalarQuantity("d1(sigma" + std::to_string(numPairs) + ")", d1Sigma);
    psMesh.addOneFormTangentVectorQuantity("sigma" + std::to_string(numPairs) + "(Whitney)", sigmaTildeGlobal, orientations);
    while(numPairs < maxSingularityPairs){
        std::pair<int, int> p = findSingularityPair(globalGeometry, gluedGeometry, gluedTimeFunction, sigmaTildeEig);
        usedFaceIndices.push_back(p.first);
        usedFaceIndices.push_back(p.second);   
        singularFaces.push_back(std::make_pair(p.first, 1));
        singularFaces.push_back(std::make_pair(p.second, -1));
        gbModel.setSingularFaceIndices(singularFaces);
        singularityPositions[globalMesh.face(p.first)] = 1; 
        singularityPositions[globalMesh.face(p.second)] = -1;
        numPairs++;
        sigmaTilde = computeOneForm(globalGeometry, gluedGeometry, gbModel, vertexMap, edgeMap, psMesh);
        Eigen::Map<Eigen::VectorXd> sigmaTildeEig(sigmaTilde.raw().data(), gluedMesh.nEdges());
        d1Sigma = dOne * sigmaTildeEig; 
        double sum = 0;
        for (int i = 0; i < d1Sigma.rows(); i++){
            sum += d1Sigma(i);
        }
        std::cout << "Sum of d1 sigma = " << sum << std::endl;
        psMesh.addFaceScalarQuantity("d1 sigma " + std::to_string(numPairs), d1Sigma);
    }
    return singularityPositions;
}

std::pair<int, int> findSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction, const Eigen::VectorXd& sigmaTilde){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    Eigen::VectorXd d1Sigma = dOne * sigmaTilde; 
  
    //find face with max absolute curl
    double maxAbsCurl = -DBL_MAX;
    int maxFace = -1;

    for (int i = 0; i < gluedMesh.nFaces(); i++){
        if (std::find(usedFaceIndices.begin(), usedFaceIndices.end(), i) != usedFaceIndices.end()) continue;//don't use faces that have already been set
        if (std::fabs(d1Sigma(i)) > maxAbsCurl){
            maxAbsCurl = d1Sigma(i);
            maxFace = i;
        }
    }
    Face maxCurlFace = gluedMesh.face(maxFace);
    double isoVal = (gluedTimeFunction[maxCurlFace.halfedge().vertex()] +  gluedTimeFunction[maxCurlFace.halfedge().next().vertex()] 
                    + gluedTimeFunction[maxCurlFace.halfedge().next().next().vertex()]) / 3.;

    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);

    Eigen::MatrixXd iV;
    Eigen::MatrixXd iE;
    std::vector<int> f;
    std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
    //auto isoLine = polyscope::registerCurveNetwork("Isoline", iV, iE);
    //isoLine->setRadius(0.001);

    //find singularity of similar curl but opposite sign
    double minDiff = DBL_MAX;
    int minFace = -1;
    for (int i = 0; i < f.size(); i++){
        if (std::find(usedFaceIndices.begin(), usedFaceIndices.end(), f[i]) != usedFaceIndices.end()) continue;//don't use faces that have already been set
        if ((d1Sigma(f[i]) * d1Sigma(maxFace)) < 0){//max abs curl and current face differ in sign
            double diffSignCurl = std::fabs(d1Sigma(f[i]));
            if (std::fabs(maxAbsCurl - diffSignCurl) < minDiff){
                minDiff = std::fabs(maxAbsCurl - diffSignCurl);
                minFace = f[i];
            }
        }
    }
    if (d1Sigma(maxFace) > 0 && d1Sigma(minFace) < 0){
        return std::make_pair(maxFace, minFace);
    }
    else{
        return std::make_pair(minFace, maxFace);
    }
}