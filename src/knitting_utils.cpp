#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//compute time function using a vector of pairs of vertex mappings instead of the map because we miss stitches then 
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::vector<std::pair<int,int>>& vertexMappingsPairs,  globalBoundaryConditions& boundaryConditions, std::map<int, int>& indexMap){

    SurfaceMesh& mesh = geometry.mesh;
    int numVertices = mesh.nVertices(); 
    int numStitches = vertexMappingsPairs.size();
    //require the cotan edge weights 
    geometry.requireEdgeCotanWeights();
    // //store mappings between index in the original mesh and index in the Laplacian matrix
    // std::map<int, int> indexMap;
    //store mappings between index in the Laplacian matrix to index in the orignal mesh 
    std::map<int, int> laplacianMatrixIndexToOriginalIndex;
    //all the pairs that have been seen so far
    std::vector<std::pair<int, int>> seenPairs;
    //number of "unique vertices" i.e., only consider one vertex per stitch
    int numUniqueVertices = 0;
   //count the number of unique entries in the map
   std::set<int> uniqueEntry; 
   for (auto entry : indexMap){
        uniqueEntry.insert(entry.second);
    }
    numUniqueVertices = uniqueEntry.size();
    Eigen::SparseMatrix<double> L(numUniqueVertices, numUniqueVertices);
    std::vector<Eigen::Triplet<double>> tripletList;
    //keep a set of indices you've already populated in the Laplacian 
    std::set<int> setIndices;//set of vertex indices we've already set in the laplacian

    for (Vertex v : mesh.vertices()){
        if (std::find(setIndices.begin(), setIndices.end(), indexMap.at(v.getIndex())) != setIndices.end()) continue;//we've handled this already
        double L_diag = 0.0;
        //iterate over the 1-ring of the vertex 
        //iterate over the one-ring of the vertex 
        for (Halfedge he : v.outgoingHalfedges()){
            //off diagonal entries 
            tripletList.emplace_back(indexMap.at(v.getIndex()), indexMap.at(he.tipVertex().getIndex()),
            -geometry.edgeCotanWeights[he.edge()]);
            setIndices.insert(indexMap.at(v.getIndex()));

            for (auto p : vertexMappingsPairs){
                if (p.first == v.getIndex()){//grab the contributes from the second in the pair
                    Vertex mappedVertex = mesh.vertex(p.second);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        tripletList.emplace_back(indexMap.at(v.getIndex()), 
                                        indexMap.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                        setIndices.insert(indexMap.at(v.getIndex()));
                    }
                }
                if (p.second == v.getIndex()){
                    Vertex mappedVertex = mesh.vertex(p.first);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        tripletList.emplace_back(indexMap.at(v.getIndex()), 
                                        indexMap.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                        setIndices.insert(indexMap.at(v.getIndex()));
                    }
                }
            }
            //handle the diagonal entries 
            //handle the diagonal entry case
            L_diag += geometry.edgeCotanWeights[he.edge()];
            for (auto p : vertexMappingsPairs){
                if (p.first == v.getIndex()){
                    Vertex mappedVertex = mesh.vertex(p.second);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        L_diag += geometry.edgeCotanWeights[mappedVertexHalfedge.edge()];
                    }
                }
                if (p.second == v.getIndex()){
                    Vertex mappedVertex = mesh.vertex(p.first);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        L_diag += geometry.edgeCotanWeights[mappedVertexHalfedge.edge()];
                    }
                }
            }
        }
        tripletList.emplace_back(indexMap.at(v.getIndex()), indexMap.at(v.getIndex()), L_diag);
        setIndices.insert(indexMap.at(v.getIndex()));
    }

    L.setFromTriplets(tripletList.begin(), tripletList.end());
    //force boundary conditions
    Eigen::VectorXd b = Eigen::VectorXd::Zero(numUniqueVertices);
    
    for (int v : boundaryConditions.courseStartBoundaryVertices){
        int updatedIndex = indexMap.at(v);
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 0.0;
    }

    for (int v : boundaryConditions.courseEndBoundaryVertices){
        int updatedIndex = indexMap.at(v);
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 1.0;
    }
    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed" << std::endl;
    }
    Eigen::VectorXd u = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed" << std::endl;
    }
    VertexData<double> timeFunction(mesh);
    for (Vertex v : mesh.vertices()){
        timeFunction[v] = u(indexMap.at(v.getIndex()));
    }
    return timeFunction;
}

//compute time function directly in the glued mesh setting 
VertexData<double> computeTimeFunction(EdgeLengthGeometry& gluedGeometry, globalBoundaryConditions& bdyConditions){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh; 
    VertexData<double> gluedTimeFunction(gluedMesh);
    gluedGeometry.requireCotanLaplacian();
    Eigen::SparseMatrix<double> L = gluedGeometry.cotanLaplacian;
    //force boundary conditions 
    Eigen::VectorXd b = Eigen::VectorXd::Zero(gluedMesh.nVertices());
    for (int v : bdyConditions.courseStartBoundaryVertices){
        L.row(v) *= 0.0;
        L.coeffRef(v, v) = 1.0;
        b(v) = 0;
    }
    for (int v : bdyConditions.courseEndBoundaryVertices){
        L.row(v) *= 0.0;
        L.coeffRef(v, v) = 1.0;
        b(v) = 1.0;
    }
    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed" << std::endl;
    }
    Eigen::VectorXd u = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed" << std::endl;
    }
    for (Vertex v : gluedMesh.vertices()){
        gluedTimeFunction[v] = u(v.getIndex());
    }
    return gluedTimeFunction;
}

//compute the gradient of a function defined as a scalar over vertices in the global mesh setting 
FaceData<Vector3> computeTimeFunctionFaceGrad(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction){
    
    SurfaceMesh& mesh = geometry.mesh;
    FaceData<Vector3> faceGradients(mesh);

    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    std::tie(V, F) = getVertexPositionsandFaceLists(geometry);

    Eigen::VectorXd U(mesh.nVertices());
    //convert time function to eigen matrix 
    for (Vertex v : mesh.vertices()){
        U(v.getIndex()) = vertexScalarFunction[v];
    }

    // Compute gradient operator: #F*3 by #V
    SparseMatrix<double> G;
    igl::grad(V,F,G);

    // Compute gradient of U
    //gradients stored per faces as (3F * 1) 
    //all xs, then all ys, all zs
    Eigen::MatrixXd GU = (G*U);

    for (int i = 0; i < GU.rows(); i+=3){
        int currFaceI = i / 3.;
        faceGradients[mesh.face(currFaceI)] = Vector3{GU(currFaceI), GU(currFaceI + F.rows()), GU(currFaceI + 2 * F.rows())};
    }
    return faceGradients;
}

//compute the gradient of the a function defined as a scalar over vertices in the glued mesh setting 
//don't think the formula I'm using in here is right
FaceData<Vector3> computeTimeFunctionFaceGrad(EdgeLengthGeometry& geometry, VertexData<double>& vertexScalarFunction){

    SurfaceMesh& mesh = geometry.mesh; 
    geometry.requireFaceAreas();
    geometry.requireEdgeCotanWeights();
    geometry.requireCornerAngles();
    FaceData<Vector3> gradients(mesh);

    //don't know if this is the right thing to do here
    for (Face f : mesh.faces()){
        double fi = vertexScalarFunction[f.halfedge().vertex()];
        double fj = vertexScalarFunction[f.halfedge().next().vertex()];
        double fk = vertexScalarFunction[f.halfedge().next().next().vertex()];
        double area = geometry.faceAreas[f];
        BarycentricVector X_ik_perp = (BarycentricVector(f, Vector3{1., 0, -1.})).rotated90(geometry);
        BarycentricVector X_ji_perp = (BarycentricVector(f, Vector3{-1., 1., 0.})).rotated90(geometry);
        BarycentricVector gradF = ((fj - fi) * X_ik_perp + (fk - fi) * X_ji_perp)/ (2. * area); 
        BarycentricVector gradFNormalized = gradF / norm(geometry, gradF);
        gradients[f] = Vector3{gradFNormalized.faceCoords[0], gradFNormalized.faceCoords[1], gradFNormalized.faceCoords[2]};
    }
    return gradients;
    

    // for (Face f : mesh.faces()){
    //     double f1 = vertexScalarFunction[f.halfedge().vertex()];
    //     double f2 = vertexScalarFunction[f.halfedge().next().vertex()];
    //     double f3 = vertexScalarFunction[f.halfedge().next().next().vertex()];
    //     Corner c1 = f.halfedge().corner();
    //     Corner c2 = f.halfedge().next().corner();
    //     Corner c3 = f.halfedge().next().next().corner();
    //     double cottheta1 = 1. / tan(geometry.cornerAngles[c1]);
    //     double cottheta2 = 1. / tan(geometry.cornerAngles[c2]);
    //     double cottheta3 = 1. / tan(geometry.cornerAngles[c3]);
    //     double area = geometry.faceAreas[f];
    //     //each of these components are mutliplied by unit vectors 
    //     //each unit vector is perpendicular to the opposite edge of the vertex
    //     gradients[f][0] = 1./(2.*area) * (((f2 - f1) * cottheta3) + ((f3 - f1) * cottheta2));
    //     gradients[f][1] = 1./(2.*area) * (((f3 - f2) * cottheta1) + ((f1 - f2) * cottheta3));
    //     gradients[f][2] = 1./(2.*area) * (((f1 - f3) * cottheta2) + ((f2 - f3) * cottheta1));
    // }
    // return gradients;
}

//compute a vector (in ambient space) per vertex that is aligned with the gradient of a scalar field
VertexData<Vector3> computeVertexValuedField(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction, double angle){

    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireVertexNormals();

    VertexData<Vector3> vertexValuedField(mesh);

    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    std::tie(V, F) = getVertexPositionsandFaceLists(geometry);

    Eigen::VectorXd U(mesh.nVertices());
    //convert time function to eigen matrix 
    for (Vertex v : mesh.vertices()){
        U(v.getIndex()) = vertexScalarFunction[v];
    }

    // Compute gradient operator: #F*3 by #V
    SparseMatrix<double> G;
    igl::grad(V,F,G);

    // Compute gradient of U
    Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G*U).eval().data(),F.rows(),3);

    Eigen::MatrixXd S;
    igl::average_onto_vertices(V, F, GU, S);

    //store as vector on each vertex 
    for (Vertex v : mesh.vertices()){
        vertexValuedField[v] = Vector3{S(v.getIndex(), 0), S(v.getIndex(), 1), S(v.getIndex(), 2)};
        //rotate by specified angle 
        vertexValuedField[v] = vertexValuedField[v].rotateAround(geometry.vertexNormals[v], angle);
    }

    return vertexValuedField;

}

//get a line field per vertex from a vertex valued vector field in ambient space
VertexData<Vector2> vertexDirectionField(VertexPositionGeometry& geometry, VertexData<Vector3>& vertexValuedField){

    SurfaceMesh& mesh = geometry.mesh;
    VertexData<Vector2> directionField(mesh);

    geometry.requireVertexNormals();

    //convert vertexValuedField to eigen matrix
    Eigen::MatrixXd vertex_valued_field(mesh.nVertices(), 3);
    for (Vertex v : mesh.vertices()){
        vertex_valued_field(v.getIndex(), 0) = vertexValuedField[v][0];
        vertex_valued_field(v.getIndex(), 1) = vertexValuedField[v][1];
        vertex_valued_field(v.getIndex(), 2) = vertexValuedField[v][2];

    }

    //angular coordinate of every halfedge
    HalfedgeData<double> angularCoordinate(mesh);

    //first compute the angle of each halfedge 
    HalfedgeData<double> angle(mesh);
  
    for (Halfedge he : mesh.halfedges()){
        Vector3 a = geometry.vertexPositions[he.next().next().tailVertex()];
        Vector3 b = geometry.vertexPositions[he.tailVertex()];
        Vector3 c = geometry.vertexPositions[he.next().tailVertex()];
        Vector3 u = (b - a).unit();
        Vector3 v = (c - a).unit();
        angle[he] = std::acos(std::max(-1.0, std::min(1.0, dot(u, v))));
    }

    //compute angular coordinate at each outgoing halfedge 
    for (Vertex v : mesh.vertices()){
        //compute the cumulative angle at each outgoing
        //halfedge, relative to the intitial halfedge
        double cumulativeAngle = 0.0;
        Halfedge he = v.halfedge();
        if (!v.isBoundary()){
            do{
                angularCoordinate[he] = cumulativeAngle;
                cumulativeAngle += angle[he.next()];
                he = he.next().next().twin();
            }
            while(he != v.halfedge());
            do{
                angularCoordinate[he] *= 2.0*PI / cumulativeAngle;
                he = he.twin().next();
            }
            while(he != v.halfedge());
        }
        else{
            //ensure the halfedge starts "on boundary" (twin is member of boundary cycle)
            while(!he.twin().face().isBoundaryLoop()){
                he = he.twin().next();
            }
            //get the halfedge on the boundary
            he = he.twin().next();
            Halfedge start = he;
            //calculate the cumulative angle
            while(!he.face().isBoundaryLoop()){
                angularCoordinate[he] = cumulativeAngle;
                cumulativeAngle += angle[he.next()];
                he = he.next().next().twin();
            }
            angularCoordinate[he] = cumulativeAngle;//assign angular coordinate for last halfedge
            //reset he
            he = start;
            //shift angle so reference halfedge is 0
            //loop over the one-ring of the boundary vertex until you find a halfedge that's part of the boundary cycle
            Halfedge bdy_he = v.halfedge();
            while(bdy_he.isInterior()){
                bdy_he = bdy_he.next().next().twin();
            }
            double shift = angularCoordinate[v.halfedge()];
            while(!he.face().isBoundaryLoop()){
                angularCoordinate[he] -= shift;
                he = he.next().next().twin();
            }
            angularCoordinate[he] -= shift;//assign angular coordinate for last halfedge
        }
    }

    for (Vertex v : mesh.vertices()){
        double alpha = angularCoordinate[v.halfedge()]; 
        std::complex<double> r(cos (2.0 * alpha), sin(2.0 * alpha));
        int i = v.getIndex();
        Vector3 n = geometry.vertexNormals[v];
        Vector3 e = geometry.vertexPositions[v.halfedge().tipVertex()] - geometry.vertexPositions[v.halfedge().tailVertex()];
        Vector2 u = projectOntoPlane(vertex_valued_field.row(i).transpose(), {n.x, n.y, n.z}, {e.x, e.y, e.z});
        double a = std::atan2(u.y, u.x);
        //for a 2-direction field
        std::complex<double> complexDirectionField(r * std::complex<double>(cos(2.0 * a), sin(2.0 * a)));
        directionField[v] = Vector2::fromComplex(complexDirectionField);
        directionField[v] = unit(directionField[v]);

        //directionField[v] = rotate90 ? directionField[v].rotate90() : directionField[v];
    }

    //way to convert a 1-direction field vector to a 2-direction field vector as defined by geometry central 
    // 1. rotate the vector 90 degrees then square
    // 2. square the vector then rotate by 180 degrees
    

    return directionField;

}

//compute a 1-form over the mesh given an input optimization problem 
EdgeData<double> computeOneForm(VertexPositionGeometry& geometry, Model& gbModel, polyscope::SurfaceMesh& globalPSMesh){

    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireDECOperators();

    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one;
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_zero;
    d_zero = geometry.d0;
    d_one = geometry.d1;
    
    //convert d_one to column major for easy updates of columns 
    Eigen::SparseMatrix<double, Eigen::ColMajor> d_oneColMajor;
    d_oneColMajor = d_one;

    EdgeData<double> oneForm(mesh);
    FaceData<double> integratedOneForm(mesh);
    EdgeData<double> difference(mesh);

    std::vector<double> omega = gbModel.getMatchingTerms();
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<int, int>> edgeMappingsPairs = gbModel.getEdgeMappingsPairs();
    double period = gbModel.getPeriod();
    std::vector<std::vector<double>> waleBdyPathConstraints = gbModel.getWaleBdyPathConstraints();

    //invert the signs to account for "stitched together" panels in d0
    for (std::pair<int, int> p : edgeMappingsPairs){
        int iE2 = p.second;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d_zero, iE2); it; ++it){
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
    //convert it back to row major format 
    d_one = d_oneColMajor;
    Eigen::VectorXd constFunc(mesh.nVertices());
    constFunc.setOnes();
    Eigen::SparseMatrix<double> prod = d_one * d_zero;
    
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        //add variable 1-form variable sigma (per edge)
        std::vector<GRBVar> sigma;

        for (size_t i = 0; i < mesh.nEdges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //add variable for integral of 1-form over each face
        std::vector<GRBVar> k;
        //add integer variable k (per face) (not including faces represented as boundaries)
        for (size_t i = 0; i < mesh.nFaces(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            k.push_back(k_i);//decision variables
        }

        //add boundary integral variable for wale direction stripes 
        std::vector<GRBVar> waleBdyIntegerConstraints;
        for (size_t i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            waleBdyIntegerConstraints.push_back(k_i);
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[bdy_edge_index] == 0.0, "Boundary Constraint");
        }

        //second constraint - (d1*sigma) == period*k
        for (int r = 0; r < d_one.outerSize(); ++r) {
            GRBLinExpr lhs = 0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d_one, r); it; ++it ) {
                lhs += it.value() * sigma[it.col()];
            }
            model.addConstr(lhs == 0, "Integral Constraint");
            //model.addConstr(lhs == period * k[r], "Integral Constraint");
        }

        //third constraint - sigma across stiched edges should be equal 
        for (int i = 0; i < edgeMappingsPairs.size(); i++){
            int iEdge1 = edgeMappingsPairs[i].first;
            int iEdge2 = edgeMappingsPairs[i].second;
            model.addConstr(sigma[iEdge2] == sigma[iEdge1], "Stitched Edge Constraint");
        }

        //fourth contraint - boundary integral in the wale direction
        for (int i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBLinExpr pathIntegral = 0;
            for (int j = 0; j < mesh.nEdges(); j++){
                pathIntegral += waleBdyPathConstraints[i][j] * sigma[j];
            }
            model.addConstr(pathIntegral == period * waleBdyIntegerConstraints[i]);
        }

        //trying to match energies
        GRBQuadExpr diff1_sum = 0;

        for (int i = 0; i < mesh.nEdges(); i++){
            diff1_sum += (sigma[i] - omega[i]) * (sigma[i] - omega[i]);
        }


        GRBQuadExpr obj = diff1_sum;

        model.setObjective(obj);
        model.optimize(); 
        std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    
        //put the computed one-form into an edge vector
        for (Edge e : mesh.edges()){
            oneForm[e] = sigma[e.getIndex()].get(GRB_DoubleAttr_X);
        }
        //put the singular faces into a face vector 
        for (Face f : mesh.faces()){
            integratedOneForm[f] = k[f.getIndex()].get(GRB_DoubleAttr_X);
        }

        //put the difference into a mesh for viz purposes 
        for (Edge e : mesh.edges()){
            difference[e] = (sigma[e.getIndex()].get(GRB_DoubleAttr_X) - omega[e.getIndex()]) * (sigma[e.getIndex()].get(GRB_DoubleAttr_X) - omega[e.getIndex()]);
        }

        globalPSMesh.addEdgeScalarQuantity("difference", difference);
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    return oneForm;
}

EdgeData<double> computeOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& vertexMap, std::map<int, int>& edgeMap,
                                        polyscope::SurfaceMesh& psMesh){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh; 
    SurfaceMesh& globalMesh = globalGeometry.mesh;

    EdgeData<double> oneForm(gluedMesh);
    //std::vector<double> omega = gbModel.getMatchingTerms();
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::array<double, 3>> gradients = gbModel.getFaceGradients();
    double period = gbModel.getPeriod();
    std::vector<double> omega = gbModel.getMatchingTerms();
    std::vector<std::vector<double>> waleBdyPathConstraints = gbModel.getWaleBdyPathConstraints();
    std::vector<std::pair<int, int>> singularFaceIndices = gbModel.getSingularFaces();
    std::vector<int> faceIndices = gbModel.getFaceIndices();
    bool hasIntegrabilityConstraint = gbModel.getIntegrabilityConstraint();
    //debug stuff 
    bool useEdgeAveraging = gbModel.useEdgeAveraging;
    bool useFaceDifferenceViz = gbModel.useFaceDifferenceViz;

    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;

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
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 180);
        //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_MIPFocus, 2);
        
        //add variable 1-form variable sigma (per edge)
        std::vector<GRBVar> sigma;

        for (size_t i = 0; i < gluedMesh.nEdges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //add variable for integral of 1-form over each face
        std::vector<GRBVar> k;
        //add integer variable k (per face) (not including faces represented as boundaries)
        for (size_t i = 0; i < gluedMesh.nFaces(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            k.push_back(k_i);//decision variables
        }

        //add boundary integral variable for wale direction stripes 
        std::vector<GRBVar> waleBdyIntegerConstraints;
        for (size_t i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            waleBdyIntegerConstraints.push_back(k_i);
        }

        //first constraint - sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[bdy_edge_index] == 0.0, "Boundary Constraint");
        }

        //regularization term 
        GRBQuadExpr diff2_sum = 0;
        //compute nP over every face 
        //add the second constraint while we're here
        //second constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (int r = 0; r < gluedMesh.nFaces(); ++r){
            GRBLinExpr nPCurr = 0;
            GRBLinExpr lhs = 0.0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d_one, r); it; ++it) {
                nPCurr += it.value() * sigma[it.col()];
                lhs += it.value() * sigma[it.col()];
            }
            nP[r] = nPCurr;
            //integrability constraint
            if (hasIntegrabilityConstraint){
                if (faceIndices.size() == 0){//no face indices specified, either set it to 0, or do MIP search
                    model.addConstr(lhs == 0, "Integral Constraint");
                    diff2_sum += ((lhs - 0) * (lhs - 0));//regularization for all faces set to 0 is dumb
                    //model.addConstr(lhs == period * k[r], "Integral Constraint");
                    //diff2_sum += ((lhs - period * k[r]) * (lhs - period * k[r]));
                }
                else{
                    model.addConstr(lhs == period * faceIndices[r], "Integral Constraint");
                    diff2_sum += ((lhs - period * faceIndices[r]) * (lhs - period * faceIndices[r]));
                }
            }
        }
        //third constraint - boundary integral in the wale direction
        for (int i = 0; i < waleBdyPathConstraints.size(); i++){
            GRBLinExpr pathIntegral = 0;
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                pathIntegral += waleBdyPathConstraints[i][j] * sigma[j];
            }
            model.addConstr(pathIntegral == period * waleBdyIntegerConstraints[i]);
        }

        //fourth constraint - place singularities at specified faces 
        for (std::pair<int, int> p : singularFaceIndices){
            int fIndex = p.first; 
            int integralVal = p.second;
            GRBLinExpr lhs = 0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d_one, fIndex); it; ++it){
                lhs += it.value() * sigma[it.col()];
            }
            model.addConstr(lhs == period * integralVal, "Specified singularity constraint");
        }

        //trying to match energies
        GRBQuadExpr diff1_sum = 0;
        for (int i = 0; i < omega.size(); i++){
            diff1_sum += (sigma[i] - omega[i]) * (sigma[i] - omega[i]);
        }

        //compute a piecewise linear function over the vertices of the mesh 
        std::vector<GRBLinExpr> u(gluedMesh.nVertices());
        std::vector<GRBLinExpr> uGluedHalfedge(gluedMesh.nHalfedges());
        std::vector<std::vector<GRBLinExpr>> gradU(gluedMesh.nFaces(), std::vector<GRBLinExpr>(3));
        for (Face f : gluedMesh.faces()){
            int signhIJ = f.halfedge().orientation() ? 1 : -1;
            int signhJK = f.halfedge().next().orientation() ? 1 : -1;
            u[f.halfedge().vertex().getIndex()] = 0.0;
            uGluedHalfedge[f.halfedge().getIndex()] = 0.0;
            u[f.halfedge().next().vertex().getIndex()] = (signhIJ * sigma[f.halfedge().edge().getIndex()]) - (nP[f.getIndex()]/3.0);
            uGluedHalfedge[f.halfedge().next().getIndex()] = (signhIJ * sigma[f.halfedge().edge().getIndex()]) - (nP[f.getIndex()]/3.0);
            u[f.halfedge().next().next().vertex().getIndex()] = (signhIJ * sigma[f.halfedge().edge().getIndex()] + signhJK * sigma[f.halfedge().next().edge().getIndex()]) - ((2.0 * nP[f.getIndex()])/3.0);
            uGluedHalfedge[f.halfedge().next().next().getIndex()] = (signhIJ * sigma[f.halfedge().edge().getIndex()] + signhJK * sigma[f.halfedge().next().edge().getIndex()]) - ((2.0 * nP[f.getIndex()])/3.0);
            
            GRBLinExpr currGradU = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + F.rows()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * F.rows()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
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

        GRBQuadExpr obj;
        //which energy we're trying to minimize
        if (useEdgeAveraging){
            obj = (diff1_sum + diff2_sum);
        }
        else{
            obj = (gradDiff + diff2_sum);
        }
        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        //put the computed one-form into an edge vector
        for (Edge e : gluedMesh.edges()){
            oneForm[e] = sigma[e.getIndex()].get(GRB_DoubleAttr_X);
        }
        std::vector<std::array<double, 3>> gradientU(gluedMesh.nFaces());
        std::vector<double> differenceViz(gluedMesh.nFaces());
    
        if (useFaceDifferenceViz){
            for (Face f : gluedMesh.faces()){
                std::array<double, 3> currGrad = {gradU[f.getIndex()][0].getValue(), gradU[f.getIndex()][1].getValue(), gradU[f.getIndex()][2].getValue()};
                gradientU[f.getIndex()] = currGrad;
                differenceViz[f.getIndex()] = difference[f.getIndex()].getValue();
            }
            if (faceIndices.size() == 0){
                psMesh.addFaceScalarQuantity("objective difference (sigma with no sings)", differenceViz);
                psMesh.addFaceVectorQuantity("gradientU from utils (sigma with no sings)", gradientU);
            }
            else{
                psMesh.addFaceScalarQuantity("objective difference (sigma with sings)", differenceViz);
                psMesh.addFaceVectorQuantity("gradientU from (sigma with sings)", gradientU);
            }
        }
        
        for (int i = 0; i < waleBdyIntegerConstraints.size(); i++){
            std::cout << "integer constraint: " << waleBdyIntegerConstraints[i].get(GRB_DoubleAttr_X) << std::endl;
        }
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return oneForm;
}


//compute \omega that is the 1-form we're trying to match over each edge 
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3>& faceGradients){
    
    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireFaceNormals();
    EdgeData<double> d0_f_avg(mesh);


    //if we're computing the matching 1-form in the wale direction, rotate all the gradients 
    if (direction == 1){
        for (Face f : mesh.faces()){
            faceGradients[f] = faceGradients[f].rotateAround(geometry.faceNormals[f], PI/2.);
        }
    }

    //evaluate the matching energy
    for (Edge edge : mesh.edges()){
        //normalize the gradients first
        faceGradients[edge.halfedge().face()] = faceGradients[edge.halfedge().face()].normalize();
        faceGradients[edge.halfedge().twin().face()] = faceGradients[edge.halfedge().twin().face()].normalize();
        //if halfedge is interior
        if (edge.halfedge().twin().isInterior()){
            if (edge.halfedge().orientation()){
                d0_f_avg[edge] = 0.5 * dot((faceGradients[edge.halfedge().face()] + faceGradients[edge.halfedge().twin().face()]),
                                    geometry.vertexPositions[edge.halfedge().tipVertex()] - geometry.vertexPositions[edge.halfedge().tailVertex()]);
            }
        }
      	//handle edges on the boundary
      	else{
      		d0_f_avg[edge] = dot(faceGradients[edge.halfedge().face()], geometry.vertexPositions[edge.halfedge().tipVertex()] - 
                                                                                geometry.vertexPositions[edge.halfedge().tailVertex()]);
      	}
    }

    return d0_f_avg;
}
//compute matching 1-form while taking into account "stitched together" edges
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3>& faceGradients, 
                                    std::vector<std::pair<int, int>>& edgeMappingsPairs){

    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireFaceNormals();
    EdgeData<double> omega(mesh);

    FaceData<Vector3> faceGradientsCopy = faceGradients;

    //if we're computing the matching 1-form in the wale direction, rotate all the gradients 
    if (direction == 1){
        for (Face f : mesh.faces()){
            //first normalize the gradients? 
            faceGradientsCopy[f] = faceGradientsCopy[f].normalize();
            faceGradientsCopy[f] = faceGradientsCopy[f].rotateAround(geometry.faceNormals[f], PI/2.);
        }
    }
    int numStitchedEdges = 0;
    //create a map from the mapped edges
    std::map<int, int> edgeMap;
    for (std::pair<int, int> pair : edgeMappingsPairs){
        edgeMap.insert({pair.first, pair.second});
    }

    //edges which we've handles already
    std::map<int, bool> seenEdges;
    for (Edge e : mesh.edges()){
        seenEdges.insert({e.getIndex(), false});
    }

    for (Edge e : mesh.edges()){
        if (seenEdges[e.getIndex()]) continue;
        if (e.halfedge().twin().isInterior()){//found an interior halfedge
            //normalize the face gradients first
            faceGradientsCopy[e.halfedge().face()] = faceGradientsCopy[e.halfedge().face()].normalize();
            faceGradientsCopy[e.halfedge().twin().face()] = faceGradientsCopy[e.halfedge().twin().face()].normalize();
            Vector3 eVector = geometry.vertexPositions[e.halfedge().tipVertex()] - geometry.vertexPositions[e.halfedge().tailVertex()];
            //eVector = eVector.normalize();
            omega[e] = 0.5 * dot((faceGradientsCopy[e.halfedge().face()] + faceGradientsCopy[e.halfedge().twin().face()]),
                                    eVector);
            seenEdges[e.getIndex()] = true;
        }
        else{//found a boundary halfedge
            if (edgeMap.find(e.getIndex()) != edgeMap.end()){//found a stitched together edge
                numStitchedEdges++;
                Vector3 faceGradientsCopy1 = faceGradientsCopy[e.halfedge().face()].normalize();
                Vector3 faceGradientsCopy2 = faceGradientsCopy[mesh.edge(edgeMap.at(e.getIndex())).halfedge().face()].normalize();
                //take the average direction vector? 
                Vector3 e1 = geometry.vertexPositions[e.halfedge().tipVertex()] - geometry.vertexPositions[e.halfedge().tailVertex()];
                Vector3 e2 = geometry.vertexPositions[mesh.edge(edgeMap.at(e.getIndex())).halfedge().tipVertex()] 
                                            - geometry.vertexPositions[mesh.edge(edgeMap.at(e.getIndex())).halfedge().tailVertex()];
                // Vector3 avgVector = (e1 + e2)/2.;
                //just pick the original edge as the "canonical" direction in the global mesh
                //I'm not really sure the polyscope Whitney interpolation scheme is the best way to visualize these
                //e1 = e1.normalize();
                omega[e] = 0.5 * dot((faceGradientsCopy1 + faceGradientsCopy2), e1);
                omega[mesh.edge(edgeMap.at(e.getIndex()))] = 0.5 * dot((faceGradientsCopy1 + (faceGradientsCopy2)), e1);
                seenEdges[e.getIndex()] = true;
                seenEdges[edgeMap.at(e.getIndex())] = true;
            }
            else{//found a boundary edge that's not stitched to anything
                faceGradientsCopy[e.halfedge().face()] = faceGradientsCopy[e.halfedge().face()].normalize();
                Vector3 eVector = geometry.vertexPositions[e.halfedge().tipVertex()] - geometry.vertexPositions[e.halfedge().tailVertex()];
                //eVector = eVector.normalize();
                omega[e] = dot(faceGradientsCopy[e.halfedge().face()], eVector);
                seenEdges[e.getIndex()] = true;
            }
        }
    }
    return omega;
}


//compute a edge-based field through the optimization
//here the singularities are placed on the vertices
HalfedgeData<double> computeVertexSingularityField(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, 
                                                    polyscope::SurfaceMesh& psMesh, std::map<int, int>& vertexMap,
                                                    std::map<int, std::vector<Halfedge>>& gluedOneRingMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> oneForm(gluedMesh);

    //query information from the model 
    std::vector<int> bdyEdges = model.getBdyEdges();
    std::vector<std::vector<double>> waleBdyPathConstraints = model.getWaleBdyPathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    std::vector<std::pair<int, int>> singularFaceIndices = model.getSingularFaces();
    std::vector<int> faceIndices = model.getFaceIndices();
    double period = model.getPeriod();
    bool hasIntegrabilityConstraint = model.getIntegrabilityConstraint();
    std::vector<std::pair<int, int>> singularEdges = model.getSingularEdges();
    std::vector<std::array<double, 3>> gradients = model.getFaceGradients();

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
            if (handled[he.edge()]) continue; //if the halfedge is already handled
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
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + F.rows()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * F.rows()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
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
        std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            oneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }

        FaceData<Vector3> gradientU(gluedMesh);
        VertexData<double> curl(globalMesh);
        for (Face f : gluedMesh.faces()){
            std::array<double, 3> currGrad = {gradU[f.getIndex()][0].getValue(), gradU[f.getIndex()][1].getValue(), gradU[f.getIndex()][2].getValue()};
            gradientU[f.getIndex()] = Vector3{currGrad[0], currGrad[1], currGrad[2]};
        }
        for (Vertex vi : globalMesh.vertices()){
            double sum = 0.0;
            for (Halfedge he : gluedOneRingMap[vi.getIndex()]){
                Halfedge hjk = he.next();
                if (!hjk.isInterior()) continue;
                Vector3 hjkVec = globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()];
                sum += dot(hjkVec, gradientU[globalMesh.face(he.face().getIndex())]);
            }
            curl[vi] = sum;
        }
        
        //psMesh.addFaceVectorQuantity("gradientU vertex singularity", gradientU);
        //psMesh.addVertexScalarQuantity("vertex curl using gradientU", curl);

    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return oneForm;
}

//@clean
//compute the per-face gradient of a 1-form (in global setting)
//the gradeints ARE NOT NORMALIZED
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

//@clean
//compute the per-face gradient of a 1-form (in global setting)
//the gradients are not normalized 
std::vector<std::array<double, 3>> vectorOneFormFaceGrad(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, HalfedgeData<double>& sigmaTilde){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    //FaceData<Vector3> gradients(globalMesh);
    std::vector<std::array<double, 3>> gradients;

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
        //gradients[f] = Vector3{soln(0), soln(1), soln(2)};
        gradients.push_back(std::array<double, 3>{soln(0), soln(1), soln(2)});
    }
    return gradients;
}



//@clean 
//compute a wale one-form s
HalfedgeData<double> computeWaleOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, 
                                        Eigen::SparseMatrix<double, Eigen::RowMajor>& G, std::map<int, int>& vertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    double period = model.getPeriod();
    std::vector<double> omega = model.getMatchingTerms();
    std::vector<int> faceIndices = model.getFaceIndices();
    std::vector<std::pair<int , int>> singularEdges = model.getSingularEdges();
    //print out the singular edges 
    std::vector<std::vector<double>> waleBdyPathConstraints = model.getWaleBdyPathConstraints();
    std::vector<std::array<double, 3>> gradients = model.getFaceGradients();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    HalfedgeData<double> oneForm(gluedMesh);

    std::cout << "Solving model for wale stripes..." << std::endl;

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        //don't log output to console 
        env.set(GRB_IntParam_OutputFlag, 0);
        env.set("LogFile", "1-form computation.log");
        env.start();
        // Create an empty model
        GRBModel model = GRBModel(env);

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

        //first constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //also compute nP 
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            nP[f.getIndex()] = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            //model.addConstr(lhs == period * faceIndices[f.getIndex()]);
            model.addConstr(lhs == 0);
        }

        //second constraint 
        //ensure that values are opposite sign across halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            //invert the sign of the constraint
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] ==  p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //third constraint 
        //add constraints in the wale direction
        for (int i = 0; i < waleBdyPathConstraints.size(); i++){
            std::vector<double> path = waleBdyPathConstraints[i];
            GRBLinExpr pathIntegral = 0;
            std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
            for (int j = 0; j < gluedMesh.nEdges(); j++){
                if (path[j] > 0){
                    hePath[gluedMesh.edge(j).halfedge().getIndex()] = path[j];
                }
                else if (path[j] < 0){
                    hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = path[j];
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
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
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

        GRBQuadExpr edgeDiff = 0;
        // for (Halfedge he : gluedMesh.halfedges()){
        //     double matchingTerm = he.orientation() ? omega[he.edge().getIndex()] : -1.0 * omega[he.edge().getIndex()];
        //     edgeDiff += (sigma[he.getIndex()] - matchingTerm) * (sigma[he.getIndex()] - matchingTerm);
        // }

        GRBQuadExpr obj;
        obj = gradDiff;

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        
        //put the computed one-form into an halfedge vector
        for (Halfedge he : gluedMesh.halfedges()){
            oneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return oneForm;

}


//@clean 
std::tuple<CornerData<double>, EdgeData<double>> computeWaleStripeInfo(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    std::vector<std::pair<int, int>>& edgeMappingsPairs, std::map<int, int>& edgeMap, 
                                                                    std::map<int, int>& vertexMap, VertexData<double>& timeFunctionGlobal, FaceData<Vector3>& courseOneFormGrad, 
                                                                    Eigen::SparseMatrix<double, Eigen::RowMajor>& G, double period, double knoppelFrequency, globalBoundaryConditions& globalBdyConditions,
                                                                    EdgeData<double>& courseSingularEdgesGlobal, polyscope::SurfaceMesh& psMesh){

    //compute a line field in the tangent space of the vertex
    VertexData<Vector3> vertexVectorField = computeVertexValuedField(globalGeometry, timeFunctionGlobal, PI/2.);
    VertexData<Vector2> lineField = vertexDirectionField(globalGeometry, vertexVectorField);
    EdgeData<double> waleSingularEdgesGlobal(globalGeometry.mesh, 0);
    VertexData<double> freq(globalGeometry.mesh, 1./(period));
    CornerData<double> stripeValues(globalGeometry.mesh);
    FaceData<int> stripeSingularities(globalGeometry.mesh);
    FaceData<int> fieldSingularities(globalGeometry.mesh);
    std::tie(stripeValues, stripeSingularities, fieldSingularities) = computeStripePattern(globalGeometry, freq, lineField); // this is a GC call
    // Do some visualization
    psMesh.addVertexVectorQuantity("vertexVectorField", vertexVectorField);
    psMesh.addFaceScalarQuantity("knoppel face singularities", stripeSingularities);
    psMesh.addFaceScalarQuantity("knoppel field singularities", fieldSingularities);
    
    std::vector<Vector3> knoppelPos; 
    std::vector<std::array<size_t, 2>> knoppelEdges; 
    std::tie(knoppelPos, knoppelEdges) = extractPolylinesFromStripePattern(globalGeometry, stripeValues, stripeSingularities,
                                            fieldSingularities, lineField, false);
    auto knoppelStripes = polyscope::registerCurveNetwork("knoppel wale stripes stripes", knoppelPos, knoppelEdges);
    knoppelStripes -> setRadius(0.001);
    knoppelStripes -> setEnabled(false);

    EdgeData<double> omegaWaleGlobal = computeMatchingOneForm(globalGeometry, 1, courseOneFormGrad, edgeMappingsPairs);

    // Fix the singularity indices
    for (Face f : globalGeometry.mesh.faces()) {
        if (stripeSingularities[f] != 0) {
            double orient = 0; // sigma . omega
            for (Halfedge he : f.adjacentHalfedges()) {
                
                double sigma = stripeValues[he.next().corner()] - stripeValues[he.corner()]; // stripe 1-form
                if (he.next() == f.halfedge())
                    sigma += 2 * stripeSingularities[f] * PI;

                double omega = omegaWaleGlobal[he.edge()]; // input vector field 1-form
                if (!he.orientation()) // if half-edge does not share orientation of its edge
                    omega *= -1;
                
                orient += omega * sigma;
            }
            if (orient > 0)
                stripeSingularities[f] *= -1;
        }
    }    

    //EdgeData<double> omegaWaleGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, omegaWaleGlobal, edgeMap);
    //Eigen::Map<Eigen::VectorXd> omegaWaleGluedEig(omegaWaleGlued.raw().data(), (gluedGeometry.mesh).nEdges());
    //std::vector<double> modelMatchingTermsWale(omegaWaleGluedEig.data(), omegaWaleGluedEig.data() + omegaWaleGluedEig.rows());
    Eigen::Map<Eigen::VectorXi> faceIndicesWaleEig(stripeSingularities.raw().data(), (gluedGeometry.mesh).nFaces());
    //place wale singularities at edges 
    std::vector<int> faceIndicesWaleModel(faceIndicesWaleEig.data(), faceIndicesWaleEig.data() + faceIndicesWaleEig.rows());
    std::vector<std::array<double, 3>> modelFaceGradients;
    std::vector<std::pair<int, int>> singularEdges;
    globalGeometry.requireFaceNormals();
    for (Face f : globalGeometry.mesh.faces()){ 
        courseOneFormGrad[f] = courseOneFormGrad[f].normalize();
        //rotate the final course gradients 
        courseOneFormGrad[f] = courseOneFormGrad[f].rotateAround(globalGeometry.faceNormals[f], PI/2);
        modelFaceGradients.push_back(std::array{courseOneFormGrad[f][0], courseOneFormGrad[f][1], 
                                                    courseOneFormGrad[f][2]});
        //place wale singularities at edges
        if (stripeSingularities[f] != 0){
            bool toSkip = false;
            for (Edge e : f.adjacentEdges()){
                //skip boundary edges and edges on faces that already have a singularity in the course or wale direction
                if (e.isBoundary() || waleSingularEdgesGlobal[e] != 0 || courseSingularEdgesGlobal[e] != 0) toSkip = true;
            }
            if (toSkip) continue;
            int edge = findSingularEdgeFromSingularFace(globalGeometry, f.getIndex(), courseOneFormGrad[f], 0., timeFunctionGlobal, 0.);
            if (edge == -1) continue;//couldn't find a very well-aligned edge
            int index = stripeSingularities[f] > 0 ? 1 : -1;
            //reduce the singular index to +/- 1 (cause of Autoknit constraints)
            singularEdges.push_back(std::make_pair(edgeMap[edge], index));
            waleSingularEdgesGlobal[edge] = index;
            
        }
    } 
    Model modelWale; 
    modelWale.setPeriod(period);
    //modelWale.setMatchingTerms(modelMatchingTermsWale);
    modelWale.setWaleBdyPathConstraints(globalBdyConditions.waleBdyPathConstraints);
    modelWale.setFaceIndices(faceIndicesWaleModel);
    modelWale.setFaceGradients(modelFaceGradients);
    modelWale.setSingularEdges(singularEdges);

    HalfedgeData<double> sigmaWaleGlued = computeWaleOneForm(globalGeometry, gluedGeometry, modelWale, G, vertexMap);
    CornerData<double> stripeValuesOneFormGlued;
    FaceData<int> stripeIndicesOneFormGlued;
    std::tie(stripeValuesOneFormGlued, stripeIndicesOneFormGlued) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, sigmaWaleGlued, period);

    return std::tie(stripeValuesOneFormGlued, waleSingularEdgesGlobal);
}

/**
//implement the course harmonic 1-form optimization 
//the below two function use the energy min cot_e ||\sigma_e||^2
std::tuple<CornerData<double>, EdgeData<double>> implCourseHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    VertexData<double>& globalTimeFunction, FaceData<Vector3>& globalTimeFunctionGradientsNormalized,
                                                                    std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, polyscope::SurfaceMesh& psMesh,
                                                                    globalBoundaryConditions& boundaryConditions, double period,
                                                                    Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double, Eigen::RowMajor>& G){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    //curl per edge
    EdgeData<double> edgeCurl(globalMesh);
    //curl per face - average the edge curl onto faces
    FaceData<double> faceCurl(globalMesh);
    //striping 1-form quantities
    HalfedgeData<double> gluedSigmaTilde(gluedMesh);
    //edge singularities 
    EdgeData<double> edgeSingularities(globalMesh, 0);
    FaceData<double> faceSingularities(globalMesh, 0);

    //face singularity pairs 
    std::pair<int, int> singFacePair;

    //hashed iso values we've already used
    //map: hashedIsoVal -> count
    std::map<int, int> hashedUsedIsoVals;
    //singular edges in the gurobi optimization
    std::vector<std::pair<int, int>> singularEdges;
    //make a map of singular edges we've seen so far 
    std::map<int, int> seenEdges;
    //number of singularity pairs  
    int numPairs = 0;
    //max number of singularity pairs to insert 
    int maxPairs = 5;
    //gurobi model we will be solving 
    Model model;
    //objective values
    double oldObj, currObj;
    //striping information 
    //global data
    CornerData<double> oldStripeValuesSigmaCourse(gluedMesh);
    CornerData<double> newStripeValuesSigmaCourse(gluedMesh);
    FaceData<int> stripeIndicesSigmaCourse(globalMesh);
    //stripe curve network information
    std::vector<Vector3> positionsCourse;
    std::vector<std::array<int, 2>> edgesCourse;
    //isoval we'll be tracing
    double isoVal;

    //path constraints in the optimization
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints; 

    //set up the boundary->boundary path information 
    //require edge lengths 
    globalGeometry.requireEdgeLengths();
    for (int i = 0; i < boundaryConditions.bdyBdyPathConstraints.size(); i++){
        //visualizing bdy-bdy edge constraints
        EdgeData<double> bdyBdyPath(globalMesh);
        double pathLength = 0.;
        for (Edge e : globalMesh.edges()){
            bdyBdyPath[e] = boundaryConditions.bdyBdyPathConstraints[i][edgeMap[e.getIndex()]];
            if (std::fabs(bdyBdyPath[e]) > 0 && !e.isBoundary()){
                pathLength += globalGeometry.edgeLengths[e];
            } 
        }
        psMesh.addEdgeScalarQuantity("bdy bdy path " + std::to_string(i), bdyBdyPath);
        edgePathConstraints.push_back(std::make_pair(boundaryConditions.bdyBdyPathConstraints[i], ((pathLength)/period)));
    }

    model.setPeriod(period);
    model.setBdyEdges(boundaryConditions.courseBdyEdges);
    //model.setEdgePathConstraints(edgePathConstraints);

    //gradients to use for curl computation 
    FaceData<Vector3> gradients = globalTimeFunctionGradientsNormalized;
    //this is very expensive, ugh
    std::vector<std::array<double, 3>> grads;
    for (Face f : gluedMesh.faces()){
        grads.push_back(std::array<double, 3>{gradients[f][0], gradients[f][1], gradients[f][2]});
    }
    model.setFaceGradients(grads);

    //solve the model without any singularities 
    std::tie(gluedSigmaTilde, currObj) = computeHarmonicCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    std::tie(newStripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
    std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, newStripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
    auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
    oldObj = currObj;
    oldStripeValuesSigmaCourse = newStripeValuesSigmaCourse;
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);

    

    while (numPairs < maxPairs){
        //this is very expensive, ugh
        std::vector<std::array<double, 3>> grads;
        for (Face f : gluedMesh.faces()){
            grads.push_back(std::array<double, 3>{gradients[f][0], gradients[f][1], gradients[f][2]});
        }
        model.setFaceGradients(grads);
        //compute curl per edge of the gradient field (in the global setting)
        edgeCurl = computeEdgeCurl(globalGeometry, gradients);
        //compute curl per face (simply average edge curl onto faces)
        faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
        if (numPairs == 0){//first pair of singular faces
            //find pair max/min curl faces over the entire mesh 
            singFacePair = findFaceSingularityPair(globalGeometry, gluedGeometry, V, F, faceCurl,
                                            globalTimeFunction, psMesh,
                                            hashedUsedIsoVals, 0.0, numPairs, true);
        }
        else{
            isoVal = findIsoValWithMaxFaceCurl(V, F, globalTimeFunction, faceCurl, hashedUsedIsoVals);
            if (std::fabs(isoVal - (-1.0) < 1e-15)){
                std::cout << "Breaking cause we can't find any more sensible isovalues " << std::endl;
                return std::tie(oldStripeValuesSigmaCourse, edgeSingularities);
            } 
            //find pair max/min curl faces over the entire mesh 
            singFacePair = findFaceSingularityPair(globalGeometry, gluedGeometry, V, F, faceCurl,
                                            globalTimeFunction, psMesh,
                                            hashedUsedIsoVals, isoVal, numPairs, false);
        }
        int posEdge = findSingularEdgeFromSingularFace(globalGeometry, singFacePair.first, globalTimeFunctionGradientsNormalized[singFacePair.first]);
        int negEdge = findSingularEdgeFromSingularFace(globalGeometry, singFacePair.second, globalTimeFunctionGradientsNormalized[singFacePair.second]);
        //don't select edges we've seen before
        if (seenEdges.count(edgeMap[posEdge]) > 0 || seenEdges.count(edgeMap[negEdge]) > 0){
            hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
            continue;
        }
        faceSingularities[singFacePair.first] = 1.0;
        faceSingularities[singFacePair.second] = -1.0;
        edgeSingularities[globalMesh.edge(posEdge)] = 1.0;
        edgeSingularities[globalMesh.edge(negEdge)] = -1.0;
        std::cout << "glued pos edge = " << edgeMap[posEdge] << std::endl;
        std::cout << "glued neg edge = " << edgeMap[negEdge] << std::endl;
        //pair of singular edges
        singularEdges.push_back(std::make_pair(edgeMap[negEdge], -1));
        singularEdges.push_back(std::make_pair(edgeMap[posEdge], 1));
        //add it to the map 
        seenEdges[edgeMap[negEdge]] = 1;
        seenEdges[edgeMap[posEdge]] = 1;
        model.setSingularEdges(singularEdges);
        numPairs++;
        //solve the optimization problem 
        //sigmaTilde is in the glued mesh setting 
        std::tie(gluedSigmaTilde, currObj) = computeHarmonicCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
        std::tie(newStripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, newStripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numPairs) + 
                                                    " singularity pairs", positionsCourse, edgesCourse);
        //normalize and update the gradients
        gradients = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
        for (Face f : globalMesh.faces()){
            gradients[f] = gradients[f].normalize();
        }
        oldObj = currObj;
        oldStripeValuesSigmaCourse = newStripeValuesSigmaCourse;
        courseStripes -> setRadius(0.001);
        courseStripes -> setEnabled(false);
    }

    return std::tie(oldStripeValuesSigmaCourse, edgeSingularities);
    
}

//solve the optimization problem for the harmonic 1-form
std::tuple<HalfedgeData<double>, double> computeHarmonicCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    HalfedgeData<double> gluedOneForm(gluedMesh);
    double objectiveVal;
    int numSingularFaces = 0;

    //generate the DEC operators over the glued mesh
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<std::array<double, 3>> comparisonGrad = gbModel.getFaceGradients();

    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();
    std::cout << "Solving model with " << singularEdges.size() / 2.0 << " singularity pairs..." << std::endl;
    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        //don't log to console 
        env.set(GRB_IntParam_OutputFlag, 0);
        //env.set("LogFile", "1-form computation.log");
        env.start();
        // Create an empty model
        GRBModel model = GRBModel(env);
        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 0.5);

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

        //constraint: sigma at boundary edges should  be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        //constraint: nP over every face = 0
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0);
        }

        //constraint: specify singular halfedges and also specify that form values 
        //are equal across non-singular halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            //invert the sign of the constraint
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] ==  p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //constraint: specify path constraints in the optimization 
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
            u[f.halfedge().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()]);
            u[f.halfedge().next().next().vertex().getIndex()] = (sigma[f.halfedge().getIndex()] + sigma[f.halfedge().next().getIndex()]);

            GRBLinExpr currGradU = 0.0;
            //GRBVar currGradU = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][2] = currGradU;
        }
        //set up the objective term
        GRBQuadExpr obj = 0;        
        //setting the objective to be min cot_e||\sigma||^2
        // for (Halfedge he : gluedMesh.halfedges()){
        //     obj +=  gluedGeometry.edgeCotanWeights[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
        // }

        //setting the objective to be min ||\delta sigma - \nabla h / ||\nabla h|| ||^2
        for (Face f : gluedMesh.faces()){
            obj +=  ((gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0]) * (gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0])) 
                    + ((gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1]) * (gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1])) 
                    + ((gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2]) * (gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2])); 
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }
    
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return std::tie(gluedOneForm, objectiveVal);
}
*/

//compute edge curl in the global setting 
EdgeData<double> computeEdgeCurl(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    EdgeData<double> globalCurl(globalMesh);
    globalGeometry.requireEdgeLengths();

    for (Edge e : globalMesh.edges()){
        if (e.isBoundary()) continue;
        Vector3 heVec = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
        Vector3 heTwinVec = globalGeometry.vertexPositions[e.halfedge().tailVertex()] - globalGeometry.vertexPositions[e.halfedge().tipVertex()];
        globalCurl[e] = (dot(globalFaceGradients[e.halfedge().face()].normalize(), heVec) + dot(globalFaceGradients[e.halfedge().twin().face()].normalize(), heTwinVec)) / 
                                globalGeometry.edgeLengths[e];
    }

    return globalCurl;
}


//compute face curl by averaging edge curl over the edges in a face 
FaceData<double> computeAverageEdgeCurlonFaces(VertexPositionGeometry& globalGeometry, EdgeData<double>& edgeCurl){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    FaceData<double> faceCurl(globalMesh);
    globalGeometry.requireFaceAreas();

    for (Face f : globalMesh.faces()){
        double sum = 0.0;
        for (Edge e : f.adjacentEdges()){
            sum += edgeCurl[e];
        }
        faceCurl[f] = sum/(3. * globalGeometry.faceAreas[f]);
    }
    return faceCurl;   
}


//find max/min curl face for a given isoline of the TIME FUNCTION 
//would probably want to change tracing the level sets of the time function with level sets of the harmonic 1-form
std::vector<std::pair<int, int>> findFaceSingularityPairs(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, FaceData<double>& curl,
                                            VertexData<double>& globalTimeFunction, polyscope::SurfaceMesh& psMesh,
                                            std::map<int, int>& hashedUsedIsoVals, FaceData<double>& faceSingularities, FaceData<int>& forbiddenFaces,
                                            double isoVal, int numPairs, bool useAllFaces){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    std::vector<std::pair<int, int>> singFacePairs;
    double eps = 1e-8;
    
    if (useAllFaces){//find max/min curl over the pairs of all isolines over the mesh
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        int maxFace, minFace = -1;
        //find max curl over entire mesh 
        //find face with max  curl
        for (Face f : globalMesh.faces()){
            if (curl[f] > maxCurl){
                maxCurl = curl[f];
                maxFace = f.getIndex();
            }
        }
        //average isovalue on that face
        double isoVal = (globalTimeFunction[globalMesh.face(maxFace).halfedge().tailVertex()] + globalTimeFunction[globalMesh.face(maxFace).halfedge().next().tailVertex()]
                            + globalTimeFunction[globalMesh.face(maxFace).halfedge().next().next().tailVertex()])/3.0;
        hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        //auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        //isoline -> setRadius(0.001);
        //isoline -> setEnabled(false);

        for (int i = 0; i < facesPerComponent.size(); i++){
            std::vector<int> facesInComponent = facesPerComponent[i];
            maxCurl = -DBL_MAX;
            minCurl = DBL_MAX;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                //ignore faces where the current face or an adjacent face is singular
                if (std::abs(faceSingularities[currFace]) > eps
                    || (std::abs(faceSingularities[currFace.halfedge().twin().face()]) > eps
                    || std::abs(faceSingularities[currFace.halfedge().next().twin().face()]) > eps
                    || std::abs(faceSingularities[currFace.halfedge().next().next().twin().face()]) > eps)) continue;
                //igore faces that are forbidden 
                if (forbiddenFaces[currFace] > 0){
                    continue;   
                }
                if (curl[currFace] < minCurl){
                    minCurl = curl[currFace];
                    minFace = currFace.getIndex();
                }
            }
            singFacePairs.push_back(std::make_pair(maxFace, minFace));
        }
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        //auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        //isoline -> setEnabled(false);
        //isoline -> setRadius(0.001);
        for (int i = 0; i < facesPerComponent.size(); i++){
            std::vector<int> facesInComponent = facesPerComponent[i];
            double maxCurl = -DBL_MAX;
            double minCurl = DBL_MAX;
            int maxFace, minFace = -1;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                //ignore faces where the current face or an adjacent face is singular
                if (std::abs(faceSingularities[currFace]) > eps
                    || (std::abs(faceSingularities[currFace.halfedge().twin().face()]) > eps
                    || std::abs(faceSingularities[currFace.halfedge().next().twin().face()]) > eps
                    || std::abs(faceSingularities[currFace.halfedge().next().next().twin().face()]) > eps)) continue;
                //igore faces that are forbidden 
                if (forbiddenFaces[currFace] > 0){
                    continue;   
                }
                if (curl[currFace] > maxCurl){
                    //only select faces where the current face and adjacent face is not singular
                    maxCurl = curl[currFace];
                    maxFace = currFace.getIndex();
                    
                }
                if (curl[currFace] < minCurl){
                    //only select faces where the current face and adjacent face is not singular
                    minCurl = curl[currFace];
                    minFace = currFace.getIndex();
                }       
            }
            
            singFacePairs.push_back(std::make_pair(maxFace, minFace));
        }
    }
    return singFacePairs;
}

//find edge singularity pair 
std::vector<std::pair<int, int>> findEdgeSingularityPairs(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, EdgeData<double>& curl,
                                            VertexData<double>& globalTimeFunction, polyscope::SurfaceMesh& psMesh,
                                            std::map<int, int>& hashedUsedIsoVals, std::map<int, int>& usedEdges, FaceData<double>& faceSingularities, 
                                            double isoVal, int numPairs, bool useAllEdges){


    SurfaceMesh& globalMesh = globalGeometry.mesh;
    std::vector<std::pair<int, int>> singEdgePairs;
    double eps = 1e-8;
    
    if (useAllEdges){//find max/min curl over the pairs of all isolines over the mesh
        double maxCurl = -DBL_MAX;
        double minCurl = DBL_MAX;
        int maxEdge, minEdge = -1;
        //find max curl over entire mesh 
        //find edge with max  curl
        for (Edge e : globalMesh.edges()){
            if (curl[e] > maxCurl && !usedEdges.count(e.getIndex()) != 1){
                maxCurl = curl[e];
                maxEdge = e.getIndex();
            }
        }
        //average isovalue on the max edge
        double isoVal = 0.5 * (globalTimeFunction[globalMesh.edge(maxEdge).halfedge().tailVertex()] +
                                     globalTimeFunction[globalMesh.edge(maxEdge).halfedge().tipVertex()]);
        hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        //auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        //isoline -> setRadius(0.001);
        //isoline -> setEnabled(false);

        for (int i = 0; i < facesPerComponent.size(); i++){
            std::vector<int> facesInComponent = facesPerComponent[i];
            maxCurl = -DBL_MAX;
            minCurl = DBL_MAX;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                for (Edge e : currFace.adjacentEdges()){
                    //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    || (isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] < minCurl && !usedEdges.count(e.getIndex()) != 1){
                            minCurl = curl[e];
                            minEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_pair(maxEdge, minEdge));
        }
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        //auto isoline = polyscope::registerCurveNetwork("Isoline for pair " + std::to_string(numPairs + 1), iV, iE);
        //isoline -> setEnabled(false);
        //isoline -> setRadius(0.001);
        for (int i = 0; i < facesPerComponent.size(); i++){
            std::vector<int> facesInComponent = facesPerComponent[i];
            double maxCurl = -DBL_MAX;
            double minCurl = DBL_MAX;
            int maxEdge, minEdge = -1;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                for (Edge e : currFace.adjacentEdges()){
                    //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                        ||(isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] > maxCurl && !usedEdges.count(e.getIndex()) > 0){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                    }
                    //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                        ||(isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] < minCurl && !usedEdges.count(e.getIndex()) > 0){
                            minCurl = curl[e];
                            minEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_pair(maxEdge, minEdge));
        }
    }
    return singEdgePairs;
}

//get the vertices (of a curve network), edges (of a curve network) and faces that a particular isovalue of the time function passes through 
//generate isolines for the time function given a specific isoVal
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, std::vector<int>> getIsoLine(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& timeFunction, double isoVal){

    Eigen::MatrixXd S(V.size(), 1);
    for (int i = 0; i < V.size(); i++){
        S(i, 0) = timeFunction[i];
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

//find the isoval with max average curl in the face setting
double findIsoValWithMaxFaceCurl(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& globalTimeFunction, 
                                FaceData<double>& curl, std::map<int, int>& hashedUsedIsoVals, FaceData<int>& forbiddenFaces, double stepSize){

    double end = 1.0;
    double curr = stepSize;
    double maxDeviation = -DBL_MAX;
    double maxDeviationIsoVal = -1.0;
    double eps = 1e-2;
    bool skipFlag = false;
    double currAvgDeviation = 0.0;
    double currDeviationSum = 0.0;

    //std::cout << "step size = " << stepSize << std::endl;

    while (curr < (end - eps)){

        if (hashedUsedIsoVals.count(hashFloatQuantized(curr)) > 0){
            curr += stepSize;
            continue;
        }
    
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, curr);

        // auto isoline = polyscope::registerCurveNetwork("sampled isoline " + std::to_string(curr) , iV, iE);
        // isoline->setRadius(0.001);
        // isoline->setEnabled(false);

        //if the isoval contains a single forbidden face
        //skip it and continue 
        for (int i = 0; i < f.size(); i++){
            if (forbiddenFaces[f[i]] > 0){
                curr += stepSize;
                continue;
            }
        }

        //reset values
        currDeviationSum = 0.0;
        currAvgDeviation = 0.0;
        int numEdges = 0;
        for (int i = 0; i < f.size(); i++){
            currDeviationSum += std::fabs(curl[f[i]]);
        }
        currAvgDeviation = currDeviationSum / f.size();
        if (currAvgDeviation > maxDeviation){
            maxDeviation = currAvgDeviation;
            maxDeviationIsoVal = curr;
        }
        curr += stepSize;
    }

    hashedUsedIsoVals[hashFloatQuantized(maxDeviationIsoVal)] = 1;
    return maxDeviationIsoVal;
}

//@debugging 
//find isoval with maximum average edge curl 
double findIsoValWithMaxAvgEdgeCurl(VertexPositionGeometry& globalGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& globalTimeFunction, 
                                    EdgeData<double>& curl, std::map<int, int>& hashedUsedIsoVals, double stepSize){

    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    double end = 1.0;
    double curr = 0.01;
    double maxDeviation = -DBL_MAX;
    double maxDeviationIsoVal = -1.0;
    double eps = 1e-5;
    bool skipFlag = false;
    double currAvgDeviation = 0.0;
    double currDeviationSum = 0.0;
    //std::cout << "step size = " << stepSize << std::endl;
    //stepSize = 0.01;
    
    while (curr < (end - eps)){

        if (hashedUsedIsoVals.count(hashFloatQuantized(curr)) > 0){
            curr += stepSize;
            continue;
        }
    
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, curr);

        // auto isoline = polyscope::registerCurveNetwork("sampled isoline " + std::to_string(curr) , iV, iE);
        // isoline->setRadius(0.001);
        // isoline->setEnabled(true);

        //reset values
        currDeviationSum = 0.0;
        currAvgDeviation = 0.0;
        int numEdges = 0;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            for (Edge e : currFace.adjacentEdges()){
                //check if current isoline crosses this edge
                if ((curr > globalTimeFunction[e.halfedge().tailVertex()] && curr < globalTimeFunction[e.halfedge().tipVertex()])
                || (curr > globalTimeFunction[e.halfedge().tipVertex()] && curr < globalTimeFunction[e.halfedge().tailVertex()])){
                    currDeviationSum += std::fabs(curl[e]);
                    numEdges++;
                }
            }
        }
        //take the average high curl per edge
        //currAvgDeviation = currDeviationSum / numEdges;
        //take the absolute value of the sum of the curl
        currAvgDeviation = currDeviationSum;
        if (currAvgDeviation > maxDeviation){
            maxDeviation = currAvgDeviation;
            maxDeviationIsoVal = curr;
        }
        curr += stepSize;
    }
    return maxDeviationIsoVal;
}




//@clean 
//the below two functions use the energy min ||\del sigma_i - \nabla \sigma_{i - 1} / ||\nabla \sigma_{i - 1}|| ||^2
//if i = 1, use \nabla h / ||\nabla h||
std::tuple<CornerData<double>, EdgeData<double>> implCourseHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    VertexData<double>& globalTimeFunction, FaceData<Vector3>& globalTimeFunctionGradientsNormalized,
                                                                    std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, polyscope::SurfaceMesh& psMesh,
                                                                    globalBoundaryConditions& boundaryConditions, double period,
                                                                    Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double, Eigen::RowMajor>& G,
                                                                    FaceData<Vector3>& courseOneFormGrad, std::map<int, std::vector<Halfedge>>& gluedOneRingMap){
    
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    
    //curl per vertex 
    VertexData<double> vertexCurl(globalMesh, 0.0);
    //curl per edge
    EdgeData<double> edgeCurl(globalMesh, 0.0);
    //curl per face - average the edge curl onto faces
    FaceData<double> faceCurl(globalMesh, 0.0);
    //striping 1-form quantities
    HalfedgeData<double> gluedSigmaTilde(gluedMesh, 0.0);
    //harmonic sigma info 
    HalfedgeData<double> harmonicSigmaTilde(gluedMesh, 0.0);
    double harmonicSigmaObj;
    //virtual sigma info
    HalfedgeData<double> virtualSigmaTilde(gluedMesh, 0.0);
    double virtualSigmaObj;
    //edge singularities 
    EdgeData<double> edgeSingularities(globalMesh, 0.0);
    FaceData<double> faceSingularities(globalMesh, 0.0);
    //forbidden faces - faces touched by some isoline can't be used again
    //this is too restrictive
    FaceData<int> forbiddenFaces(globalMesh, 0.0);
    //set the boundary loops to 0 as well
    for (BoundaryLoop b : globalMesh.boundaryLoops()){
        Halfedge he = b.halfedge();
        faceSingularities[he.face()] = 0.0;
    }
    //gradient of the one form used in each iterative optimization
    FaceData<Vector3> gradSigmaTilde(globalMesh, Vector3{0.0, 0.0, 0.0});
    //gradient of the one form after subtracting off the impluse function 
    FaceData<Vector3> adjustedGradSigmaTilde(globalMesh, Vector3{0.0, 0.0, 0.0});
    //face singularity pairs 
    std::vector<std::pair<int, int>> singFacePairs;
    //hashed iso values we've already used
    //map: hashedIsoVal -> count
    std::map<int, int> hashedUsedIsoVals;
    //singular edges in the gurobi optimization
    std::vector<std::pair<int, int>> singularEdges;
    //sing edge pairs (max/min curl edges)
    std::vector<std::pair<int, int>> singEdgePairs;
    //make a map of singular edges we've used so far 
    //so we don't re-use edges
    std::map<int, int> usedEdges;
    //path constraints
    std::vector<double> globalPath;
    std::vector<double> gluedPath;
    //number of runs of the optimization
    int numRuns = 0;
    //max number of singularity pairs to insert 
    int maxRuns = 2;
    //threshold for alignment with the gradient 
    double threshold = 0.85;
    //gurobi model we will be solving 
    Model model;
    //harmonic 1-form model that doesn't collapse 
    Model harmonicModel;
    //edge path constraints in the optimization 
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints;
    //objective values
    double oldObj, currObj;
    //distance from unit norm 
    double oldDistance = DBL_MAX;
    double newDistance = 0.;
    //striping information 
    //global data
    CornerData<double> stripeValuesSigmaCourse(gluedMesh);
    FaceData<int> stripeIndicesSigmaCourse(globalMesh);
    //stripe curve network information
    std::vector<Vector3> positionsCourse;
    std::vector<std::array<int, 2>> edgesCourse;
    //unique vertices and unique edges
    std::vector<Vector3> uniquePos;
    std::vector<std::array<int, 2>> uniqueEdges;

    //isoval we'll be tracing
    double isoVal;
    //skip flag - to skip the current iteration
    bool skipFlag = false;

    //set the period for the model 
    model.setPeriod(period);
    model.setBdyEdges(boundaryConditions.courseBdyEdges);
    //gradients to use for curl computation 
    gradSigmaTilde = globalTimeFunctionGradientsNormalized;
    //this is very expensive, ugh
    std::vector<std::array<double, 3>> grads;
    for (Face f : gluedMesh.faces()){
        grads.push_back(std::array<double, 3>{gradSigmaTilde[f][0], gradSigmaTilde[f][1], gradSigmaTilde[f][2]});
    }
    model.setFaceGradients(grads);

    //find the step size to sample level sets at
    double avgSum = 0.;
    double stepSize = 0.;
    for (int i = 0; i < boundaryConditions.bdyBdyPathConstraints.size(); i++){
        //visualizing bdy-bdy edge constraints
        EdgeData<double> bdyBdyPath(globalMesh);
        double sum = 0.;
        int ctr = 0;
        for (Edge e : globalMesh.edges()){
            bdyBdyPath[e] = boundaryConditions.bdyBdyPathConstraints[i][edgeMap[e.getIndex()]];
            if (std::fabs(bdyBdyPath[e]) > 1e-10){
                Vertex startVertex, endVertex;
                startVertex = globalTimeFunction[e.halfedge().tipVertex()] > globalTimeFunction[e.halfedge().tailVertex()] ? 
                            e.halfedge().tailVertex() : e.halfedge().tipVertex();
                endVertex = globalTimeFunction[e.halfedge().tipVertex()] > globalTimeFunction[e.halfedge().tailVertex()] ?
                            e.halfedge().tipVertex() :  e.halfedge().tailVertex();
                sum += globalTimeFunction[endVertex] - globalTimeFunction[startVertex];
                ctr++;
            }
        }
        avgSum += sum/ctr;
        //psMesh.addEdgeScalarQuantity("bdy bdy path " + std::to_string(i), bdyBdyPath);
    }
    stepSize = avgSum / boundaryConditions.bdyBdyPathConstraints.size();


    //set up the boundary->boundary path information 
    //require edge lengths 
    // globalGeometry.requireEdgeLengths();
    // std::vector<std::pair<std::vector<double>, double>> harmonicPathConstraints;
    // for (int i = 0; i < boundaryConditions.bdyBdyPathConstraints.size(); i++){
    //     //visualizing bdy-bdy edge constraints
    //     EdgeData<double> bdyBdyPath(globalMesh);
    //     double pathLength = 0.;
    //     for (Edge e : globalMesh.edges()){
    //         bdyBdyPath[e] = boundaryConditions.bdyBdyPathConstraints[i][edgeMap[e.getIndex()]];
    //         if (std::fabs(bdyBdyPath[e]) > 0 && !e.isBoundary()){
    //             pathLength += globalGeometry.edgeLengths[e];
    //         } 
    //     }
    //     psMesh.addEdgeScalarQuantity("bdy bdy path " + std::to_string(i), bdyBdyPath);
    //     harmonicPathConstraints.push_back(std::make_pair(boundaryConditions.bdyBdyPathConstraints[i], ((pathLength)/period)));
    // }
    // //set non-collapse constraint
    // harmonicModel.setEdgePathConstraints(harmonicPathConstraints);


    //construct edge weights for halfedge path
    FaceData<Vector3> rotatedFaceGradients = clockWiseRotatedGradients(globalGeometry, globalTimeFunctionGradientsNormalized);
    double maxDotProd = maximumDotProduct(globalGeometry, rotatedFaceGradients);
    HalfedgeData<double> gluedHeWeights = constructGluedHalfedgeWeights(globalGeometry, gluedGeometry, rotatedFaceGradients, maxDotProd);


    //solve the model without any singularities 
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    numRuns++;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, 1.0 * period);
    std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, 1.0 * period);
    //remove the duplicate vertices in the curve network
    std::tie(uniquePos, uniqueEdges) = removeCurveNetworkDuplicatedVertices(globalGeometry, positionsCourse, edgesCourse);
    auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after duplicate removal " + std::to_string(numRuns - 1) + 
                                                    " singularity insertions", uniquePos, uniqueEdges);
    //find the connected components in the unique graph
    findCurveNetworkConnectedComponents(globalGeometry, uniquePos, uniqueEdges);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);
    gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    //update the gradients for the next iteration of the model 
    model.setFaceGradients(gradSigmaTilde);
    std::cout << "distance from unit norm after " << std::to_string(numRuns - 1) << " singularity insertions " << computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde) << std::endl;
    
    //compute curl quantities without impulse function
    vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                    gradSigmaTilde, gluedOneRingMap);
    edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (before subtracting)", edgeCurl);
    //compute virtual sigma
    std::tie(virtualSigmaTilde, virtualSigmaObj) = computeVirtualSigma(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    //subtract off the impulse function
    gluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
    adjustedGradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    //compute curl quantities after impulse function
    vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                    adjustedGradSigmaTilde, gluedOneRingMap);
    edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (after subtracting)", edgeCurl);
    
    while(numRuns < 1){
        //new singularity pair iteration
        isoVal = findIsoValWithMaxAvgEdgeCurl(globalGeometry, V, F, globalTimeFunction, edgeCurl, hashedUsedIsoVals, stepSize);
        std::cout << "next isoVal = " << isoVal << std::endl;
        //find edge singularity pair 
        singEdgePairs = findEdgeSingularityPairs(globalGeometry, gluedGeometry,  V,  F,  edgeCurl,
                                            globalTimeFunction, psMesh,
                                            hashedUsedIsoVals, usedEdges, faceSingularities, 
                                            isoVal, numRuns, false);
        for (std::pair<int, int> singEdgePair : singEdgePairs){
            Eigen::MatrixXd iV;
            Eigen::MatrixXd iE;
            std::vector<int> f;
            std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
            auto isoline = polyscope::registerCurveNetwork("Isoline for " + std::to_string(numRuns) + " singularity", iV, iE);
            isoline -> setRadius(0.001);
            isoline -> setEnabled(false);
            edgeSingularities[globalMesh.edge(singEdgePair.first)] = 1.0;
            edgeSingularities[globalMesh.edge(singEdgePair.second)] = -1.0;
            //pair of singular edges
            singularEdges.push_back(std::make_pair(edgeMap[singEdgePair.first], -1));
            singularEdges.push_back(std::make_pair(edgeMap[singEdgePair.second], 1));
            std::tie(globalPath, gluedPath) = constructEdgePath(globalGeometry, gluedGeometry, globalMesh.edge(singEdgePair.first), globalMesh.edge(singEdgePair.second),
                                        vertexMap, edgeMap, globalTimeFunctionGradientsNormalized, gluedHeWeights);
            //don't retake edges from another path
            updateGluedHalfedgeWeights(globalGeometry, gluedGeometry, gluedPath, gluedHeWeights);
            psMesh.addEdgeScalarQuantity("path for " + std::to_string(numRuns) + " singularity", globalPath);
            edgePathConstraints.push_back(std::make_pair(gluedPath, 0.));
            model.setEdgePathConstraints(edgePathConstraints);
            model.setSingularEdges(singularEdges);
            //solve the model with 1 pair of  singularities 
            std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
            numRuns++;
            std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, 1.0 * period);
            std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, 1.0 * period);
            //remove the duplicates from the curve network
            std::tie(uniquePos, uniqueEdges) = removeCurveNetworkDuplicatedVertices(globalGeometry, positionsCourse, edgesCourse);
            courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after duplicate removal " + std::to_string(numRuns - 1) + 
                                                    " singularity insertions", uniquePos, uniqueEdges);
            courseStripes -> setRadius(0.001);
            courseStripes -> setEnabled(false);
            gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
            //update the gradients for the next round of the model 
            model.setFaceGradients(gradSigmaTilde);
            std::cout << "distance from unit norm after " << std::to_string(numRuns - 1) << " singularity insertions " << computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde) << std::endl;
            //find the connected components in the unique graph
            findCurveNetworkConnectedComponents(globalGeometry, uniquePos, uniqueEdges);

            //compute curl quantities without accounting for impulse function
            vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                    gradSigmaTilde, gluedOneRingMap);
            edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
            psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (before subtracting)", edgeCurl);
            //compute virtual sigma
            std::tie(virtualSigmaTilde, virtualSigmaObj) = computeVirtualSigma(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
            //subtract off the impulse function
            gluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
            adjustedGradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
            //compute curl quantities after accounting for impulse function
            vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                    adjustedGradSigmaTilde, gluedOneRingMap);
            edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
            psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (after subtracting)", edgeCurl);
        }
        //add the isovalue we've used to the map
        hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});

    }

    // solve the harmonic model without any singularities 
    // std::tie(harmonicSigmaTilde, harmonicSigmaObj) = computeHarmonicCourseOneForm(globalGeometry, gluedGeometry, harmonicModel, vertexMap, G, psMesh);
    // std::tie(newStripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, harmonicSigmaTilde, period);
    // std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, newStripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
    // auto harmonicCourseStripes = polyscope::registerCurveNetwork("harmonic sigma tilde stripes after " + std::to_string(numRuns - 1) + 
    //                                                 " singularity insertions", positionsCourse, edgesCourse);
    // harmonicCourseStripes -> setRadius(0.001);
    // harmonicCourseStripes -> setEnabled(false);

    /** 
     * face-averaged edge singularity setting
    while(numRuns < 5){
        skipFlag = false;
        //non-harmonic sigma tilde
        gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
        vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                        gradSigmaTilde, gluedOneRingMap);
        edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
        //edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
        faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
        psMesh.addFaceVectorQuantity("non-harmonic grad(sigma) after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", gradSigmaTilde);
        psMesh.addEdgeScalarQuantity("non-harmonic edge curl after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", edgeCurl);
        //psMesh.addFaceScalarQuantity("non-harmonic face curl after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", faceCurl);
        
        //subtract off the virtual sigma from the non-harmonic sigma
        gluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
        gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
        vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                        gradSigmaTilde, gluedOneRingMap);
        edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
        //edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
        faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
        psMesh.addFaceVectorQuantity("non-harmonic grad(sigma) after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", gradSigmaTilde);
        psMesh.addEdgeScalarQuantity("non-harmonic edge curl after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", edgeCurl);
        //psMesh.addFaceScalarQuantity("non-harmonic face curl after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", faceCurl);

        // harmonic sigma tilde
        // gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, harmonicSigmaTilde);
        // edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
        // faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
        // psMesh.addFaceVectorQuantity("harmonic grad(sigma) after " + std::to_string(numRuns) + " runs before subtracting", gradSigmaTilde);
        // psMesh.addEdgeScalarQuantity("harmonic edge curl after " + std::to_string(numRuns) + " runs before subtracting", edgeCurl);
        // psMesh.addFaceScalarQuantity("harmonic face curl after " + std::to_string(numRuns) + " runs before subtracting", faceCurl);
        // //subtract off the virtual sigma from the harmonic sigma
        // harmonicSigmaTilde = harmonicSigmaTilde - virtualSigmaTilde;
        // gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, harmonicSigmaTilde);
        // edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
        // faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
        // psMesh.addFaceVectorQuantity("harmonic grad(sigma) after " + std::to_string(numRuns) + " runs after subtracting", gradSigmaTilde);
        // psMesh.addEdgeScalarQuantity("harmonic edge curl after " + std::to_string(numRuns) + " runs after subtracting", edgeCurl);
        // psMesh.addFaceScalarQuantity("harmonic face curl after " + std::to_string(numRuns) + " runs after subtracting", faceCurl);

        newDistance = 0.;
        //this is very expensive, ugh
        for (Face f : globalMesh.faces()){
            //compute the squared distance from unit norm
            newDistance += pow(gradSigmaTilde[f].norm() - 1., 2.);
            //normalize before we pass to the model
            gradSigmaTilde[f] = gradSigmaTilde[f].normalize();
            grads.push_back(std::array<double, 3>{gradSigmaTilde[f][0], gradSigmaTilde[f][1], gradSigmaTilde[f][2]});
        }
        // if (newDistance > oldDistance){
        //     std::cout << "number of runs = " << numRuns << std::endl;
        //     std::cout << "(breaking) oldDistance = " << oldDistance << std::endl;
        //     std::cout << "(breaking) newDistance = " << newDistance << std::endl;
        //     return std::tie(oldStripeValuesSigmaCourse, edgeSingularities);
        // }
        std::cout << "number of runs = " << numRuns << std::endl;
        std::cout << "(not breaking) oldDistance = " << oldDistance << std::endl;
        std::cout << "(not breaking) newDistance = " << newDistance << std::endl;
        std::cout << "-----------------------------" << std::endl;
        oldDistance = newDistance;
        model.setFaceGradients(grads);
        isoVal = findIsoValWithMaxFaceCurl(V, F, globalTimeFunction, faceCurl, hashedUsedIsoVals, forbiddenFaces, stepSize);
        if (std::fabs(isoVal - (-1.0) < 1e-15)){//our function couldn't find anymore isovals
            std::cout << "Breaking cause we can't find any more sensible isovalues " << std::endl;
            courseOneFormGrad = gradSigmaTilde;
            return std::tie(oldStripeValuesSigmaCourse, edgeSingularities);
        } 
        singFacePairs = findFaceSingularityPairs(globalGeometry, gluedGeometry, V, F, faceCurl,
                                                globalTimeFunction, psMesh,
                                                hashedUsedIsoVals, faceSingularities, forbiddenFaces,
                                                isoVal, numRuns, false);
        //handle every pair returned 
        for (std::pair<int, int> singFacePair : singFacePairs){
            if (singFacePair.first == -1 || singFacePair.second == -1){
                std::cout << "skipping face " << singFacePair.first << std::endl;
                std::cout << "skipping face " << singFacePair.second << std::endl;
                continue;
            }
            //ensure to not pick edges that pass through another isoline
            int posEdge = findSingularEdgeFromSingularFace(globalGeometry, singFacePair.first, globalTimeFunctionGradientsNormalized[singFacePair.first], threshold, globalTimeFunction, isoVal);
            int negEdge = findSingularEdgeFromSingularFace(globalGeometry, singFacePair.second, globalTimeFunctionGradientsNormalized[singFacePair.second], threshold, globalTimeFunction, isoVal);
            //don't select edges we've seen before
            if (usedEdges.count(edgeMap[posEdge]) > 0 || usedEdges.count(edgeMap[negEdge]) > 0){
                std::cout << "skipping edge " << edgeMap[posEdge] << std::endl;
                std::cout << "skipping edge " << edgeMap[negEdge] << std::endl;
                hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
                skipFlag = true;
                continue;
            }
            if (posEdge == negEdge){//could happen if two singular faces are adjacent to each other
                hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
                skipFlag = true;
                continue;
            }
            if (posEdge == -1 || negEdge == -1){
                hashedUsedIsoVals.insert({hashFloatQuantized(isoVal), 1});
                skipFlag = true;
                continue;//couldn't find a very well aligned edge
            }
            //we're definitely using this isoval so update forbidden faces 
            //updateForbiddenFaces(V, F, globalTimeFunction, isoVal, forbiddenFaces);
            //view the isoline 
            Eigen::MatrixXd iV;
            Eigen::MatrixXd iE;
            std::vector<int> f;
            std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
            auto isoline = polyscope::registerCurveNetwork("Isoline for " + std::to_string(numRuns) + " singularity", iV, iE);
            isoline -> setRadius(0.001);
            isoline -> setEnabled(false);
            faceSingularities[singFacePair.first] = 1.0;
            faceSingularities[singFacePair.second] = -1.0;
            edgeSingularities[globalMesh.edge(posEdge)] = 1.0;
            edgeSingularities[globalMesh.edge(negEdge)] = -1.0;
            //pair of singular edges
            singularEdges.push_back(std::make_pair(edgeMap[negEdge], -1));
            singularEdges.push_back(std::make_pair(edgeMap[posEdge], 1));
            std::tie(globalPath, gluedPath) = constructEdgePath(globalGeometry, gluedGeometry, globalMesh.edge(negEdge), globalMesh.edge(posEdge),
                                        vertexMap, edgeMap, globalTimeFunctionGradientsNormalized, gluedHeWeights);
            //don't retake edges from another path
            updateGluedHalfedgeWeights(globalGeometry, gluedGeometry, gluedPath, gluedHeWeights);
            psMesh.addEdgeScalarQuantity("path for " + std::to_string(numRuns) + " singularity", globalPath);
            edgePathConstraints.push_back(std::make_pair(gluedPath, 0.));
            //add to the harmonic model
            //harmonicPathConstraints.push_back(std::make_pair(gluedPath, 0.));
            //add it to the map 
            usedEdges[edgeMap[negEdge]] = 1;
            usedEdges[edgeMap[posEdge]] = 1;
        }
        if (skipFlag) continue;//if any of the skipping conditions are met, skip the outer iteration as well
        psMesh.addFaceScalarQuantity("face singularities after " + std::to_string(numRuns) + " singularities", faceSingularities);
        psMesh.addEdgeScalarQuantity("edge singularities after " + std::to_string(numRuns) + " singularities", edgeSingularities);
        model.setEdgePathConstraints(edgePathConstraints);
        model.setSingularEdges(singularEdges);
        //harmonicModel.setEdgePathConstraints(harmonicPathConstraints);
        //harmonicModel.setSingularEdges(singularEdges);
        numRuns++;

        //solve the optimization problem 
        //sigmaTilde is in the glued mesh setting 
        std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
        std::tie(newStripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, newStripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numRuns - 1) + 
                                                    " singularity insertions", positionsCourse, edgesCourse);
        oldObj = currObj;
        oldStripeValuesSigmaCourse = newStripeValuesSigmaCourse;
        courseStripes -> setRadius(0.001);
        courseStripes -> setEnabled(false);

        //compute harmonic sigma
        //std::tie(harmonicSigmaTilde, harmonicSigmaObj) = computeHarmonicCourseOneForm(globalGeometry, gluedGeometry, harmonicModel, vertexMap, G, psMesh);
        
        //compute virtual sigma
        std::tie(virtualSigmaTilde, virtualSigmaObj) = computeVirtualSigma(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
        std::tie(newStripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, virtualSigmaTilde, period);
        std::tie(positionsCourse, edgesCourse) = generateIsoLines(globalGeometry, newStripeValuesSigmaCourse, stripeIndicesSigmaCourse, period);
        auto virtualStripes = polyscope::registerCurveNetwork("virtual sigma tilde stripes after " + std::to_string(numRuns - 1) + 
                                                    " runs", positionsCourse, edgesCourse);
        virtualStripes -> setRadius(0.001);
        virtualStripes -> setEnabled(false);
        // gradient of the harmonic virtual sigma
        FaceData<Vector3> gradVirtualSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, virtualSigmaTilde);
        psMesh.addFaceVectorQuantity("grad virtual sigma tilde after " + std::to_string(numRuns) + " runs", gradVirtualSigmaTilde);

    }
    
    //psMesh.addFaceScalarQuantity("forbidden faces", forbiddenFaces);
    //non-harmonic sigma tilde
    gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                        gradSigmaTilde, gluedOneRingMap);
    edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
    //edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
    psMesh.addFaceVectorQuantity("non-harmonic grad(sigma) after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", gradSigmaTilde);
    psMesh.addEdgeScalarQuantity("non-harmonic edge curl after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", edgeCurl);
    //psMesh.addFaceScalarQuantity("non-harmonic face curl after " + std::to_string(numRuns - 1) + " singularities (before subtracting)", faceCurl);

    //subtract off the virtual sigma from the non-harmonic sigma
    gluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
    gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    vertexCurl = computeVertexCurl(globalGeometry, gluedGeometry, 
                                        gradSigmaTilde, gluedOneRingMap);
    edgeCurl = computeEdgeCurl(globalGeometry, gradSigmaTilde);
    //edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    faceCurl = computeAverageEdgeCurlonFaces(globalGeometry, edgeCurl);
    psMesh.addFaceVectorQuantity("non-harmonic grad(sigma) after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", gradSigmaTilde);
    psMesh.addEdgeScalarQuantity("non-harmonic edge curl after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", edgeCurl);
    //psMesh.addFaceScalarQuantity("non-harmonic face curl after " + std::to_string(numRuns - 1) + " singularities (after subtracting)", faceCurl);
    */

    courseOneFormGrad = gradSigmaTilde;
    return std::tie(stripeValuesSigmaCourse, edgeSingularities);
    

}

std::tuple<HalfedgeData<double>, double> computeCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);;
    double objectiveVal;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<std::array<double, 3>> comparisonGrad = gbModel.getFaceGradients();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();

    
    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        //model.getEnv().set(GRB_DoubleParam_TimeLimit, 0.5);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_SolutionLimit, 2);

        //sigma defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //constraint: sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        
        //compute nP over every face 
        //add the second constraint while we're here
        //constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0);
            nP[f.getIndex()] = 0.;
        }

        //constraint: 
        //specify singular halfedges and also specify that form values 
        //are equal across non-singular halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] ==  p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //constraint: add bdy-bdy path constraint 
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
                model.addConstr(pathIntegral == 0.);
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
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][0] = currGradU;
            currGradU = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][1] = currGradU;
            currGradU = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * gluedMesh.nFaces()); it; ++it){
                currGradU += it.value() * u[vertexMap[it.col()]];
            }
            gradU[f.getIndex()][2] = currGradU;
        }


        //set up the objective term
        GRBQuadExpr obj = 0;        
    

        //setting the objective to be min ||\delta sigma - \nabla h / ||\nabla h|| ||^2
        for (Face f : gluedMesh.faces()){
            obj +=  ((gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0]) * (gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0])) 
                    + ((gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1]) * (gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1])) 
                    + ((gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2]) * (gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2]));  
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
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

    return std::tie(gluedOneForm, objectiveVal);
}


//@debugging 
//compute a virtual sigma
std::tuple<HalfedgeData<double>, double> computeVirtualSigma(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh){


    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);;
    double objectiveVal;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();

    
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
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_SolutionLimit, 2);

        //sigma defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //constraint: sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        //add the second constraint while we're here
        //constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0);
        }

        //constraint: 
        //specify singular halfedges and also specify that form values 
        //are equal across non-singular halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] ==  p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

    

        //set up the objective term
        GRBQuadExpr obj = 0;        
    
        //setting the objective to be min ||\sigma||^2
        for (Halfedge he : gluedMesh.halfedges()){
            obj += gluedGeometry.edgeCotanWeights[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
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

    return std::tie(gluedOneForm, objectiveVal);
    
}

//@debugging
//compute a harmonic course 1-form
std::tuple<HalfedgeData<double>, double> computeHarmonicCourseOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);;
    double objectiveVal;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = gbModel.getEdgePathConstraints();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();

    
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
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_SolutionLimit, 2);

        //sigma defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //constraint: sigma at boundary edges should be 0
        for (int bdy_edge_index : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().getIndex()] == 0.0, "Boundary Constraint");
            model.addConstr(sigma[gluedMesh.edge(bdy_edge_index).halfedge().twin().getIndex()] == 0.0, "Boundary Constraint");
        }

        
        //compute nP over every face 
        //add the second constraint while we're here
        //constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        //here it will (d1 * sigma) == 0
        std::vector<GRBLinExpr> nP(gluedMesh.nFaces());
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0);
            nP[f.getIndex()] = 0.;
        }

        //constraint: 
        //specify singular halfedges and also specify that form values 
        //are equal across non-singular halfedges
        EdgeData<int> handled(gluedMesh, 0);
        for (std::pair<int, int> p : singularEdges){
            model.addConstr(sigma[gluedMesh.edge(p.first).halfedge().getIndex()] 
                                    + sigma[gluedMesh.edge(p.first).halfedge().twin().getIndex()] ==  p.second * period);
            handled[gluedMesh.edge(p.first)] = 1;
        }
        for (Halfedge he : gluedMesh.halfedges()){
            if (handled[he.edge()]) continue; //if the halfedge is already handled, continue
            model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
        }

        //constraint: add bdy-bdy path constraint 
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
                model.addConstr(pathIntegral == 0.);
        }

        //set up the objective term
        GRBQuadExpr obj = 0;        
    

        //setting the objective to be min ||\sigma||^2
        for (Halfedge he : gluedMesh.halfedges()){
            obj += gluedGeometry.edgeCotanWeights[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
        }


        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
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

    return std::tie(gluedOneForm, objectiveVal);

}

//update forbidden faces 
//in particular, if some isoline passes through a set of faces ensure we don't select another pair of faces 
//on the same isoline
void updateForbiddenFaces(Eigen::MatrixXd& V, Eigen::MatrixXi& F, VertexData<double>& timeFunction, double isoVal, FaceData<int>& forbiddenFaces){

    Eigen::MatrixXd S(V.size(), 1);
    for (int i = 0; i < V.size(); i++){
        S(i, 0) = timeFunction[i];
    }

    //std::vector<int> passes;
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
            //passes.push_back(f);
            forbiddenFaces[f] = 1;
        }
    }
}

//@debugging
//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes (in global setting)
VertexData<double> computeVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                        FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> curl(globalMesh);
    globalGeometry.requireFaceAreas();
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        for (Halfedge he : gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()];
            //always normalize the field
            sum += dot(hjkVec, field[he.face()].normalize());
        }
        curl[vi] = sum;
    }

    return curl;
}

//@debugging
//average vertex curl onto edges
EdgeData<double> computeVertexAveragedEdgeCurl(VertexPositionGeometry& globalGeometry, VertexData<double>& vertexCurl){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    globalGeometry.requireEdgeLengths();
    EdgeData<double> edgeCurl(globalMesh);
    for (Edge e : globalMesh.edges()){
        double c1 = vertexCurl[e.halfedge().tailVertex()];
        double c2 = vertexCurl[e.halfedge().tipVertex()];
        edgeCurl[e] = (c1 + c2) /(2. * globalGeometry.edgeLengths[e]);
    }
    return edgeCurl;

}

//@debugging
//compute distance from unit norm of a per-face vector field
double computeDistanceFromUnitNorm(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& gradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    double sum = 0.0;

    for (Face f : globalMesh.faces()){
        sum += pow(gradients[f].norm() - 1., 2.);
    }

    return sum;
}