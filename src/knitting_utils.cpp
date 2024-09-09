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
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, polyscope::SurfaceMesh& psMesh, int direction, FaceData<Vector3>& faceGradients, std::vector<std::pair<int, int>>& edgeMappingsPairs){

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
                Vector3 avgVector = (e1 + e2)/2.;
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


//compute a face-based field through the optimization 
//here the singularities are placed on the vertices
FaceData<Vector3> computeFaceBasedField(VertexPositionGeometry& globalGeometry, Model& model){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    FaceData<Vector3> field(globalMesh);
    globalGeometry.requireFaceAreas();

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
        std::vector<std::array<GRBVar, 3>> W(globalMesh.nFaces());
        for (Face f : globalMesh.faces()){
            std::array<GRBVar, 3> grad = {model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS), model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS),
                                            model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS)};
            
        }

    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }



    return field;
}
