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
    // gradients stored per faces as (3F * 1) 
    // all xs, then all ys, all zs
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
VertexData<Vector2> vertexDirectionField(VertexPositionGeometry& geometry, VertexData<Vector3>& vertexValuedField, 
                                         VertexData<Vector2>& usedRoot){

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
        //storing which root we're using 
        usedRoot[v] = unit(u);
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
//compute a wale one-form 
std::tuple<HalfedgeData<double>, double> computeWaleOneForm(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& model, 
                                        Eigen::SparseMatrix<double, Eigen::RowMajor>& G, std::map<int, int>& vertexMap){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    double period = model.getPeriod();
    std::vector<double> omega = model.getMatchingTerms();
    std::vector<int> faceIndices = model.getFaceIndices();
    std::vector<std::array<double, 3>> gradients = model.getFaceGradients();
    std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = model.getEdgePathConstraints();
    std::vector<int> edgeIndices = model.getEdgeIndices();
    std::vector<int> bdyEdges = model.getBdyEdges();
    std::vector<std::vector<double>> homologyGenerators = model.getHomologyGenerators();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    HalfedgeData<double> oneForm(gluedMesh);
    double objectiveVal = 0;

    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        //don't log output to console
        //env.set(GRB_IntParam_OutputFlag, 0);
        env.start();
        // Create an empty model
        GRBModel model = GRBModel(env);
        //model.set(GRB_DoubleParam_TimeLimit, 5.0);

        //model.set(GRB_IntParam_Method, 2);

        std::vector<GRBVar> generatorIntegers;
        for (size_t i = 0; i < homologyGenerators.size(); i++){
            GRBVar k_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            //add gurobi vars
            generatorIntegers.push_back(k_i);//decision variables
        }

        //defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        //constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
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

        //constraint: alignment constraints in the wale direction
        for (int e : bdyEdges){
            model.addConstr(sigma[gluedMesh.edge(e).halfedge().getIndex()] == 0);
            model.addConstr(sigma[gluedMesh.edge(e).halfedge().twin().getIndex()] == 0);
        }

        //constraint: add integer variables for homology generators
        for (int i = 0; i < homologyGenerators.size(); i++){
            std::vector<double> path = homologyGenerators[i];
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
            //add the constraints for the homology generators
            model.addConstr(pathIntegral == period * generatorIntegers[i]);
        }

        // constraint: 
        // specify singular halfedges and also specify that form values 
        for (Halfedge he : gluedMesh.halfedges()){
            if (edgeIndices[he.edge().getIndex()] != 0){
                model.addConstr(sigma[he.getIndex()] + sigma[he.twin().getIndex()] == edgeIndices[he.edge().getIndex()] * period);
            }
            else{
                model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
            }
        }

        //third constraint 
        //add boundary constraints in the wale direction
        for (int i = 0; i < edgePathConstraints.size(); i++){
            std::vector<double> path = edgePathConstraints[i].first;
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
        GRBQuadExpr obj = 0;
        std::vector<std::array<double, 3>> grads;
        for (Face f : gluedMesh.faces()){
            GRBQuadExpr diffX = (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]) * (gradU[f.getIndex()][0] - gradients[f.getIndex()][0]);
            GRBQuadExpr diffY = (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]) * (gradU[f.getIndex()][1] - gradients[f.getIndex()][1]);
            GRBQuadExpr diffZ = (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]) * (gradU[f.getIndex()][2] - gradients[f.getIndex()][2]);
            double currArea = gluedGeometry.faceAreas[f];
            obj += currArea * (diffX + diffY + diffZ);
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 

        //put the computed one-form into an halfedge vector
        for (Halfedge he : gluedMesh.halfedges()){
            oneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }

        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return std::tie(oneForm, objectiveVal);

}

//compute the function that we'll be using to measure the curl in the wale direction
//this will return the result in the glued setting 
VertexData<double> computeWaleCurlFunction(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                FaceData<Vector3>& courseOneFormGrad, Eigen::SparseMatrix<double, Eigen::RowMajor>& G,
                                                std::map<int, int>& vertexMap, globalBoundaryConditions& globalBdyConditions){

    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh;

    VertexData<double> result(gluedMesh, 0);
    gluedGeometry.requireFaceAreas();

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

        //set the timeout
        //model.getEnv().set(GRB_DoubleParam_TimeLimit, 1.0);
        
        //sigma defined over halfedges
        std::vector<GRBVar> h;
        for (size_t i = 0; i < gluedMesh.nVertices(); i++){
            GRBVar h_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            h.push_back(h_i);//decision variables
        }
        
        std::vector<std::vector<GRBLinExpr>> gradH(gluedMesh.nFaces(), std::vector<GRBLinExpr>(3));
        for (Face f : gluedMesh.faces()){

            GRBLinExpr currGradH = 0.0;
            //X component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex()); it; ++it){
                currGradH += it.value() * h[vertexMap[it.col()]];
            }
            gradH[f.getIndex()][0] = currGradH;
            currGradH = 0.0;
            //Y component
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + gluedMesh.nFaces()); it; ++it){
                currGradH += it.value() * h[vertexMap[it.col()]];
            }
            gradH[f.getIndex()][1] = currGradH;
            currGradH = 0.0;
            //Z component 
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G, f.getIndex() + 2 * gluedMesh.nFaces()); it; ++it){
                currGradH += it.value() * h[vertexMap[it.col()]];
            }
            gradH[f.getIndex()][2] = currGradH;

        }
        
        for (int v : globalBdyConditions.courseStartBoundaryVertices){
            model.addConstr(h[v] == 0);
        }
        for (int v : globalBdyConditions.courseEndBoundaryVertices){
            model.addConstr(h[v] == 1);
        }        

        //set up the objective term
        GRBQuadExpr obj = 0;   

        for (Face f : gluedMesh.faces()){
            GRBQuadExpr diffX = (gradH[f.getIndex()][0] - courseOneFormGrad[f.getIndex()][0]) * (gradH[f.getIndex()][0] - courseOneFormGrad[f.getIndex()][0]);
            GRBQuadExpr diffY = (gradH[f.getIndex()][1] - courseOneFormGrad[f.getIndex()][1]) * (gradH[f.getIndex()][1] - courseOneFormGrad[f.getIndex()][1]);
            GRBQuadExpr diffZ = (gradH[f.getIndex()][2] - courseOneFormGrad[f.getIndex()][2]) * (gradH[f.getIndex()][2] - courseOneFormGrad[f.getIndex()][2]);
            double currArea = gluedGeometry.faceAreas[f];
            obj += currArea * (diffX + diffY + diffZ);
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();

        //put the result in a vertex data
        for (Vertex v : gluedMesh.vertices()){
            result[v] = h[v.getIndex()].get(GRB_DoubleAttr_X);
        } 

    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return result;
    
}

//measure L1 distance between two vector fields 
double L1Distance(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& V1, FaceData<Vector3>& V2){

    SurfaceMesh& mesh = globalGeometry.mesh;
    double L1Distance = 0.;
    for (Face f : mesh.faces()){
        for (int j = 0; j < 3; j++)
            L1Distance += std::fabs(V1[f][j] - V2[f][j]);
    }

    return L1Distance;
}


//@clean 
std::tuple<CornerData<double>, EdgeData<double>> computeWaleStripeInfo(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    std::vector<std::pair<int, int>>& edgeMappingsPairs, std::map<int, int>& edgeMap, 
                                                                    std::map<int, int>& vertexMap, VertexData<double>& timeFunctionGlobal, VertexData<double>& timeFunctionGlued, FaceData<Vector3>& courseOneFormGrad, 
                                                                    Eigen::SparseMatrix<double, Eigen::RowMajor>& G, double period, double knoppelFrequency, globalBoundaryConditions& globalBdyConditions,
                                                                    EdgeData<double>& courseSingularEdgesGlobal, std::map<int, std::vector<Halfedge>> gluedOneRingMap, polyscope::SurfaceMesh& psMesh, 
                                                                    std::vector<std::vector<double>>& allSaddleLoops, std::vector<std::vector<double>>& homologyGenerators){
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    SurfaceMesh& globalMesh = globalGeometry.mesh;

    // Create the Heat Method solver
    HeatMethodDistanceSolver heatSolver(gluedGeometry, 1.0);
    //objectives of the wale solve 
    double oldObj = 0;
    double currObj = 0;
    double oldDistance, newDistance;         
    //which root we used in specifying the direction field
    VertexData<Vector2> usedRoot(globalGeometry.mesh);
    //compute a line field in the tangent space of the vertex
    VertexData<Vector3> vertexVectorField = computeVertexValuedField(globalGeometry, timeFunctionGlobal, PI/2.);
    VertexData<Vector2> lineField = vertexDirectionField(globalGeometry, vertexVectorField, usedRoot);
    EdgeData<double> waleSingularEdgesGlobal(globalGeometry.mesh, 0);
    //doing (1/2.5 * period) just to reduce the number of wale singularities
    VertexData<double> freq(globalGeometry.mesh, 1./(period));
    CornerData<double> stripeValues(globalGeometry.mesh);
    FaceData<int> stripeSingularities(globalGeometry.mesh);
    FaceData<int> fieldSingularities(globalGeometry.mesh);
    //stripe curve network information
    //unique vertices and unique edges
    std::vector<Vector3> uniquePos;
    std::vector<std::array<int, 2>> uniqueEdges;
    //conencted components identified by their component id 
    std::vector<StripeConnectedComponent> components;

    std::tie(stripeValues, stripeSingularities, fieldSingularities) = computeStripePattern(globalGeometry, freq, lineField); // this is a GC call
    
    // Do some visualization
    psMesh.addVertexVectorQuantity("vertexVectorField", vertexVectorField);
    psMesh.addFaceScalarQuantity("knoppel face singularities", stripeSingularities);
    psMesh.addFaceScalarQuantity("knoppel field singularities", fieldSingularities);

    // compute stripe values along integrals
    HalfedgeData<double> formValueHalfedges(globalGeometry.mesh, 0.0);
    EdgeData<double> formValueEdges(globalGeometry.mesh, 0.0);
    //std::vector<std::pair<std::vector<double>, double>> waleBdyEdgePathConstraints; 

    //adjust so that the signs make sense i.e., we consider values on the same sheet
    //THIS IS TO EVALUATE KNÖPPEL INTEGRALS
    // for (Face f : globalGeometry.mesh.faces()){
    //     Vector2 X = Vector2::fromAngle(lineField[f.halfedge().tailVertex()].arg() / 2);
    //     Vector2 root = usedRoot[f.halfedge().tailVertex()];
    //     int sign = dot(X, root) < 0 ? -1 : 1;
    //     for (Halfedge he : f.adjacentHalfedges()){
    //         formValueHalfedges[he] = sign * (stripeValues[he.next().corner()] - stripeValues[he.corner()]); // stripe 1-form
    //         if (he.next() == f.halfedge())
    //             formValueHalfedges[he] += 2 * stripeSingularities[f] * PI;
    //     }
    // }

    // for (Edge e : globalGeometry.mesh.edges()){
    //     formValueEdges[e] = formValueHalfedges[e.halfedge()];
    // }

    //for a model with n boundaries, we need (n - 1) bdy integral constraints
    // for (int i = 0; i < globalBdyConditions.waleBdyPathConstraints.size(); i++){
    //     std::vector<double> path = globalBdyConditions.waleBdyPathConstraints[i];
    //     double sum = 0;
    //     int k = 0;
    //     for (int j = 0; j < path.size(); j++){
    //         if (path[j] > 0){
    //             Edge e = gluedGeometry.mesh.edge(j);
    //             sum += formValueEdges[globalGeometry.mesh.edge(j)];
    //         }
    //         if (path[j] < 0){
    //             Edge e = gluedGeometry.mesh.edge(j);
    //             sum += -1.0 * formValueEdges[globalGeometry.mesh.edge(j)];
    //         }
    //     }
    //     // std::cout << "sum / (2 * PI) = " << sum / (2. * PI) << std::endl;
    //     // std::cout << "number of stripes on boundary " << i << " = " << std::round(sum/ (2. * PI)) << std::endl;
    //     k = std::round(sum / (2. * PI));
    //     waleBdyEdgePathConstraints.push_back(std::make_pair(path, k));

    // }

    std::vector<Vector3> knoppelPos; 
    std::vector<std::array<size_t, 2>> knoppelEdges; 
    std::tie(knoppelPos, knoppelEdges) = extractPolylinesFromStripePattern(globalGeometry, stripeValues, stripeSingularities,
                                            fieldSingularities, lineField, false);
    auto knoppelStripes = polyscope::registerCurveNetwork("knoppel wale stripes stripes", knoppelPos, knoppelEdges);
    knoppelStripes -> setRadius(0.001);
    knoppelStripes -> setEnabled(false);

    //average the final course 1-form gradient onto edges
    EdgeData<double> omegaWaleGlobal = computeMatchingOneForm(globalGeometry, 1, courseOneFormGrad, edgeMappingsPairs);
    EdgeData<double> omegaWaleGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, omegaWaleGlobal, edgeMap);
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

    //before rotating the course face gradients compute the "psuedo time function" we'll use for curl in the wale direction 
    //VertexData<double> pseudoTimeFunctionGlued = computeWaleCurlFunction(globalGeometry, gluedGeometry, courseOneFormGrad,
    //                                                                        G, vertexMap, globalBdyConditions);
    //VertexData<double> pseudoTimeFunctionGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, pseudoTimeFunctionGlued, vertexMap);
    psMesh.addVertexScalarQuantity("wale time function", timeFunctionGlobal);
    //FaceData<Vector3> waleCurlFunctionGrad = computeTimeFunctionFaceGrad(globalGeometry, pseudoTimeFunctionGlobal);
    FaceData<Vector3> waleCurlFunctionGrad = computeTimeFunctionFaceGrad(globalGeometry, timeFunctionGlobal);
    psMesh.addFaceVectorQuantity("wale time function grad", waleCurlFunctionGrad);    


    Eigen::Map<Eigen::VectorXi> faceIndicesWaleEig(stripeSingularities.raw().data(), (gluedGeometry.mesh).nFaces());
    //place wale singularities at edges 
    std::vector<int> faceIndicesWaleModel(faceIndicesWaleEig.data(), faceIndicesWaleEig.data() + faceIndicesWaleEig.rows());
    std::vector<std::array<double, 3>> modelFaceGradients;
    //std::vector<std::pair<int, int>> singularEdges;
    //singular edges in the gurobi optimization (in a better representation)
    std::vector<int> edgeIndices(gluedGeometry.mesh.nEdges(), 0);

    globalGeometry.requireFaceNormals();
    for (Face f : globalGeometry.mesh.faces()){ 
        courseOneFormGrad[f] = courseOneFormGrad[f].normalize();
        //rotate the final course gradients 
        courseOneFormGrad[f] = courseOneFormGrad[f].rotateAround(globalGeometry.faceNormals[f], PI/2.);
        //also rotate the wale curl function gradient
        waleCurlFunctionGrad[f] = waleCurlFunctionGrad[f].normalize();
        waleCurlFunctionGrad[f] = waleCurlFunctionGrad[f].rotateAround(globalGeometry.faceNormals[f], PI/2.);

        modelFaceGradients.push_back(std::array{waleCurlFunctionGrad[f][0], waleCurlFunctionGrad[f][1], 
                                                    waleCurlFunctionGrad[f][2]});
        //place wale singularities at edges
        // if (stripeSingularities[f] != 0){
        //     bool toSkip = false;
        //     for (Edge e : f.adjacentEdges()){
        //         //skip boundary edges and edges on faces that already have a singularity in the course or wale direction
        //         if (e.isBoundary() || waleSingularEdgesGlobal[e] != 0 || courseSingularEdgesGlobal[e] != 0) toSkip = true;
        //     }
        //     if (toSkip) continue;
        //     int edge = findSingularEdgeFromSingularFace(globalGeometry, f.getIndex(), courseOneFormGrad[f], 0., timeFunctionGlobal, 0.);
        //     if (edge == -1) continue;//couldn't find a very well-aligned edge
        //     //reduce the singular index to +/- 1 (cause of Autoknit constraints)
        //     //int index = stripeSingularities[f] > 0 ? 1 : -1;
        //     //if we don't want to reduce the index constraints to +/-1 
        //     int index = stripeSingularities[f];
        //     singularEdges.push_back(std::make_pair(edgeMap[edge], index));
        //     waleSingularEdgesGlobal[edge] = index;
            
        // }
    }

    psMesh.addFaceVectorQuantity("rotated wale time function grad", modelFaceGradients);
    
    //convert the global singular edges to glued singular in the course direction 
    EdgeData<double> courseSingularEdgesGlued = convertGlobalToGluedEdgeFunction(globalGeometry, gluedGeometry, courseSingularEdgesGlobal, edgeMap);
    std::vector<int> gluedEdgeSingularities(gluedGeometry.mesh.nEdges(), 0);
    for (Edge e : gluedGeometry.mesh.edges()){
        if (std::fabs(courseSingularEdgesGlued[e]) > 0)
            gluedEdgeSingularities[e.getIndex()] = 1;
    }

    std::vector<Vertex> heatSourceVerts;
    heatSourceVerts = getSaddleVertices(gluedGeometry, timeFunctionGlued);
    //keep singualrities away from all saddle loops
    for (int i = 0; i < allSaddleLoops.size(); i++){
        std::vector<double> path = allSaddleLoops[i];
        for (int j = 0; j < path.size(); j++){
            if (std::fabs(path[j]) > 0){
                Edge e = gluedGeometry.mesh.edge(edgeMap[j]);
                heatSourceVerts.push_back(e.halfedge().tailVertex());
                heatSourceVerts.push_back(e.halfedge().tipVertex());
            }
        }
    }

    for (Edge e : globalGeometry.mesh.edges()){
        if (std::fabs(courseSingularEdgesGlobal[e]) > 0){
            heatSourceVerts.push_back(gluedGeometry.mesh.vertex(vertexMap[e.halfedge().tailVertex().getIndex()]));
            heatSourceVerts.push_back(gluedGeometry.mesh.vertex(vertexMap[e.halfedge().tipVertex().getIndex()]));
        }
    }


    //add wale alignment constraints 
    // EdgeData<double> constrainedEdges(globalGeometry.mesh, 0.0);
    // double alignment = 0.6;
    // std::vector<int> waleBdyEdges;
    // for (Edge e : globalGeometry.mesh.edges()){
    //     if (!e.isBoundary()) continue;
    //     Vector3 edgeVector = globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()];
    //     edgeVector = edgeVector.normalize();
    //     if (std::fabs(dot(edgeVector, waleCurlFunctionGrad[e.halfedge().face()])) < alignment){
    //         waleBdyEdges.push_back(edgeMap[e.getIndex()]);
    //         constrainedEdges[e] = 1;
    //     }
    // }
    // psMesh.addEdgeScalarQuantity("constrained edges", constrainedEdges);

    VertexData<double> waleVertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, waleCurlFunctionGrad, gluedOneRingMap,
                                                         gluedEdgeSingularities, heatSolver, vertexMap);
    VertexData<double> allDist(gluedGeometry.mesh, 1.0);
    VertexData<double> allDistGlobal(globalGeometry.mesh);
    VertexData<double> waleWeighting(globalGeometry.mesh);
    if (heatSourceVerts.size()){
        allDist = heatSolver.computeDistance(heatSourceVerts);
        double maxVal = std::numeric_limits<double>::min();
        double maxSourceVal = std::numeric_limits<double>::min();
        //for all the source vertices, find the max value 
        for (Vertex v : heatSourceVerts){
            maxSourceVal = std::max(maxSourceVal, allDist[v]);
        }
        //find the maximum value over all the distances
        for (Vertex v : gluedGeometry.mesh.vertices()){
            maxVal = std::max(maxVal, allDist[v]);
        }
        //shift down all the source vertex values
        for (Vertex v : heatSourceVerts){
            allDist[v] -= maxSourceVal;
        }
        //clip all values to 0
        for (Vertex v : gluedGeometry.mesh.vertices()){
            allDist[v] = std::max(allDist[v], 0.0);
        }
        allDistGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, allDist, vertexMap);
        for (Vertex v : globalGeometry.mesh.vertices()){
            waleWeighting[v] = (allDistGlobal[v] > 2*period);
            //waleWeighting[v] = (1 - exp(-pow(allDistGlobal[v], 2.) / (2 * pow(2*period, 2))));
            waleVertexCurl[v] = waleWeighting[v] * waleVertexCurl[v];
        }
    }

    psMesh.addVertexScalarQuantity("initial wale weighting", waleWeighting);
    psMesh.addVertexScalarQuantity("initial all distance wale", allDistGlobal);
    EdgeData<double> waleEdgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, waleVertexCurl);
    psMesh.addVertexScalarQuantity("initial wale vertex curl", waleVertexCurl);
    psMesh.addEdgeScalarQuantity("wale edge curl without any wale singularities", waleEdgeCurl);


    double oldL1Distance = 0.;
    double newL1Distance = 0.;
    Model modelWale; 
    modelWale.setPeriod(period);
    modelWale.setFaceGradients(modelFaceGradients);
    //modelWale.setEdgeIndices(edgeIndices);
    HalfedgeData<double> sigmaWaleGlued(gluedGeometry.mesh);
    //solve the model with out any singularities
    // std::tie(sigmaWaleGlued, currObj) = computeWaleOneForm(globalGeometry, gluedGeometry, modelWale, G, vertexMap);
    // FaceData<Vector3> gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, sigmaWaleGlued);
    //newL1Distance = L1Distance(globalGeometry, gradSigmaTilde, waleCurlFunctionGrad);
    //newDistance = computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde);
    //std::cout << "distance from unit norm without any wale singularities = " << newDistance << std::endl;

    // std::cout << "objective without any wale singularities = " << currObj << std::endl;
    // psMesh.addFaceVectorQuantity("wale gradient after without any wale singularities", gradSigmaTilde);
    //oldDistance = newDistance;
    //oldObj = currObj;
    //oldL1Distance = newL1Distance;

    CornerData<double> stripeValuesOneFormGlued;
    FaceData<int> stripeIndicesOneFormGlued;
    // //std::tie(stripeValuesOneFormGlued, stripeIndicesOneFormGlued) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, sigmaWaleGlued, period);

    // int numWaleSingularities = 0;
    // // int maxWaleSingularities = 5;
    // int topPairs = 10;

    // int numPos = 0;
    // int numNeg = 0;

    // EdgeData<double> checked(globalGeometry.mesh, 0.0);

    /** 
    bool toBreak = false;
    while(true){
        std::vector<std::pair<int, int>> waleSingularityEdges = findWaleSingularityEdges(globalGeometry, gluedGeometry,
                                                                waleEdgeCurl, topPairs);
        
        std::cout << "size of waleSingularityEdges = " << waleSingularityEdges.size() << std::endl;
        
        //DEBUG
        //print out the candidate edges
        EdgeData<int> candidateEdges(globalGeometry.mesh, 0);
        std::vector<Vector3> pointCloud;
        for (int i = 0; i < waleSingularityEdges.size(); i++){
            candidateEdges[waleSingularityEdges[i].first] = waleSingularityEdges[i].second;
            pointCloud.push_back(globalGeometry.vertexPositions[gluedGeometry.mesh.edge(waleSingularityEdges[i].first).halfedge().tailVertex()]);
            pointCloud.push_back(globalGeometry.vertexPositions[gluedGeometry.mesh.edge(waleSingularityEdges[i].first).halfedge().tipVertex()]);
        }
        polyscope::registerPointCloud("top edges after " + std::to_string(numWaleSingularities) + "wale singularities", pointCloud)->setEnabled(false);
        psMesh.addEdgeScalarQuantity("candidate edges after " + std::to_string(numWaleSingularities) + "wale singularities", candidateEdges);
        int numSkips = 0;
        int maxSkips = topPairs;
        for (auto s : waleSingularityEdges){
            if (s.second > 0){
                numPos++;
            }
            else{
                numNeg++;
            }
            //test model 
            //to test a new pair of singularities
            Model testModel = modelWale; 
            std::vector<int> testEdgeIndices = edgeIndices;
            testEdgeIndices[edgeMap[s.first]] = -1 * s.second;
            testModel.setEdgeIndices(testEdgeIndices);
            std::tie(sigmaWaleGlued, currObj) = computeWaleOneForm(globalGeometry, gluedGeometry, testModel, G, vertexMap);
            gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, sigmaWaleGlued);
            checked[s.first] = 1.0;
            FaceData<Vector3> diffFromTarget(gluedGeometry.mesh);
            for (Face f : gluedGeometry.mesh.faces())
                for (int i = 0; i < 3; i++)
                    diffFromTarget[f][i] = gradSigmaTilde[f][i] - modelFaceGradients[f.getIndex()][i];
            psMesh.addFaceVectorQuantity("diffFromTarget", diffFromTarget);
            psMesh.addFaceVectorQuantity("candidate gradSigmaTilde", gradSigmaTilde);

            //newDistance = computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde);            
            if (currObj >= oldObj){
                std::cout << "new objective = " << currObj << std::endl;
                std::cout << "old objective = " << oldObj << std::endl;
                std::cout << "norm is not improving trying next singular edge " << std::endl;
                numSkips++;
                std::cout << "size of waleSingularityEdges = " << waleSingularityEdges.size() << std::endl;
                std::cout << "numSkips = " << numSkips << std::endl;
                if (numSkips >= maxSkips){
                    std::cout << "breaking after " << numWaleSingularities << " singularity insertions " << std::endl;
                    std::cout << "new objective = " << currObj << std::endl;
                    std::cout << "old objective = " << oldObj << std::endl;
                    toBreak = true;
                    break;
                }
                continue;
            }
            else{//we have a valid pair of singularities that improves the objective
                numWaleSingularities++;
                std::cout << "not breaking after " << std::to_string(numWaleSingularities) << " singularity insertions " << std::endl;
                std::cout << "new objective = " << currObj << std::endl;
                std::cout << "old objective = " << oldObj << std::endl;
                //oldDistance = newDistance;
                oldObj = currObj;
                //oldL1Distance = newL1Distance;
                if (s.second > 0){
                    numPos++;
                    std::cout << "INDEX IS POSITIVE " << std::endl;
                }
                else{
                    numNeg++;
                    std::cout << "INDEX IS NEGATIVE " << std::endl;
                }
                edgeIndices[edgeMap[s.first]] = -1 * s.second;
                waleSingularEdgesGlobal[s.first] = s.second;
                psMesh.addEdgeScalarQuantity("wale singular edge after " + std::to_string(numWaleSingularities), waleSingularEdgesGlobal);
                gluedEdgeSingularities[edgeMap[s.first]] = 1;
                //set the singular edges 
                modelWale.setEdgeIndices(edgeIndices);
                std::tie(stripeValuesOneFormGlued, stripeIndicesOneFormGlued) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, sigmaWaleGlued, period);


                // Plot stripes
                std::vector<Vector3> positionsWale;
                std::vector<std::array<int, 2>> edgesWale;
                std::tie(positionsWale, edgesWale) = generateIsoLines(globalGeometry, stripeValuesOneFormGlued, stripeIndicesOneFormGlued, period);
                auto waleStripes = polyscope::registerCurveNetwork("wale stripes after " + std::to_string(numWaleSingularities) + "wale singularities", positionsWale, edgesWale);
                waleStripes -> setRadius(0.001);
                waleStripes -> setEnabled(false);

                //compute the impulse function 
                HalfedgeData<double> waleVirtualSigma = computeWaleVirtualSigma(globalGeometry, gluedGeometry, modelWale);
                //view the Green's function 
                FaceData<Vector3> impulseGrad = computeOneFormFaceGrad(globalGeometry, gluedGeometry, waleVirtualSigma);
                psMesh.addFaceVectorQuantity("Green's function after " + std::to_string(numWaleSingularities), impulseGrad);
                //subtracts off the impulse function
                sigmaWaleGlued = sigmaWaleGlued - waleVirtualSigma;
                FaceData<Vector3> adjustedGradSigmaWaleGlued = computeOneFormFaceGrad(globalGeometry, gluedGeometry, sigmaWaleGlued);
                waleVertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, adjustedGradSigmaWaleGlued, gluedOneRingMap,
                                                    gluedEdgeSingularities, heatSolver, vertexMap);
                

                FaceData<double> stoppingCondn(globalGeometry.mesh, 0);
                for (Face f : globalGeometry.mesh.faces()){
                    stoppingCondn[f] = pow(norm(adjustedGradSigmaWaleGlued[f] - waleCurlFunctionGrad[f]), 2);
                }
                psMesh.addFaceScalarQuantity("stopping condn after " + std::to_string(numWaleSingularities), stoppingCondn);

                //keep wale singularities away from each other
                heatSourceVerts.push_back(gluedGeometry.mesh.edge(edgeMap[s.first]).halfedge().tailVertex());
                heatSourceVerts.push_back(gluedGeometry.mesh.edge(edgeMap[s.first]).halfedge().tipVertex());
                allDist = heatSolver.computeDistance(heatSourceVerts);
                
                double maxVal = std::numeric_limits<double>::min();
                double maxSourceVal = std::numeric_limits<double>::min();
                //for all the source vertices, find the max value 
                for (Vertex v : heatSourceVerts){
                    maxSourceVal = std::max(maxSourceVal, allDist[v]);
                }
                //find the maximum value over all the distances
                for (Vertex v : gluedGeometry.mesh.vertices()){
                    maxVal = std::max(maxVal, allDist[v]);
                }
                //shift down all the source vertex values
                for (Vertex v : heatSourceVerts){
                    allDist[v] -= maxSourceVal;
                }
                //clip all values to 0
                for (Vertex v : gluedGeometry.mesh.vertices()){
                    allDist[v] = std::max(allDist[v], 0.0);
                }
                VertexData<double> allDistGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, allDist, vertexMap);
                VertexData<double> gaussianMask(globalGeometry.mesh);
                //use a hard mask in the wale direction
                for (Vertex v : globalGeometry.mesh.vertices()) {
                    gaussianMask[v] = (allDistGlobal[v] > 2*period);
                    //gaussianMask[v] = 1 - exp(-pow(allDistGlobal[v], 2.) / (2 * pow(2*period, 2)));
                    waleVertexCurl[v] = gaussianMask[v] * waleVertexCurl[v];
                }
                waleEdgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, waleVertexCurl);
                for (Edge e : globalGeometry.mesh.edges()){
                    if (std::fabs(checked[e]) > 0) waleEdgeCurl[e] = 0.;
                }
                psMesh.addFaceVectorQuantity("wale gradient after " + std::to_string(numWaleSingularities) + "wale singularities ", adjustedGradSigmaWaleGlued);
                psMesh.addVertexScalarQuantity("gaussian mask after " + std::to_string(numWaleSingularities) + " wale singularities", gaussianMask);
                psMesh.addEdgeScalarQuantity("wale edge curl after " + std::to_string(numWaleSingularities) + " wale singularities", waleEdgeCurl);
                break;
            }
        }
        if (toBreak) break;//we are no longer improving the distance to unit norm 
    }
    */

    //-------------------------TESTING----------------------//
    //Attempting to find the center of distributions using Vector Heat Method
    VertexData<double> curlMeasure = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                    waleCurlFunctionGrad, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);

    psMesh.addVertexScalarQuantity("initial curl measure ", curlMeasure);

    // Cap curl measure to avoid high concentration of singularities
    for (Vertex v : globalMesh.vertices()) {
        curlMeasure[v] = fmin(curlMeasure[v], period / (3*globalGeometry.vertexDualAreas[v]));
        curlMeasure[v] = fmax(curlMeasure[v], -period / (3*globalGeometry.vertexDualAreas[v]));
    }
    psMesh.addVertexScalarQuantity("initial curl measure (capped)", curlMeasure);
    //now divide the measure into positive and negative 
    VertexData<double> posMeasure(globalMesh, 0.0);
    VertexData<double> negMeasure(globalMesh, 0.0);
    double totalPosMeasure = 0;
    double totalNegMeasure = 0;
    for (Vertex v : globalMesh.vertices()){
        if (curlMeasure[v] > 0){
            posMeasure[v] = curlMeasure[v];
            totalPosMeasure += globalGeometry.vertexDualAreas[v] * posMeasure[v];
        }
        else{
            negMeasure[v] = std::fabs(curlMeasure[v]);
            totalNegMeasure += globalGeometry.vertexDualAreas[v] * negMeasure[v];
        }
    }
    H(totalPosMeasure);
    H(totalNegMeasure);

    psMesh.addVertexScalarQuantity("initial positive measure ", posMeasure);
    psMesh.addVertexScalarQuantity("initial negative measure ", negMeasure);
    std::unique_ptr<ManifoldSurfaceMesh> manifoldGlobalMesh = globalMesh.toManifoldMesh();
    int numSings = 0.;
    VoronoiOptions posOptions = defaultVoronoiOptions;
    posOptions.nSites = std::round(totalPosMeasure / period);
    posOptions.useDelaunay = false;
    posOptions.computeDistributions = true;
    posOptions.iterations = 500;
    VoronoiResult posVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(*manifoldGlobalMesh, globalGeometry, posOptions, posMeasure, psMesh);
    std::vector<Vector3> positiveCenters;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
       positiveCenters.push_back(posVoronoiCenters.siteLocations[i].interpolate(globalGeometry.vertexPositions));
    }
    polyscope::registerPointCloud("positive voronoi sites (wale)", positiveCenters);
    std::vector<VertexData<double>> posSiteDistributions = posVoronoiCenters.siteDistributions;

    //print the masses
    for (size_t i = 0; i < posSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){
            if (std::fabs(posSiteDistributions[i][v]) > 1e-8){
                mass += posSiteDistributions[i][v];
            }
        }
        std::cout << "positive mass at site " << i << " = " << mass << std::endl;
        psMesh.addVertexScalarQuantity("positve site distribution (wale)" + std::to_string(i), posSiteDistributions[i]);
    }

   
    VoronoiOptions negOptions = defaultVoronoiOptions;
    negOptions.nSites = std::round(totalNegMeasure / period); // we don't care if this is different that posOptions.nSites
    negOptions.useDelaunay = false;
    negOptions.computeDistributions = true;
    negOptions.iterations = 500;
    std::cout << "# positive sites " << posOptions.nSites << std::endl;
    std::cout << "# negative sites " << negOptions.nSites << std::endl;
    VoronoiResult negVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(*manifoldGlobalMesh, globalGeometry, negOptions, negMeasure, psMesh);
    std::vector<Vector3> negativeCenters;
    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
       negativeCenters.push_back(negVoronoiCenters.siteLocations[i].interpolate(globalGeometry.vertexPositions));
    }
    polyscope::registerPointCloud("negative voronoi sites (wale)", negativeCenters);
    std::vector<VertexData<double>> negSiteDistributions = negVoronoiCenters.siteDistributions;

    //print the masses
    for (size_t i = 0; i < negSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){
            if (std::fabs(negSiteDistributions[i][v]) > 1e-8){
                mass += negSiteDistributions[i][v];
            }
        }
        std::cout << "negative mass at site " << i << " = " << mass << std::endl;
        psMesh.addVertexScalarQuantity("negative site distribution (wale)" + std::to_string(i), negSiteDistributions[i]);
    }

    //store a vector of vertex ids and singularity indices (for optimal matching)
    std::vector<std::pair<Vertex, int>> singularities;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
        SurfacePoint facePoint = posVoronoiCenters.siteLocations[i].inSomeFace();
        Face f = facePoint.face;
        Edge singularEdge;
        double maxDotProd = -DBL_MAX;
        // if (facePoint.faceCoords.x == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().vertex(), 1));
        // else if (facePoint.faceCoords.y == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().vertex(), 1));
        // else if (facePoint.faceCoords.z == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().next().vertex() , 1));
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - 
                                globalGeometry.vertexPositions[he.tailVertex()]).normalize();
            if (std::fabs(dot(heVec, waleCurlFunctionGrad[f])) > maxDotProd){
                maxDotProd = std::fabs(dot(heVec, waleCurlFunctionGrad[he.face()]));
                singularEdge = he.edge();
            }
        }
        //increment the indices really close so indices cancel out
        edgeIndices[singularEdge.getIndex()] = -1.0;
        waleSingularEdgesGlobal[singularEdge] = 1.0;
    }

    // for (int i = 0; i < posVoronoiCenters.steps.size(); i++){
    //     std::vector<SurfacePoint> steps = posVoronoiCenters.steps[i]; //steps for this site
    //     std::vector<Vector3> stepPos;
    //     for (int i = 0; i < steps.size(); i++){
    //         stepPos.push_back(steps[i].interpolate(globalGeometry.vertexPositions));
    //     }
    //     polyscope::registerPointCloud("pos site " + std::to_string(i) + " steps ", stepPos);
    // }

    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
        SurfacePoint facePoint = negVoronoiCenters.siteLocations[i].inSomeFace();
        Face f = facePoint.face;
        Edge singularEdge;
        double maxDotProd = -DBL_MAX;
        // if (facePoint.faceCoords.x == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().vertex(), -1));
        // else if (facePoint.faceCoords.y == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().vertex(), -1));
        // else if (facePoint.faceCoords.z == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().next().vertex() , -1));
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - 
                                globalGeometry.vertexPositions[he.tailVertex()]).normalize();
            if (std::fabs(dot(heVec, waleCurlFunctionGrad[f])) > maxDotProd){
                maxDotProd = std::fabs(dot(heVec, waleCurlFunctionGrad[he.face()]));
                singularEdge = he.edge();
            }
        }
        //increment the indices really close so indices cancel out
        edgeIndices[singularEdge.getIndex()] = 1.0;
        waleSingularEdgesGlobal[singularEdge] = -1.0;
    }

    //set the edge indices from the singularities found using the OT method 
    //select the outgoing halfedge that is most aligned with the gradient of the time function
    // for (auto p : singularities){
    //     //find the halfedge that is most aligned with the gradient of the time function
    //     Edge singularEdge;
    //     double maxDotProd = DBL_MIN;
    //     for (Halfedge he : p.first.outgoingHalfedges()){
    //         Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - globalGeometry.vertexPositions[he.tailVertex()]).normalize();
    //         if (std::fabs(dot(heVec, waleCurlFunctionGrad[he.face()])) > maxDotProd){
    //             maxDotProd = std::fabs(dot(heVec, waleCurlFunctionGrad[he.face()]));
    //             singularEdge = he.edge();
    //         }
    //     }
    //     //increment the indices really close so indices cancel out
    //     edgeIndices[singularEdge.getIndex()] += (-1.0 * p.second);
    //     waleSingularEdgesGlobal[singularEdge] = p.second;
    // }


    modelWale.setEdgeIndices(edgeIndices);
    // Compute the final non-integer 1-form
    std::tie(sigmaWaleGlued, currObj) = computeWaleOneForm(globalGeometry, gluedGeometry, modelWale, G, vertexMap);
    
    std::tie(stripeValuesOneFormGlued, stripeIndicesOneFormGlued) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, sigmaWaleGlued, period);


    // Compute (integer) sum of indices
    int sumIndices = 0;
    for (int i = 0; i < gluedGeometry.mesh.nEdges(); i++)
        sumIndices += edgeIndices[i]; // I think the sign convention in edgeIndices is wrong, which is why we have to put a negative here

    // Compute (non-integer) number of stripes on boundaries
    std::vector<double> boundaryStripes;
    for (std::vector<double> path : globalBdyConditions.waleBdyPathConstraints) {
        double integral = 0;
        for (int i = 0; i < path.size(); i++) if (path[i] != 0) {
            Edge e = gluedGeometry.mesh.edge(i);
            integral += sigmaWaleGlued[e.halfedge()] * path[i];
        }
        boundaryStripes.push_back(integral / period);
    }

    // Round it up
    std::vector<int> boundaryStripesInt = roundWithSum(boundaryStripes, sumIndices);

    // Set up the edge path constraints 
    std::vector<std::pair<std::vector<double>, double>> waleEdgePathConstraints;
    for (int i = 0; i < globalBdyConditions.waleBdyPathConstraints.size(); i++){
        std::vector<double> path = globalBdyConditions.waleBdyPathConstraints[i];
        waleEdgePathConstraints.push_back(std::make_pair(path, -boundaryStripesInt[i]));
    }

    // solve the model with all the path constraints and boundary constraints 
    //modelWale.setBdyEdges(waleBdyEdges);
    modelWale.setEdgeIndices(edgeIndices);
    modelWale.setEdgePathConstraints(waleEdgePathConstraints);
    modelWale.setHomologyGenerators(homologyGenerators);
    std::tie(sigmaWaleGlued, currObj) = computeWaleOneForm(globalGeometry, gluedGeometry, modelWale, G, vertexMap);
    std::tie(stripeValuesOneFormGlued, stripeIndicesOneFormGlued) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, sigmaWaleGlued, period);
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesOneFormGlued, stripeIndicesOneFormGlued, period, edgeMap);
    auto waleStripes = polyscope::registerCurveNetwork("our wale stripes ", uniquePos, uniqueEdges);
    waleStripes -> setRadius(0.001);
    waleStripes -> setEnabled(false);

    std::cout << "boundary sum = " << sumIndices << std::endl;
    // std::cout << "Number of positive edges sampled = " << numPos << std::endl;
    // std::cout << "Number of negative edges sampled = " << numNeg << std::endl;


    return std::tie(stripeValuesOneFormGlued, waleSingularEdgesGlobal);

    //---------------------------TESTING END------------------------------------//

}

HalfedgeData<double> computeWaleVirtualSigma(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model &model){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);
    double period = model.getPeriod();
    std::vector<int> edgeIndices = model.getEdgeIndices();
    //std::vector<std::pair<std::vector<double>, double>> edgePathConstraints = model.getEdgePathConstraints();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();
    //require halfedge cotan weights 
    gluedGeometry.requireHalfedgeCotanWeights();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    
    try {
        //reformulate the problem in terms of halfedges 
        // Create an environment
        GRBEnv env = GRBEnv(true);
        //don't log output to console 
        //env.set(GRB_IntParam_OutputFlag, 0);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        //model.getEnv().set(GRB_DoubleParam_TimeLimit, 1.0);
        
        //sigma defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++){
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            //add gurobi vars
            sigma.push_back(sigma_i);//decision variables
        }

        // //add the second constraint while we're here
        // //constraint - (d1*sigma) == kP, P is period of optimization, k \in \mathbb{Z}
        // //here it will (d1 * sigma) == 0
        for (Face f : gluedMesh.faces()){
            GRBLinExpr lhs = 0.0;
            Halfedge hij = f.halfedge();
            Halfedge hjk = hij.next();
            Halfedge hki = hjk.next();
            lhs = sigma[hij.getIndex()] + sigma[hjk.getIndex()] + sigma[hki.getIndex()];
            model.addConstr(lhs == 0);
        }

        // constraint: 
        // specify singular halfedges and also specify that form values 
        for (Halfedge he : gluedMesh.halfedges()){
            if (edgeIndices[he.edge().getIndex()] != 0){
                model.addConstr(sigma[he.getIndex()] + sigma[he.twin().getIndex()] == edgeIndices[he.edge().getIndex()] * period);
            }
            else{
                model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
            }
        }

        // //add constraint 
        // //add boundary constraints in the wale direction
        // for (int i = 0; i < edgePathConstraints.size(); i++){
        //     std::vector<double> path = edgePathConstraints[i].first;
        //     GRBLinExpr pathIntegral = 0;
        //     std::vector<double> hePath(gluedMesh.nHalfedges(), 0.0);
        //     for (int j = 0; j < gluedMesh.nEdges(); j++){
        //         if (path[j] > 0){
        //             hePath[gluedMesh.edge(j).halfedge().getIndex()] = path[j];
        //         }
        //         else if (path[j] < 0){
        //             hePath[gluedMesh.edge(j).halfedge().twin().getIndex()] = path[j];
        //         }
        //     }
        //     for (int k = 0; k < gluedMesh.nHalfedges(); k++){
        //         pathIntegral += hePath[k] * sigma[k];
        //     }
        //     model.addConstr(pathIntegral == period * edgePathConstraints[i].second);
        // }

        //set up the objective term
        GRBQuadExpr obj = 0;        
    
        //setting the objective to be min ||\sigma||^2
        for (Halfedge he : gluedMesh.halfedges()){
            if (gluedGeometry.edgeCotanWeights[he.edge()] > 0)
                obj += gluedGeometry.edgeCotanWeights[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
            // obj += gluedGeometry.edgeLengths[he.edge()] * sigma[he.getIndex()] * sigma[he.getIndex()];
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
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
    return gluedOneForm;
}

//find singular edges in the wale direction in the global setting
//first entry is the edge index, second entry is the index of the singularity
std::vector<std::pair<int, int>> findWaleSingularityEdges(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry,
                                          EdgeData<double>& edgeCurl, int topPairs){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    std::vector<std::pair<int, int>> singularEdges;
    
    int numEdges = 0;
    //stores in descending order
    std::map<double, int, std::greater<double>> curlToEdgeIndex;
    for (Edge e : globalMesh.edges()){
        curlToEdgeIndex[std::fabs(edgeCurl[e])] = e.getIndex();
    }

    
    // for (auto &entry : curlToEdgeIndex){
    //     if (globalGeometry.mesh.edge(entry.second).isBoundary()) continue;//don't place singularities on boundary edges
    //     if (edgeCurl[entry.second] > 0) singularEdges.push_back(std::make_pair(entry.second, 1));
    //     else singularEdges.push_back(std::make_pair(entry.second, -1));
    //     numEdges++;
    //     if (numEdges == topPairs) break;
    // }

    // NEW: we force to have half positive, half negative candidate edges
    int nPos = 0, nNeg = 0; // number of positive and negative edge candidates
    for (auto &entry : curlToEdgeIndex){
        if (globalGeometry.mesh.edge(entry.second).isBoundary()) continue;//don't place singularities on boundary edges
        if (nPos < topPairs/2 && edgeCurl[entry.second] > 0) {
            singularEdges.push_back(std::make_pair(entry.second, 1));
            nPos++;
            numEdges++;
        } 
        if (nNeg < topPairs/2 && edgeCurl[entry.second] < 0) {
            singularEdges.push_back(std::make_pair(entry.second, -1));
            nNeg++;
            numEdges++;
        }
        if (numEdges == topPairs) break;
    }

    H(singularEdges.size());

    // // for (auto &s : singularEdges){
    // //     std::cout << "singular edges = " << s.first << " " << s.second << std::endl;
    // // }

    // return singularEdges;

    // for (auto &s : singularEdges){
    //     std::cout << "singular edges = " << s.first << " " << s.second << std::endl;
    // }

    return singularEdges;

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

//find edge singularity pairs 
//this method finds the max curl edges and finds an edge on the same isoline of the time function 
//of the opposite sign
std::vector<std::pair<int, int>> findEdgeSingularityPairsUsingMaxCurls(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                            Eigen::MatrixXd& V, Eigen::MatrixXi& F, EdgeData<double>& curl,
                                            VertexData<double>& globalTimeFunction, std::map<int, int>& hashedUsedIsoVals, 
                                            int numSingularityPairs){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    std::vector<std::pair<int, int>> singEdgePairs;
    
    //stores in ascending order
    std::map<double, int> curlToEdgeIndex;
    for (Edge e : globalMesh.edges()){
        curlToEdgeIndex[curl[e]] = e.getIndex();
    }

    int ctr = 0;
    for (auto entry : curlToEdgeIndex){
        Edge minCurlEdge = globalMesh.edge(entry.second);
        double minCurlIsoVal = 0.5 * (globalTimeFunction[minCurlEdge.halfedge().tailVertex()] + 
                                    globalTimeFunction[minCurlEdge.halfedge().tipVertex()]);
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, minCurlIsoVal);
        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        if (facesPerComponent.size() == 1){//there's only one connected component in this isoline
            std::vector<int> facesInComponent = facesPerComponent[0];
            double maxCurl = -DBL_MAX;
            int maxEdge;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                for (Edge e : currFace.adjacentEdges()){
                    //check if current isoline crosses this edge
                    if ((minCurlIsoVal > globalTimeFunction[e.halfedge().tailVertex()] && minCurlIsoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    || (minCurlIsoVal > globalTimeFunction[e.halfedge().tipVertex()] && minCurlIsoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] > maxCurl){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_pair(maxEdge, minCurlEdge.getIndex()));
        }
        else{//need to find only the faces in this particular connected component
            std::vector<int> facesInComponent;
            bool shouldBreak = false;
            for (int i = 0; i < facesPerComponent.size(); i++){
                std::vector<int> currComponent = facesPerComponent[i];
                for (int j = 0; j < currComponent.size(); j++){
                    Face currFace = globalMesh.face(currComponent[j]);
                    for (Edge e : currFace.adjacentEdges()){
                        if (e == minCurlEdge){
                            facesInComponent = currComponent;
                            shouldBreak = true;
                            break;
                        }
                    }
                }
                if (shouldBreak) break;
            }
            double maxCurl = -DBL_MAX;
            int maxEdge;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                for (Edge e : currFace.adjacentEdges()){
                    //check if current isoline crosses this edge
                    if ((minCurlIsoVal > globalTimeFunction[e.halfedge().tailVertex()] && minCurlIsoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    || (minCurlIsoVal > globalTimeFunction[e.halfedge().tipVertex()] && minCurlIsoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] > maxCurl){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_pair(maxEdge, minCurlEdge.getIndex()));
        }
        ctr++;
        if (ctr == numSingularityPairs) break;
    }
    return singEdgePairs;
}

std::vector<std::tuple<std::pair<int, int>, double>> findEdgeSingularityPairsUsingTimeFunctionIsoVals(VertexPositionGeometry& globalGeometry, 
                                                            EdgeLengthGeometry& gluedGeometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
                                                            EdgeData<double>& curl, VertexData<double>& globalTimeFunction, 
                                                            double stepSize, std::map<int, int>& hashedUsedIsoVals, 
                                                            FaceData<Vector3>& timeFunctionGrad, int numSingularityPairs){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    //first map avg curl values to iso values 
    //sort the keys in descending order
    std::map<double, double, std::greater<double>> curlToIsoVal;

    double curr = stepSize;
    double eps = 1e-5;
    double end = 1.0 - eps;
    bool skipFlag = false;
    double currAvgDeviation = 0.0;
    double currDeviationSum = 0.0;
    std::vector<std::tuple<std::pair<int, int>, double>> singEdgePairs;
    //should figure out some constant step size
    stepSize = 0.05;
    //changing the alignment strongly affects the helicing condition due to the path constraints
    double alignment = 0.99;

    while (curr < end){

        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, curr);

        //reset values
        currDeviationSum = 0.0;
        currAvgDeviation = 0.0;
        int numEdges = 0;
        for (int i = 0; i < f.size(); i++){
            Face currFace = globalMesh.face(f[i]);
            Vector3 grad = timeFunctionGrad[currFace];
            for (Edge e : currFace.adjacentEdges()){
                Vector3 heVector = (globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()]).normalize();
                Vector3 heTwinVector = (globalGeometry.vertexPositions[e.halfedge().tailVertex()] - globalGeometry.vertexPositions[e.halfedge().tipVertex()]).normalize();
                if (dot(heVector, grad) < alignment && dot(heTwinVector, grad) < alignment) continue; //skip edges that are poorly aligned with the time function gradient
                //check if current isoline crosses this edge
                if ((curr > globalTimeFunction[e.halfedge().tailVertex()] && curr < globalTimeFunction[e.halfedge().tipVertex()])
                || (curr > globalTimeFunction[e.halfedge().tipVertex()] && curr < globalTimeFunction[e.halfedge().tailVertex()])){
                    currDeviationSum += std::fabs(curl[e]);
                    numEdges++;
                }
            }
        }
        //currAvgDeviation = currDeviationSum / numEdges;
        currAvgDeviation = currDeviationSum;
        curlToIsoVal[currAvgDeviation] = curr;
        curr += stepSize;
    }

    int ctr = 0;
    for (auto entry : curlToIsoVal){
        double isoVal = entry.second;
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getIsoLine(V, F, globalTimeFunction, isoVal);
        // auto isoline = polyscope::registerCurveNetwork("isoline at " + std::to_string(isoVal), iV, iE);
        // std::cout << "average curl at isoval " << isoVal << " is " << entry.first << std::endl;
        // isoline->setEnabled(false);
        // isoline->setRadius(0.001);

        std::vector<std::vector<int>> facesPerComponent = findConnectedComponents(gluedGeometry, f);
        if (facesPerComponent.size() == 1){//there's only one connected component in this isoline
            std::vector<int> facesInComponent = facesPerComponent[0];
            double maxCurl = -DBL_MAX;
            double minCurl = DBL_MAX;
            int maxEdge, minEdge;
            for (int i = 0; i < facesInComponent.size(); i++){
                Face currFace = globalMesh.face(facesInComponent[i]);
                Vector3 grad = timeFunctionGrad[currFace];
                for (Edge e : currFace.adjacentEdges()){
                    Vector3 heVector = (globalGeometry.vertexPositions[e.halfedge().tipVertex()] - globalGeometry.vertexPositions[e.halfedge().tailVertex()]).normalize();
                    Vector3 heTwinVector = (globalGeometry.vertexPositions[e.halfedge().tailVertex()] - globalGeometry.vertexPositions[e.halfedge().tipVertex()]).normalize();
                    if (dot(heVector, grad) < alignment && dot(heTwinVector, grad) < alignment) continue; //skip edges that are poorly aligned with the time function gradient
                        //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                        || (isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] > maxCurl){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                        else if (curl[e] < minCurl){
                            minCurl = curl[e];
                            minEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_tuple(std::make_pair(maxEdge, minEdge), isoVal));
        }
        else{
            //first need to find component with max curl
            std::vector<int> componentWithMaxAvgCurl;
            double maxAvgCurl = -DBL_MAX;
            for (int i = 0; i < facesPerComponent.size(); i++){
                std::vector<int> currComponent = facesPerComponent[i];
                double currAvg = 0.;
                for (int j = 0; j < currComponent.size(); j++){
                    Face currFace = globalMesh.face(currComponent[j]);
                    for (Edge e : currFace.adjacentEdges()){
                        //check if current isoline crosses this edge
                        if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                            || (isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                            currAvg += std::fabs(curl[e]);
                        }
                    }
                }
                //currAvg /= currComponent.size();
                if (currAvg > maxAvgCurl){
                    maxAvgCurl = currAvg;
                    componentWithMaxAvgCurl = currComponent;
                }
            }
            double maxCurl = -DBL_MAX;
            double minCurl = DBL_MAX;
            int maxEdge, minEdge;
            for (int i = 0; i < componentWithMaxAvgCurl.size(); i++){
                Face currFace = globalMesh.face(componentWithMaxAvgCurl[i]);
                for (Edge e : currFace.adjacentEdges()){
                    //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                    || (isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] > maxCurl){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                        else if (curl[e] < minCurl){
                            minCurl = curl[e];
                            minEdge = e.getIndex();
                        }
                    }
                }
            }
            singEdgePairs.push_back(std::make_tuple(std::make_pair(maxEdge, minEdge), isoVal));
        }
        ctr++;
        if (ctr == numSingularityPairs) break;
    }
    return singEdgePairs;

    

}

//find edge singularity pairs 
//this method samples isolines of the time function 
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
                        if (curl[e] > maxCurl && (!usedEdges.count(e.getIndex()) > 0)){
                            maxCurl = curl[e];
                            maxEdge = e.getIndex();
                        }
                    }
                    //check if current isoline crosses this edge
                    if ((isoVal > globalTimeFunction[e.halfedge().tailVertex()] && isoVal < globalTimeFunction[e.halfedge().tipVertex()])
                        ||(isoVal > globalTimeFunction[e.halfedge().tipVertex()] && isoVal < globalTimeFunction[e.halfedge().tailVertex()])){
                        if (curl[e] < minCurl && (!usedEdges.count(e.getIndex()) > 0)){
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

void revealCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, std::map<int, int>& vertexMap,  Eigen::SparseMatrix<double, Eigen::RowMajor>& G) {
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    std::vector<int> bdyEdges = gbModel.getBdyEdges();

    try {
        GRBEnv env(true);
        env.start();

        GRBModel model = GRBModel(env);
        
        // 1-form defined over halfedges
        std::vector<GRBVar> sigma;
        for (size_t i = 0; i < gluedMesh.nHalfedges(); i++) {
            GRBVar sigma_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
            sigma.push_back(sigma_i);
        }

        // Constraint: sigma at boundary edges should be 0
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

        // Constraint: ||grad(sigma)||^2  = 1 on each face
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

            GRBQuadExpr gradUnorm2 = 0;
            for (int i = 0; i < 3; i++)
                gradUnorm2 += gradU[f.getIndex()][i] * gradU[f.getIndex()][i];
            model.addQConstr(gradUnorm2 == 1);
        }

        // Objective: minimize curl on edges
        GRBQuadExpr obj = 0;
        for (size_t i = 0; i < gluedMesh.nEdges(); i++) {
            Halfedge he = gluedMesh.edge(i).halfedge();
            GRBLinExpr curl = sigma[he.getIndex()] + sigma[he.twin().getIndex()];
            obj += curl*curl;
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();

    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }


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
//check if an edge pair is valid or not given a list of edge singularity pairs 
bool isValidEdgePair(VertexPositionGeometry& globalGeometry, EdgeData<double>& edgeSingularities, std::pair<int, int>& p,
                    std::map<int, int>& hashedUsedIsoVals, double isoVal){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    if (edgeSingularities[p.first] != 0 || edgeSingularities[p.second] != 0){//don't select edges we've seen before 
        return false;
    }

    if (hashedUsedIsoVals.count(hashFloatQuantized(isoVal)) > 0){//don't reuse isovals
        return false;
    }


    //not really sure about whether we want to allow boundary edges to be singular or not
    // if (globalMesh.edge(p.first).isBoundary() || globalMesh.edge(p.second).isBoundary()){//don't select boundary edges
    //     return false;
    // }
    return true;
}




//@clean 
//the below two functions use the energy min ||\del sigma_i - \nabla \sigma_{i - 1} / ||\nabla \sigma_{i - 1}|| ||^2
//if i = 1, use \nabla h / ||\nabla h||
std::tuple<CornerData<double>, EdgeData<double>> implCourseHarmonic1Form(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                                                    VertexData<double>& globalTimeFunction, FaceData<Vector3>& globalTimeFunctionGradientsNormalized,
                                                                    std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, polyscope::SurfaceMesh& psMesh,
                                                                    globalBoundaryConditions& boundaryConditions, double period,
                                                                    Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double, Eigen::RowMajor>& G,
                                                                    FaceData<Vector3>& courseOneFormGrad, std::map<int, std::vector<Halfedge>>& gluedOneRingMap,
                                                                    std::vector<std::vector<double>> allSaddleLoops, std::vector<std::vector<double>> homologyGenerators,
                                                                    Options &opts){
    
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    //gluedGeometry.requireCotanLaplacian();
    //gluedGeometry.requireVertexGalerkinMassMatrix();
    //Eigen::SparseMatrix<double> L = gluedGeometry.cotanLaplacian;
    //Eigen::SparseMatrix<double> M = gluedGeometry.vertexGalerkinMassMatrix;
    double sum = 0;
    // Create the Heat Method solver
    HeatMethodDistanceSolver heatSolver(gluedGeometry);
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
    //singular edges in the gurobi optimization (in a better representation)
    std::vector<int> edgeIndices(gluedMesh.nEdges(), 0);
    //heat source vertices 
    std::vector<Vertex> heatSourceVerts;
    //sing edge pairs (max/min curl edges)
    std::vector<std::pair<int, int>> singEdgePairs;
    //accepted sing edge pairs 
    std::vector<std::pair<int, int>> acceptedSingEdgePairs;
    //make a map of singular edges we've used so far 
    //so we don't re-use edges
    std::map<int, int> usedEdges;
    //path constraints
    std::vector<double> globalPath;
    std::vector<double> gluedPath;
    //number of runs of the optimization
    int numRuns = 0;
    //max number of singularity pairs to check for insertion
    int topPairs = 10;
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
    //epsilon values 
    double eps = 1e-8;
    //distance from unit norm 
    double oldDistance = 0.;
    double newDistance = 0.;
    //striping information 
    //global data
    CornerData<double> stripeValuesSigmaCourse(gluedMesh);
    FaceData<int> stripeIndicesSigmaCourse(globalMesh);
    //stripe curve network information
    //unique vertices and unique edges
    std::vector<Vector3> uniquePos;
    std::vector<std::array<int, 2>> uniqueEdges;
    //conencted components identified by their component id 
    std::vector<StripeConnectedComponent> components;
    //isoval we'll be tracing
    double isoVal;
    //break flag - to stop searching for singularities
    bool toBreak = false;

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
    model.setEdgeIndices(edgeIndices);

    //find the step size to sample level sets at
    double avgSum = 0.;
    double stepSize = 0.01;
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

    //masked vertex for bent cylinder masking experiment/figure
    // std::vector<int> maskedVerticesIds = {2624, 2620, 624, 2824, 2820, 2816, 2812, 528, 3016, 3012, 3008, 525,
    //                               526, 527, 52, 620, 621, 622, 623, 64, 716, 717, 2643, 2639, 619, 2843, 2839, 2835, 
    //                               2831, 523, 3035, 3031, 3027, 3026, 3030, 3034, 522, 2830, 2834, 2838, 2842, 618, 2638, 2642, 
    //                               2641, 2637, 617, 2841, 2837, 2833, 2820, 521, 3033, 3029, 3024, 3028, 3032, 520, 2828, 2832, 2836, 2840, 616,
    //                               2636, 2640, 2829, 3025};
    // VertexData<double> maskedVertices(globalMesh, 0.0);
    // for (int v : maskedVerticesIds){
    //     maskedVertices[v] = 1.0;
    //     //add the masked vertices to the heat sources
    //     heatSourceVerts.push_back(gluedMesh.vertex(vertexMap[v]));
    // }
    // psMesh.addVertexScalarQuantity("masked vertices", maskedVertices)


    //handle saddle loops in the intrinsic setting
    for(int i = 0; i < allSaddleLoops.size(); i++){
        edgePathConstraints.push_back(std::make_pair(allSaddleLoops[i], 0.0));
    }
    for (int i = 0; i < allSaddleLoops.size(); i++){
        std::vector<double> path = allSaddleLoops[i];
        for (int j = 0; j < path.size(); j++){
            if (std::fabs(path[j]) > 0){
                Edge e = gluedGeometry.mesh.edge(j);
                heatSourceVerts.push_back(e.halfedge().tailVertex());
                heatSourceVerts.push_back(e.halfedge().tipVertex());
            }
        }
    }
    

    //solve the model without any singularities 
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    std::cout << "objective without placing any singularities = " << currObj << std::endl;
    numRuns++;
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);

    // Plot gluedSigmaTilde and its stripe values
    //psMesh.addHalfedgeScalarQuantity("gluedSigmaTilde after " + std::to_string(numRuns-1) + " runs", gluedSigmaTilde);
    //psMesh.addCornerScalarQuantity("stripeValuesSigmaCourse after " + std::to_string(numRuns-1) + " runs", prepareCornerData(stripeValuesSigmaCourse));

    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, edgeMap);
    auto courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numRuns - 1) + 
                                                    " singularity insertions", uniquePos, uniqueEdges);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);
    gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    courseOneFormGrad = gradSigmaTilde;


    //update the gradients for the next iteration of the model 
    //model.setFaceGradients(gradSigmaTilde);
    //newDistance = computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde);
    oldObj = currObj;
    //std::cout << "distance from unit norm after " << std::to_string(numRuns - 1) << " singularity insertion " << newDistance << std::endl;
    //oldDistance = newDistance;


    //---------------------Testing-------------------------//

    //ONLY IMPLEMENTED IN THE GLOBAL SETTING FOR NOW

    //rotate the gradients so that we measure curl in the wale direction
    // globalGeometry.requireFaceNormals();
    // //rotate the gradients to measure the wale curl 
    // for (Face f : globalMesh.faces()){
    //     globalTimeFunctionGradientsNormalized[f] = globalTimeFunctionGradientsNormalized[f].rotateAround(globalGeometry.faceNormals[f], M_PI/2.);
    // }

    //Attempting to find the center of distributions using Vector Heat Method
    VertexData<double> curlMeasure = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                    globalTimeFunctionGradientsNormalized, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);

    psMesh.addVertexScalarQuantity("initial curl measure ", curlMeasure);



    // Cap curl measure to avoid high concentration of singularities
    for (Vertex v : globalMesh.vertices()) {
        curlMeasure[v] = fmin(curlMeasure[v], period / (3*globalGeometry.vertexDualAreas[v]));
        curlMeasure[v] = fmax(curlMeasure[v], -period / (3*globalGeometry.vertexDualAreas[v]));
    }
                                    
    //flag of faces we've placed singularities on 
    std::vector<int> usedFaces(globalMesh.nFaces(), 0);
    
    psMesh.addVertexScalarQuantity("initial curl measure (capped)", curlMeasure);
    //now divide the measure into positive and negative 
    VertexData<double> posMeasure(globalMesh, 0.0);
    VertexData<double> negMeasure(globalMesh, 0.0);
    double totalPosMeasure = 0;
    double totalNegMeasure = 0;
    for (Vertex v : globalMesh.vertices()){
        if (curlMeasure[v] > 0){
            posMeasure[v] = curlMeasure[v];
            totalPosMeasure += globalGeometry.vertexDualAreas[v] * posMeasure[v];
        }
        else{
            negMeasure[v] = std::fabs(curlMeasure[v]);
            totalNegMeasure += globalGeometry.vertexDualAreas[v] * negMeasure[v];
        }
    }
    double avgTotalMeasure = (totalPosMeasure + totalNegMeasure) / 2;

    //debugging on a simple square
    // for (Vertex v : globalMesh.vertices()){
    //     //if (globalGeometry.vertexPositions[v].z < 0 && globalGeometry.vertexPositions[v].x < 0) posMeasure[v] = 1.0;
    //     // posMeasure[v] = std::fabs(globalGeometry.vertexPositions[v].x);
    //     posMeasure[v] = globalTimeFunction[v];
    //     // posMeasure[v] = globalGeometry.vertexPositions[v].x * globalGeometry.vertexPositions[v].x;
    //     // posMeasure[v] = std::exp(globalGeometry.vertexPositions[v].x);
    // }

    psMesh.addVertexScalarQuantity("initial positive measure ", posMeasure);
    psMesh.addVertexScalarQuantity("initial negative measure ", negMeasure);
    std::unique_ptr<ManifoldSurfaceMesh> manifoldGlobalMesh = globalMesh.toManifoldMesh();
    int numSings = 0.;
    VoronoiOptions posOptions = defaultVoronoiOptions;
    posOptions.nSites = std::round(avgTotalMeasure / period);
    posOptions.useDelaunay = false;
    posOptions.computeDistributions = true;
    posOptions.iterations = 500;
    VoronoiResult posVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(*manifoldGlobalMesh, globalGeometry, posOptions, posMeasure, psMesh);
    std::vector<Vector3> positiveCenters;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
       positiveCenters.push_back(posVoronoiCenters.siteLocations[i].interpolate(globalGeometry.vertexPositions));
    }
    polyscope::registerPointCloud("positive voronoi sites (course)", positiveCenters);
    std::vector<VertexData<double>> posSiteDistributions = posVoronoiCenters.siteDistributions;

    //print the masses
    for (size_t i = 0; i < posSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){
            if (std::fabs(posSiteDistributions[i][v]) > 1e-8){
                mass += posSiteDistributions[i][v];
            }
        }
        std::cout << "positive mass at site " << i << " = " << mass << std::endl;
        psMesh.addVertexScalarQuantity("positve site distribution (course)" + std::to_string(i), posSiteDistributions[i]);
    }

   
    VoronoiOptions negOptions = defaultVoronoiOptions;
    negOptions.nSites = posOptions.nSites;
    negOptions.useDelaunay = false;
    negOptions.computeDistributions = true;
    negOptions.iterations = 500;
    std::cout << "# positive sites " << posOptions.nSites << std::endl;
    std::cout << "# negative sites " << negOptions.nSites << std::endl;
    VoronoiResult negVoronoiCenters = computeGeodesicCentroidalVoronoiTessellationWithWeights(*manifoldGlobalMesh, globalGeometry, negOptions, negMeasure, psMesh);
    std::vector<Vector3> negativeCenters;
    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
       negativeCenters.push_back(negVoronoiCenters.siteLocations[i].interpolate(globalGeometry.vertexPositions));
    }
    polyscope::registerPointCloud("negative voronoi sites (course)", negativeCenters);
    std::vector<VertexData<double>> negSiteDistributions = negVoronoiCenters.siteDistributions;

    //print the masses
    for (size_t i = 0; i < negSiteDistributions.size(); i++){
        double mass = 0;
        for (Vertex v : globalMesh.vertices()){
            if (std::fabs(negSiteDistributions[i][v]) > 1e-8){
                mass += negSiteDistributions[i][v];
            }
        }
        std::cout << "negative mass at site " << i << " = " << mass << std::endl;
        psMesh.addVertexScalarQuantity("negative site distribution (course)" + std::to_string(i), negSiteDistributions[i]);
    }

    //store a vector of vertex ids and singularity indices (for optimal matching)
    std::vector<std::pair<Vertex, int>> singularities;
    for (int i = 0; i < posVoronoiCenters.siteLocations.size(); i++){
        SurfacePoint facePoint = posVoronoiCenters.siteLocations[i].inSomeFace();
        Face f = facePoint.face;
        Edge singularEdge;
        double maxDotProd = -DBL_MAX;
        // if (facePoint.faceCoords.x == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().vertex(), 1));
        // else if (facePoint.faceCoords.y == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().vertex(), 1));
        // else if (facePoint.faceCoords.z == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().next().vertex() , 1));
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - 
                                globalGeometry.vertexPositions[he.tailVertex()]).normalize();
            if (std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[f])) > maxDotProd){
                maxDotProd = std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[he.face()]));
                singularEdge = he.edge();
            }
        }
        //increment the indices really close so indices cancel out
        edgeIndices[singularEdge.getIndex()] = -1.0;
        edgeSingularities[singularEdge] = 1.0;

    }

    // for (int i = 0; i < posVoronoiCenters.steps.size(); i++){
    //     std::vector<SurfacePoint> steps = posVoronoiCenters.steps[i]; //steps for this site
    //     std::vector<Vector3> stepPos;
    //     for (int i = 0; i < steps.size(); i++){
    //         stepPos.push_back(steps[i].interpolate(globalGeometry.vertexPositions));
    //     }
    //     polyscope::registerPointCloud("pos site " + std::to_string(i) + " steps ", stepPos)->setEnabled(false);
    // }

    for (int i = 0; i < negVoronoiCenters.siteLocations.size(); i++){
        SurfacePoint facePoint = negVoronoiCenters.siteLocations[i].inSomeFace();
        Face f = facePoint.face;
        Edge singularEdge;
        double maxDotProd = -DBL_MAX;
        // if (facePoint.faceCoords.x == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().vertex(), -1));
        // else if (facePoint.faceCoords.y == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().vertex(), -1));
        // else if (facePoint.faceCoords.z == std::max({facePoint.faceCoords.x, facePoint.faceCoords.y, facePoint.faceCoords.z})) singularities.push_back(std::make_pair(f.halfedge().next().next().vertex() , -1));
        for (Halfedge he : f.adjacentHalfedges()){
            Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - 
                                globalGeometry.vertexPositions[he.tailVertex()]).normalize();
            if (std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[f])) > maxDotProd){
                maxDotProd = std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[he.face()]));
                singularEdge = he.edge();
            }
        }
        //increment the indices really close so indices cancel out
        edgeIndices[singularEdge.getIndex()] = 1.0;
        edgeSingularities[singularEdge] = -1.0;
    }

    
    //perform optimal matching 
    // HalfedgeData<double> heWeights = gluedHeWeights;
    // std::vector<std::pair<Vertex, Vertex>> matchedVertices = performOptimalMatching(globalGeometry, heWeights, singularities);  

    // for (int i = 0; i < matchedVertices.size(); i++){
    //     std::pair<Vertex, Vertex> p = matchedVertices[i];
    //     std::vector<Vector3> pC; 
    //     pC.push_back(globalGeometry.vertexPositions[p.first]);
    //     pC.push_back(globalGeometry.vertexPositions[p.second]);
    //     polyscope::registerPointCloud("matched pair " + std::to_string(i), pC)->setEnabled(false);
    // }

    //set the edge indices from the singularities found using the OT method 
    //select the outgoing halfedge that is most aligned with the gradient of the time function
    // for (auto p : singularities){
    //     // if (usedFaces[p.first.halfedge().face().getIndex()] == 0 &&
    //     //     usedFaces[p.first.halfedge().twin().face().getIndex()] == 0){
            
    //     //     edgeIndices[p.first.halfedge().edge().getIndex()] = -1.0 * p.second;
    //     //     usedFaces[p.first.halfedge().face().getIndex()] = 1;
    //     //     usedFaces[p.first.halfedge().twin().face().getIndex()] = 1;
    //     // }

    //     //find the halfedge that is most aligned with the gradient of the time function
    //     Edge singularEdge;
    //     double maxDotProd = -DBL_MAX;
    //     for (Halfedge he : p.first.outgoingHalfedges()){
    //         Vector3 heVec = (globalGeometry.vertexPositions[he.tipVertex()] - 
    //                             globalGeometry.vertexPositions[he.tailVertex()]).normalize();
    //         if (std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[he.face()])) > maxDotProd){
    //             maxDotProd = std::fabs(dot(heVec, globalTimeFunctionGradientsNormalized[he.face()]));
    //             singularEdge = he.edge();
    //         }
    //     }
    //     //increment the indices really close so indices cancel out
    //     edgeIndices[singularEdge.getIndex()] += (-1.0 * p.second);
    //     edgeSingularities[singularEdge] = p.second;
    // }


    model.setEdgeIndices(edgeIndices);
    model.setHomologyGenerators(homologyGenerators);
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    //plot the stripes
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, edgeMap);
    courseStripes = polyscope::registerCurveNetwork("our course stripes ", uniquePos, uniqueEdges);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);

    //-------------------End of testing-------------------//

    /** 
    //compute curl quantities without impulse function
    vertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                    gradSigmaTilde, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);
    psMesh.addVertexScalarQuantity("vertex curl after " + std::to_string(numRuns - 1) + " singularity insertions (before subtracting)", vertexCurl);
    
    
    if (heatSourceVerts.size() != 0){//keep singularities away from saddle vertices
        allDist = heatSolver.computeDistance(heatSourceVerts);
        //for all the source vertices, find the max value 
        for (Vertex v : heatSourceVerts){
            maxSourceVal = std::max(maxSourceVal, allDist[v]);
        }
        //find the maximum value over all the distances
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
        allDistGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, allDist, vertexMap);
        psMesh.addVertexScalarQuantity("initial all distance course", allDistGlobal);
        for (Vertex v : globalMesh.vertices()){
            //vertexCurl[v] = (1 - exp(-pow(allDistGlobal[v], 2.) / (2 * pow(period, 2)))) * vertexCurl[v];
            courseWeighting[v] = (allDistGlobal[v] > 2*period);
            vertexCurl[v] = courseWeighting[v] * vertexCurl[v];
        }
    }
    edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (before subtracting)", edgeCurl);
    //compute virtual sigma
    std::tie(virtualSigmaTilde, virtualSigmaObj) = computeCourseVirtualSigma(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    //subtract off the impulse function
    gluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
    adjustedGradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
    //compute curl quantities after impulse function
    vertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                    adjustedGradSigmaTilde, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);
    if (heatSourceVerts.size() != 0){
        maxVal = std::numeric_limits<double>::min();
        maxSourceVal = std::numeric_limits<double>::min();
        //for all the source vertices, find the max value 
        for (Vertex v : heatSourceVerts){
            maxSourceVal = std::max(maxSourceVal, allDist[v]);
        }
        //find the maximum value over all the distances
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
        edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
        psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertion (after subtracting w/o mask)", edgeCurl);
        allDistGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, allDist, vertexMap);
        //adjust with the Gaussian mask
        //have a hard mask in the course direction
        for (Vertex v : globalMesh.vertices()){
            //vertexCurl[v] = (1 - exp(-pow(allDistGlobal[v], 2.) / (2 * pow(period, 2)))) * vertexCurl[v];
            courseWeighting[v] = (allDistGlobal[v] > 2*period);
            //courseWeighting[v] = (allDistGlobal[v] > period);//reduce the radius of the Gaussian (for the figure)
            vertexCurl[v] = courseWeighting[v] * vertexCurl[v];
        }
    }
    edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
    psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertions (after subtracting)", edgeCurl);
    int numSingularities = 0;
    while(true){
        //finding edge singularity pairs by sampling time function isolines
        //return type is a a vector of tuples
        std::vector<std::tuple<std::pair<int, int>, double>> singEdgePairs = findEdgeSingularityPairsUsingTimeFunctionIsoVals(globalGeometry,  gluedGeometry, V, F, 
                                                                            edgeCurl, globalTimeFunction, stepSize, hashedUsedIsoVals, globalTimeFunctionGradientsNormalized, 
                                                                            topPairs);

        //finding edge singularity pairs by just looking for max curl edges
        // singEdgePairs = findEdgeSingularityPairsUsingMaxCurls(globalGeometry, gluedGeometry, V, F, edgeCurl,
        //                 globalTimeFunction, hashedUsedIsoVals, topPairs);

        
        //finding edge singularity pair using stripe isolines
        // doesn't return a tuple
        // singEdgePairs = findEdgeSingularityPairFromStripeIsoVals(globalGeometry, gluedGeometry, 
        //                                                         edgeCurl, edgeMap, components, topPairs);
        int numSkips = 0;
        int maxSkips = topPairs;
        
        for (auto s : singEdgePairs){
            
            //if the return type is not a tuple
            //std::pair<int, int> singEdgePair = s;

            //if the return type is a tuple 
            std::pair<int, int> singEdgePair = std::get<0>(s);
            isoVal = std::get<1>(s);

            if (!isValidEdgePair(globalGeometry, edgeSingularities, singEdgePair,
                                hashedUsedIsoVals, isoVal)){
                std::cout << "pair = " << singEdgePair.first << ", " << singEdgePair.second << std::endl;
                std::cout << "Invalid edge pair. Trying next pair...." << std::endl;
                //try the top pairs
                numSkips++;
                std::cout << "numSkips = " << numSkips << std::endl;
                if (numSkips >= maxSkips){//we've exhausted all possible pairs and there is no valid pair left
                    toBreak = true;
                    break;
                }
                else{
                    continue;//skip edge singularity pairs that are invalid
                }
            }
            //test model 
            //to test a new pair of singularities
            Model testModel = model; 
            std::vector<std::pair<std::vector<double>, double>> testEdgePathConstraints = edgePathConstraints;
            std::vector<int> testEdgeIndices = edgeIndices;
            std::vector<std::pair<int, int>> testAcceptedSingEdgePairs = acceptedSingEdgePairs;
            //try a new pair of singularities
            testEdgeIndices[edgeMap[singEdgePair.first]] = -1;//positive edge (sign gets flipped)
            testEdgeIndices[edgeMap[singEdgePair.second]] = 1;//negative edge (sign gets flipped)
            testAcceptedSingEdgePairs.push_back(singEdgePair);
            bool connectSaddles = false;
            std::tie(globalPath, gluedPath) = constructEdgePath(globalGeometry, gluedGeometry, globalMesh.edge(singEdgePair.first), globalMesh.edge(singEdgePair.second),
                                                                vertexMap, edgeMap, globalTimeFunctionGradientsNormalized, gluedHeWeights, gluedSigmaTilde, connectSaddles);
            testModel.setSingularEdges(testAcceptedSingEdgePairs);
            testEdgePathConstraints.push_back(std::make_pair(gluedPath, 0.));
            testModel.setEdgePathConstraints(testEdgePathConstraints);
            testModel.setEdgeIndices(testEdgeIndices);
            //solve the model with singularities
            std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, testModel, vertexMap, G, psMesh);
            numRuns++;
            gradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, gluedSigmaTilde);
            //newDistance = computeDistanceFromUnitNorm(globalGeometry, gradSigmaTilde);
            if (currObj > oldObj){
                std::cout << "currObj = " << currObj << std::endl;
                std::cout << "oldObj = " << oldObj << std::endl;
                std::cout << "objective is not improving. Trying the next top pair... " << std::endl;
                numSkips++;
                std::cout << "numSkips = " << numSkips << std::endl;
                if (numSkips >= maxSkips){//we've exhausted all possible pairs and nothing is improving the objective
                    std::cout << "breaking after " << std::to_string(numRuns - 2) << " singularity attempts " << std::endl;
                    std::cout << "currObj = " << currObj << std::endl;
                    std::cout << "oldObj = " << oldObj << std::endl;
                    toBreak = true;
                    break;
                }
                continue;
            }
            else{
                hashedUsedIsoVals[hashFloatQuantized(isoVal)] = 1;
                numSingularities++;
                //gotten to a valid pair of edge indices that improves the objective
                //append to accepted singEdgePairs
                acceptedSingEdgePairs.push_back(singEdgePair);
                //pair of singular edges
                edgeIndices[edgeMap[singEdgePair.first]] = -1;
                edgeIndices[edgeMap[singEdgePair.second]] = 1;
                heatSourceVerts.push_back(gluedMesh.edge(edgeMap[singEdgePair.first]).halfedge().tailVertex());
                heatSourceVerts.push_back(gluedMesh.edge(edgeMap[singEdgePair.first]).halfedge().tipVertex());
                heatSourceVerts.push_back(gluedMesh.edge(edgeMap[singEdgePair.second]).halfedge().tailVertex());
                heatSourceVerts.push_back(gluedMesh.edge(edgeMap[singEdgePair.second]).halfedge().tipVertex());
                edgeSingularities[globalMesh.edge(singEdgePair.first)] = 1.0;
                edgeSingularities[globalMesh.edge(singEdgePair.second)] = -1.0;
                //make a curve network to viz singular edges
                std::vector<Vector3> singPos; 
                std::vector<std::array<int, 2>> edges;
                singPos.push_back(globalGeometry.vertexPositions[globalMesh.edge(singEdgePair.first).halfedge().tailVertex()]);
                singPos.push_back(globalGeometry.vertexPositions[globalMesh.edge(singEdgePair.first).halfedge().tipVertex()]);
                edges.push_back(std::array<int, 2>{0, 1});
                //don't retake edges from another path
                updateGluedHalfedgeWeights(globalGeometry, gluedGeometry, gluedPath, gluedHeWeights);
                if (opts.showAllIterations) psMesh.addEdgeScalarQuantity("path for " + std::to_string(numSingularities) + " singularity", globalPath);
                edgePathConstraints.push_back(std::make_pair(gluedPath, 0.));
                model.setSingularEdges(acceptedSingEdgePairs);
                model.setEdgePathConstraints(edgePathConstraints);
                model.setEdgeIndices(edgeIndices);
                std::cout << "not breaking after " << std::to_string(numSingularities) << " singularity insertions " << std::endl;
                std::cout << "newDistance = " << newDistance << std::endl;
                std::cout << "oldDistance = " << oldDistance << std::endl;
                //oldDistance = newDistance;
                oldObj = currObj;
                courseOneFormGrad = gradSigmaTilde;
                std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
                std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, 
                                                edgeMap);
                //std::cout << "Number of components after " << std::to_string(numSingularities) << " singularity insertions is " << components.size() << std::endl;
                if (opts.showAllIterations) {
                    psMesh.addCornerScalarQuantity("stripeValuesSigmaCourse after " + std::to_string(numRuns-1) + " runs", prepareCornerData(stripeValuesSigmaCourse));
                    psMesh.addHalfedgeScalarQuantity("sigma values after " + std::to_string(numRuns-1) + " runs", gluedSigmaTilde);
                    courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numSingularities) + 
                                                    " singularity insertions", uniquePos, uniqueEdges);
                    courseStripes -> setRadius(0.001);
                    courseStripes -> setEnabled(false);
                }
                
                //update the gradients for the next round of the model 
                //model.setFaceGradients(gradSigmaTilde);
                //compute curl quantities without accounting for impulse function
                vertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                    gradSigmaTilde, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);
                edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
                if (opts.showAllIterations) psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numSingularities) + " singularity insertions (before subtracting)", edgeCurl);
                //compute virtual sigma
                std::tie(virtualSigmaTilde, virtualSigmaObj) = computeCourseVirtualSigma(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
                //psMesh.addHalfedgeScalarQuantity("virtualSigmaTilde after " + std::to_string(numRuns-1) + " runs", virtualSigmaTilde);
                std::cout << "virtual sigma objective = " << virtualSigmaObj << std::endl;
                //subtract off the impulse function
                HalfedgeData<double> adjustedGluedSigmaTilde = gluedSigmaTilde - virtualSigmaTilde;
                adjustedGradSigmaTilde = computeOneFormFaceGrad(globalGeometry, gluedGeometry, adjustedGluedSigmaTilde);
                //compute curl quantities after accounting for impulse function
                vertexCurl = computeCourseVertexCurl(globalGeometry, gluedGeometry, 
                                                adjustedGradSigmaTilde, gluedOneRingMap, edgeIndices, heatSolver, vertexMap);
                
                allDist = heatSolver.computeDistance(heatSourceVerts);
                maxVal = std::numeric_limits<double>::min();
                maxSourceVal = std::numeric_limits<double>::min();
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
                edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
                if (opts.showAllIterations) psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numRuns - 1) + " singularity insertion (after subtracting w/o mask)", edgeCurl);
                allDistGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, allDist, vertexMap);
                //have a hard mask in the course direction
                for (Vertex v : globalMesh.vertices()){
                    //vertexCurl[v] = (1 - exp(-pow(allDistGlobal[v], 2.) / (2 * pow(period, 2)))) * vertexCurl[v];
                    courseWeighting[v] = (allDistGlobal[v] > 2*period);
                    //courseWeighting[v] = (allDistGlobal[v] > period);//reduce the radius of the Gaussian (for fig.)
                    vertexCurl[v] = courseWeighting[v] * vertexCurl[v];
                }
                psMesh.addVertexScalarQuantity("all distance", allDistGlobal);
                edgeCurl = computeVertexAveragedEdgeCurl(globalGeometry, vertexCurl);
                if (opts.showAllIterations) psMesh.addEdgeScalarQuantity("edge curl after " + std::to_string(numSingularities) + " singularity insertions (after subtracting)", edgeCurl);
                if (opts.showAllIterations) psMesh.addEdgeScalarQuantity("edge singularities after " + std::to_string(numSingularities) + " singularity insertion", edgeSingularities);
                //polyscope::registerCurveNetwork("edge singularities network after " + std::to_string(numSingularities) + " singularitiy insertion", singPos, edges);
                break;
            }
        }
        if (toBreak) break;//we are no longer improving the distance to unit norm

    }

    //solve the final model with the singularities and homology generators 
    std::cout << "solving final course stripes...." << std::endl;
    model.setHomologyGenerators(homologyGenerators);

    // Without helicing correction (for comparison)
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, edgeMap);
    courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numSingularities) + 
                                                    " singularity insertions (no correction)", uniquePos, uniqueEdges);

    // With helicing correction
    model.useHelicingCorrection = true;
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, edgeMap);
    //std::cout << "Number of components after " << std::to_string(numSingularities) << " singularity insertions is " << components.size() << std::endl;
    psMesh.addCornerScalarQuantity("stripeValuesSigmaCourse after " + std::to_string(numRuns-1) + " runs", prepareCornerData(stripeValuesSigmaCourse))
        ->setIsolinesEnabled(true)
        ->setIsolineWidth(period/2, false)
        ->setIsolineDarkness(0);
    
    psMesh.addHalfedgeScalarQuantity("1-form after " + std::to_string(numSingularities) + " singularity insertions", gluedSigmaTilde);
    
    courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes after " + std::to_string(numSingularities) + 
                                                    " singularity insertions", uniquePos, uniqueEdges);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);
    registerShortRows("short rows after " + std::to_string(numSingularities) + " singularity insertions", components);    

    model.useHelicingCorrection = true;
    std::tie(gluedSigmaTilde, currObj) = computeCourseOneForm(globalGeometry, gluedGeometry, model, vertexMap, G, psMesh);
    std::tie(stripeValuesSigmaCourse, stripeIndicesSigmaCourse) = computeStripeValuesFromOneForm(globalGeometry, gluedGeometry, gluedSigmaTilde, period);
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesSigmaCourse, stripeIndicesSigmaCourse, period, edgeMap);
    
    //visualize the saddle vertices
    std::vector<Vector3> saddleVertexPositions;
    for (Vertex v : globalMesh.vertices()){
        int ctr = 0;
        for (Halfedge he : v.outgoingHalfedges()){
            if (sgn(gluedSigmaTilde[he]) != sgn(gluedSigmaTilde[he.next().next().twin()])) ctr++;
        }
        int index = (2 - ctr) / 2;
        if (index == -1) saddleVertexPositions.push_back(globalGeometry.vertexPositions[v]);
    }
    polyscope::registerPointCloud("saddle vertices", saddleVertexPositions);

    // Plot stripe values with offset (to debug knit graph)
    CornerData<double> stripeValuesWithOffset = stripeValuesSigmaCourse;
    for (Corner co : globalMesh.corners())
        stripeValuesWithOffset[co] -= period/4;
    std::tie(uniquePos, uniqueEdges, components) = findStripeConnectedComponents(globalGeometry, gluedGeometry, stripeValuesWithOffset, stripeIndicesSigmaCourse, period, edgeMap);
    //std::cout << "Number of components after " << std::to_string(numSingularities) << " singularity insertions is " << components.size() << std::endl;
    courseStripes = polyscope::registerCurveNetwork("sigma tilde stripes with offset", uniquePos, uniqueEdges);
    courseStripes -> setRadius(0.001);
    courseStripes -> setEnabled(false);

    std::cout << "Number of singularities inserted = " << numSingularities << std::endl;
    */

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
    std::vector<int> edgeIndices = gbModel.getEdgeIndices();
    std::vector<std::array<double, 3>> comparisonGrad = gbModel.getFaceGradients();
    std::vector<std::vector<double>> homologyGenerators = gbModel.getHomologyGenerators();

    std::cout << "size of singular edges = " << singularEdges.size() << std::endl;
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
        //model.getEnv().set(GRB_DoubleParam_TimeLimit, 30);
        //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_IntParam_SolutionLimit, 2);
        model.getEnv().set(GRB_IntParam_NumericFocus, 3);

        //add integer variables for all the generators
        std::vector<GRBVar> generatorIntegers;

        //add integer variable for all the generators
        for (size_t i = 0; i < homologyGenerators.size(); i++){
            GRBVar gen_i = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_INTEGER);
            generatorIntegers.push_back(gen_i);
        }

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
        for (Halfedge he : gluedMesh.halfedges()){
            if (edgeIndices[he.edge().getIndex()] != 0){
                model.addConstr(sigma[he.getIndex()] + sigma[he.twin().getIndex()] == edgeIndices[he.edge().getIndex()] * period);
                // int globalCornerTail = he.face().halfedge().corner().getIndex();
                // int globalCornerTip = he.face().halfedge().twin().corner().getIndex();
                // Vector3 heVector = globalGeometry.vertexPositions[globalMesh.corner(globalCornerTail).halfedge().tipVertex()] - 
                //                            globalGeometry.vertexPositions[globalMesh.corner(globalCornerTail).halfedge().tailVertex()];
                // Vector3 heTwinVector = globalGeometry.vertexPositions[globalMesh.corner(globalCornerTip).halfedge().tipVertex()] - 
                //                            globalGeometry.vertexPositions[globalMesh.corner(globalCornerTip).halfedge().tailVertex()];
                //just test it in the global setting first
                // Vector3 heVector = globalGeometry.vertexPositions[globalMesh.halfedge(he.getIndex()).tipVertex()] - 
                //                             globalGeometry.vertexPositions[globalMesh.halfedge(he.getIndex()).tailVertex()];
                // Vector3 heTwinVector = globalGeometry.vertexPositions[globalMesh.halfedge(he.getIndex()).tailVertex()] - 
                //                             globalGeometry.vertexPositions[globalMesh.halfedge(he.getIndex()).tipVertex()];
                // Vector3 gradientVector = Vector3{comparisonGrad[he.face().getIndex()][0], comparisonGrad[he.face().getIndex()][1], 
                //                                 comparisonGrad[he.face().getIndex()][2]};
                // if (dot(gradientVector, heVector) > dot(gradientVector, heTwinVector)){
                //     model.addConstr(sigma[he.getIndex()] >= sigma[he.twin().getIndex()]); 
                // }
                // else{
                //     model.addConstr(sigma[he.twin().getIndex()] >= sigma[he.getIndex()]);
                // }
            }
            else{
                model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
            }
        }
        
        //constraint: 
        //equality constraint across singular edges 
        if (gbModel.forceSaddle) {
            for (std::pair<int, int> s : singularEdges){
                Halfedge he1 = gluedMesh.edge(s.first).halfedge();
                Halfedge he1Twin = he1.twin();
                Halfedge he2 = gluedMesh.edge(s.second).halfedge();
                Halfedge he2Twin = he2.twin();

                Vector3 he1Vector = globalGeometry.vertexPositions[he1.tipVertex()] - globalGeometry.vertexPositions[he1.tailVertex()];
                he1Vector = he1Vector.normalize();
                Vector3 he1TwinVector = globalGeometry.vertexPositions[he1Twin.tipVertex()] - globalGeometry.vertexPositions[he1Twin.tailVertex()];
                he1TwinVector = he1TwinVector.normalize();
                Vector3 he2Vector = globalGeometry.vertexPositions[he2.tipVertex()] - globalGeometry.vertexPositions[he2.tailVertex()];
                he2Vector = he2Vector.normalize();
                Vector3 he2TwinVector = globalGeometry.vertexPositions[he2Twin.tipVertex()] - globalGeometry.vertexPositions[he2Twin.tailVertex()];
                he2TwinVector = he2TwinVector.normalize();

                //compute the average gradient across faces on both sides of the singular edges
                Vector3 gradientVector1 = 0.5 * (Vector3{comparisonGrad[he1.face().getIndex()][0], comparisonGrad[he1.face().getIndex()][1], 
                                                    comparisonGrad[he1.face().getIndex()][2]} + Vector3{comparisonGrad[he1Twin.face().getIndex()][0], comparisonGrad[he1Twin.face().getIndex()][1], 
                                                    comparisonGrad[he1Twin.face().getIndex()][2]});
                gradientVector1 = gradientVector1.normalize();
                Vector3 gradientVector2 = 0.5 * (Vector3{comparisonGrad[he2.face().getIndex()][0], comparisonGrad[he2.face().getIndex()][1], 
                                                    comparisonGrad[he2.face().getIndex()][2]} + Vector3{comparisonGrad[he2Twin.face().getIndex()][0], comparisonGrad[he2Twin.face().getIndex()][1], 
                                                    comparisonGrad[he2Twin.face().getIndex()][2]});
                gradientVector2 = gradientVector2.normalize();

                //something going wrong in these constraints? 
                if ((dot(he1Vector, gradientVector1) > dot(he1TwinVector, gradientVector1)) && (dot(he2Vector, gradientVector2) > dot(he2TwinVector, gradientVector2))){
                    //add constraints on singular edges 
                    // model.addConstr(sigma[he1.getIndex()] == sigma[he2.twin().getIndex()]);
                    // model.addConstr(sigma[he1.getIndex()] >= sigma[he1.twin().getIndex()]);
                    // model.addConstr(sigma[he2.twin().getIndex()] <= sigma[he2.getIndex()]);
                    // //add sign constraints for the rest of the halfedges on the face
                    // model.addConstr(sigma[he1.next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he1.next().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he2.twin().next().getIndex()] <= 0.);
                    // model.addConstr(sigma[he2.twin().next().next().getIndex()] <= 0.);
                    
                    //just place the constraints for the positive singularity face for now 
                    model.addConstr(sigma[he1.twin().getIndex()] <= 0.);
                    model.addConstr(sigma[he1.next().getIndex()] >= 0.);
                    model.addConstr(sigma[he1.next().next().getIndex()] >= 0.);

                    //add the constraints for the negative singularity face 
                    model.addConstr(sigma[he2.twin().getIndex()] >= 0.);
                    model.addConstr(sigma[he2.twin().next().getIndex()] <= 0.);
                    model.addConstr(sigma[he2.twin().next().next().getIndex()] <= 0.);

                }
                else if ((dot(he1Vector, gradientVector1) > dot(he1TwinVector, gradientVector1)) && (dot(he2TwinVector, gradientVector2) > dot(he2Vector, gradientVector2))){
                    //add constraints on singular edges
                    // model.addConstr(sigma[he1.getIndex()] == sigma[he2.getIndex()]);
                    // model.addConstr(sigma[he1.getIndex()] >= sigma[he1.twin().getIndex()]);
                    // model.addConstr(sigma[he2.getIndex()] <= sigma[he2.twin().getIndex()]);
                    // //add sign constraints for the rest of the halfedges on the face 
                    // model.addConstr(sigma[he1.next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he1.next().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he2.next().getIndex()] <= 0.);
                    // model.addConstr(sigma[he2.next().next().getIndex()] <= 0.);

                    //just place the constraints for the positive singularity face for now  
                    model.addConstr(sigma[he1.twin().getIndex()] <= 0.);
                    model.addConstr(sigma[he1.next().getIndex()] >= 0.);
                    model.addConstr(sigma[he1.next().next().getIndex()] >= 0.);

                    //add the constraints for the negative singularity face now 
                    model.addConstr(sigma[he2.getIndex()] >= 0.);
                    model.addConstr(sigma[he2.next().getIndex()] <= 0.);
                    model.addConstr(sigma[he2.next().next().getIndex()] <= 0.);


                }
                else if((dot(he1TwinVector, gradientVector1) > dot(he1Vector, gradientVector1)) && (dot(he2Vector, gradientVector2) > dot(he2TwinVector, gradientVector2))){
                    //add constraints on singular edges
                    // model.addConstr(sigma[he1.twin().getIndex()] == sigma[he2.twin().getIndex()]);
                    // model.addConstr(sigma[he1.twin().getIndex()] >= sigma[he1.getIndex()]);
                    // model.addConstr(sigma[he2.twin().getIndex()] <= sigma[he2.getIndex()]);
                    // //add sign constraints for the rest of the halfedges on the face
                    // model.addConstr(sigma[he1.twin().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he1.twin().next().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he2.twin().next().getIndex()] <= 0.);
                    // model.addConstr(sigma[he2.twin().next().next().getIndex()] <= 0.);

                    //just place the constraints for the positive singularity face for now  
                    model.addConstr(sigma[he1.getIndex()] <= 0.);
                    model.addConstr(sigma[he1.twin().next().getIndex()] >= 0.);
                    model.addConstr(sigma[he1.twin().next().next().getIndex()] >= 0.);

                    //add the constraints for the negative singularity face 
                    model.addConstr(sigma[he2.twin().getIndex()] >= 0.);
                    model.addConstr(sigma[he2.twin().next().getIndex()] <= 0.);
                    model.addConstr(sigma[he2.twin().next().next().getIndex()] <= 0.);

                }
                else if((dot(he1TwinVector, gradientVector1) > dot(he1Vector, gradientVector1)) && (dot(he2TwinVector, gradientVector2) > dot(he2Vector, gradientVector2))){
                    //add constraints on singular edges
                    // model.addConstr(sigma[he1.twin().getIndex()] == sigma[he2.getIndex()]);
                    // model.addConstr(sigma[he1.twin().getIndex()] >= sigma[he1.getIndex()]);
                    // model.addConstr(sigma[he2.getIndex()] <= sigma[he2.twin().getIndex()]);
                    // //add sign constraints for the rest of the halfedges on the face
                    // model.addConstr(sigma[he1.twin().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he1.twin().next().next().getIndex()] >= 0.);
                    // model.addConstr(sigma[he2.next().getIndex()] <= 0.);
                    // model.addConstr(sigma[he2.next().next().getIndex()] <= 0.);

                    //just place the constraints for the positive singularity face for now  
                    model.addConstr(sigma[he1.getIndex()] <= 0.);
                    model.addConstr(sigma[he1.twin().next().getIndex()] >= 0.);
                    model.addConstr(sigma[he1.twin().next().next().getIndex()] >= 0);

                    //add the constraints for the negative singularity face now 
                    model.addConstr(sigma[he2.getIndex()] >= 0.);
                    model.addConstr(sigma[he2.next().getIndex()] <= 0.);
                    model.addConstr(sigma[he2.next().next().getIndex()] <= 0.);
                }
            }
        }

        //constraint: add integer variables for homology generators
        for (int i = 0; i < homologyGenerators.size(); i++){
            std::vector<double> path = homologyGenerators[i];
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
            //add the constraints for the homology generators
            //model.addConstr(pathIntegral == period * generatorIntegers[i]);
            model.addConstr(pathIntegral == period * 0.0);  
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
    

        //setting the objective to be min ||\delta sigma - \nabla h||^2
        for (Face f : gluedMesh.faces()){
            //this is comparing to the non-normalized previous iterate 
            //how is it getting closer to unit norm?
            double area = gluedGeometry.faceAreas[f];
            obj +=  area * (((gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0]) * (gradU[f.getIndex()][0] - comparisonGrad[f.getIndex()][0])) 
                    + ((gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1]) * (gradU[f.getIndex()][1] - comparisonGrad[f.getIndex()][1])) 
                    + ((gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2]) * (gradU[f.getIndex()][2] - comparisonGrad[f.getIndex()][2])));  
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize(); 
        objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
        //put the computed one-form into an edge vector
        for (Halfedge he : gluedMesh.halfedges()){
            gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
        }


        // Fix helicing
        if (gbModel.useHelicingCorrection) {
            P("Fixing helicity...");
            EdgeData<int> bigonIndex(gluedMesh); // 0 = regular, 1 = source. No other values are possible!
            int counter = 0;
            for (auto [e1index,e2index] : singularEdges) {
                Edge e1 = gluedMesh.edge(e1index), e2 = gluedMesh.edge(e2index);
                Halfedge he1 = e1.halfedge(), he2 = e2.halfedge();

                if (abs(gluedOneForm[he1.twin()]) > abs(gluedOneForm[he1]))
                    he1 = he1.twin();
                if (abs(gluedOneForm[he2.twin()]) > abs(gluedOneForm[he2]))
                    he2 = he2.twin();

                bigonIndex[e1] = (sgn(gluedOneForm[he1]) == sgn(gluedOneForm[he1.twin()]));
                bigonIndex[e2] = (sgn(gluedOneForm[he2]) == sgn(gluedOneForm[he2.twin()]));

                model.addConstr(sigma[he1.getIndex()] == -sigma[he2.getIndex()]);

                // if (bigonIndex[e1] == 0 && bigonIndex[e2] == 0) { // if both are regular, freeze them to their average
                    // double frozen1form = (abs(gluedOneForm[he1]) + abs(gluedOneForm[he2])) / 2;
                    // model.addConstr(sigma[he1.getIndex()] == frozen1form * sgn(gluedOneForm[he1]));
                    // model.addConstr(sigma[he2.getIndex()] == frozen1form * sgn(gluedOneForm[he2]));
                //     counter++;
                // }
            }
            psMesh.addEdgeScalarQuantity("bigon index", bigonIndex);
            std::cout << "Added " << 2*counter << " equality constraints to ensure non-helicing." << std::endl;

            // Reoptimize
            model.optimize(); 
            objectiveVal = model.get(GRB_DoubleAttr_ObjVal);
            //put the computed one-form into an edge vector
            for (Halfedge he : gluedMesh.halfedges()){
                gluedOneForm[he] = sigma[he.getIndex()].get(GRB_DoubleAttr_X);
            }
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
std::tuple<HalfedgeData<double>, double> computeCourseVirtualSigma(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Model& gbModel, 
                                                        std::map<int, int>& vertexMap, Eigen::SparseMatrix<double, Eigen::RowMajor>& G, polyscope::SurfaceMesh& psMesh){


    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> gluedOneForm(gluedMesh);
    double objectiveVal;
    //query information from the model 
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    gluedGeometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one = gluedGeometry.d1;
    double period = gbModel.getPeriod();
    //std::vector<std::pair<int, int>> singularEdges = gbModel.getSingularEdges();
    std::vector<int> edgeIndices = gbModel.getEdgeIndices();
    //require the face areas
    gluedGeometry.requireFaceAreas();
    //require edge lengths 
    gluedGeometry.requireEdgeLengths();
    //require edge cotan weights
    gluedGeometry.requireEdgeCotanWeights();
    //require the corner angles 
    gluedGeometry.requireCornerAngles();


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

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 1.0);
        
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
        for (Halfedge he : gluedMesh.halfedges()){
            if (edgeIndices[he.edge().getIndex()] != 0){
                model.addConstr(sigma[he.getIndex()] + sigma[he.twin().getIndex()] == edgeIndices[he.edge().getIndex()] * period);
            }
            else{
                model.addConstr(sigma[he.getIndex()] == -1.0 * sigma[he.twin().getIndex()]);
            }
        }

        //set up the objective term
        GRBQuadExpr obj = 0;        
    
        //setting the objective to be min ||\sigma||^2
        for (Halfedge he : gluedMesh.halfedges()){
            if (gluedGeometry.edgeCotanWeights[he.edge()] > 0)
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
    std::cout << "---------------------------------" << std::endl;
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

//compute curl per vertex using the curl discretization from De Goes SIGGRAPH notes (in global setting)
VertexData<double> computeCourseVertexCurl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, 
                                    FaceData<Vector3>& field, std::map<int, std::vector<Halfedge>>& gluedOneRingMap, 
                                    std::vector<int>& gluedEdgeSingularities, HeatMethodDistanceSolver& heatSolver, std::map<int, int>& vertexMap){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    VertexData<double> curl(globalMesh);
    // VertexData<double> distToSourceGlued(gluedMesh, 1.0);
    // VertexData<double> distToSourceGlobal(globalMesh, 1.0);

    // std::vector<Vertex> gluedSourceVerts;
    // for (Edge e : gluedMesh.edges()){
    //     if (std::fabs(gluedEdgeSingularities[e.getIndex()]) == 1 || e.isBoundary()){
    //         gluedSourceVerts.push_back(e.halfedge().tailVertex());
    //         gluedSourceVerts.push_back(e.halfedge().tipVertex());
    //     }
    // }

    // if (gluedSourceVerts.size() > 0){
    //     //compute distance in the glued setting 
    //     distToSourceGlued = heatSolver.computeDistance(gluedSourceVerts);
    //     //convert distance to global setting
    //     distToSourceGlobal = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, distToSourceGlued, vertexMap);
    // }

    globalGeometry.requireFaceAreas();
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        double area = 0.0; // area of the 1-ring of faces
        for (Halfedge he : gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()];
            field[he.face()] = field[he.face()].normalize(); //always normalize the field
            area += globalGeometry.faceArea(he.face());
            sum += dot(hjkVec, field[he.face()]);
        }

        //if you want to do a visualization of curl correction with the Green's function
        //don't multiply by the heat distance as this messes up the curl!!
        if (vi.isBoundary()) curl[vi] = 0.;
        //multiply curl by distance as well
        // else curl[vi] = distToSourceGlobal[vi] * sum / area;
        // else curl[vi] = sum;
        //else curl[vi] = distToSourceGlobal[vi] * sum;
        //else curl[vi] = exp(-pow(distToSourceGlobal[vi], 2)/ (2 * pow(0.107, 2))) * sum;
        else curl[vi] = sum / area;
    }

    return curl;
}

//@debugging
//average vertex curl onto edges
//this should probably be in the glued mesh setting
EdgeData<double> computeVertexAveragedEdgeCurl(VertexPositionGeometry& globalGeometry, VertexData<double>& vertexCurl){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    globalGeometry.requireEdgeLengths();
    EdgeData<double> edgeCurl(globalMesh, 0.);
    for (Edge e : globalMesh.edges()){
        //if (e.isBoundary()) continue;//lens complex at boundary edges is just 0
        double c1 = vertexCurl[e.halfedge().tailVertex()];
        double c2 = vertexCurl[e.halfedge().tipVertex()];
        edgeCurl[e] = (c1 + c2) / 2;
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

//repair a knit graph vertex that's missing connections
//write out GUI to repair a knitgraph vertex
std::vector<int> repairKnitGraphVertex(){

  int id = -1;
  int row_in = -1;
  int row_out = -1;
  int col_in_1 = -1;
  int col_in_2 = -1;
  int col_out_1 = -1;
  int col_out_2 = -1;
  //int isVirtual = 0;

  std::vector<int> result;

  // Register the callback which creates the UI and does the hard work
  auto focusedPopupUI = [&]() {
      
      static bool showWindow = true;
      ImGui::SetNextWindowSize(ImVec2(400, 0), ImGuiCond_Once);
      ImGui::Begin("Repair Knit Graph Vertex", &showWindow);

      ImGui::PushItemWidth(100);
      ImGui::InputInt("Id", &id); 
      ImGui::InputInt("Row in", &row_in);
      ImGui::InputInt("Row out", &row_out);
      ImGui::InputInt("Column in 0", &col_in_1);
      ImGui::InputInt("Column in 1", &col_in_2);
      ImGui::InputInt("Column out 0", &col_out_1);
      ImGui::InputInt("Column out 1", &col_out_2);
      //ImGui::InputInt("isVirtual", &isVirtual);

      ImGui::Separator();
      if (ImGui::Button("Ok")) {
        polyscope::popContext();
      }

      ImGui::End();
  };
  polyscope::pushContext(focusedPopupUI);

  result.push_back(id);
  result.push_back(row_in);
  result.push_back(row_out);
  result.push_back(col_in_1);
  result.push_back(col_in_2);
  result.push_back(col_out_1);
  result.push_back(col_out_2);
  //result.push_back(isVirtual);

  return result;

}


//given a time function over the mesh, extract the saddle vertices from it 
//only works in the global setting for now 
std::vector<Vertex> getSaddleVertices(IntrinsicGeometryInterface& geometry, VertexData<double>& timeFunction){
    
    SurfaceMesh& mesh = geometry.mesh;
    std::vector<Vertex> saddleVertices;
    geometry.requireDECOperators();
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_not;
    d_not = geometry.d0;
    Eigen::VectorXd f(mesh.nVertices());

    //copy time function values into f
    for (Vertex v : mesh.vertices()){
        f(v.getIndex()) = timeFunction[v];
    }

    Eigen::VectorXd d0_f(mesh.nEdges());

    d0_f = d_not * f;
    HalfedgeData<double> oriented_d0_f(mesh);

    for(Edge e : mesh.edges()){
        oriented_d0_f[e.halfedge()] = d0_f(e.getIndex());
        oriented_d0_f[e.halfedge().twin()] = -1.0 * d0_f(e.getIndex());
    }

    VertexData<double> saddleVerticesCopy(mesh);
    
    for (Vertex v : mesh.vertices()){
        saddleVerticesCopy[v] = 0.0;
        //number of sign changes around vertex v
        int ctr = 0;
        for (Halfedge he : v.outgoingHalfedges()){
            assert(he.vertex() == v); // true
            if (sgn(oriented_d0_f[he]) != sgn(oriented_d0_f[he.next().next().twin()])) ctr++;
        }
        if (ctr == 4) saddleVertices.push_back(v);
        
    }

    return saddleVertices;
}

/*
 * Find all of the halfedges in the neighborhood of a a vertex that
 * contain the isoline the harmonic time function of a saddle point.
 */
std::vector<Halfedge> findAllHalfedgesInSaddleLoop(const Vertex &v, const VertexData<double> &timeFunc) {
    std::vector<Halfedge> saddleHalfedges;

    for (const Halfedge &he : v.outgoingHalfedges()) {
        const Halfedge nextHe = he.twin().next();
        const Vertex v1 = he.tipVertex();
        const Vertex v2 = nextHe.tipVertex();
        const double saddleValue = timeFunc[v];
        const double min = std::min(timeFunc[v1], timeFunc[v2]); 
        const double max = std::max(timeFunc[v1], timeFunc[v2]); 
        if (saddleValue > min && saddleValue < max) {
           saddleHalfedges.push_back(nextHe.next());
        }
    }
    return saddleHalfedges;
}

/*
 * This will only ever find one particular halfedge in the saddle loop,
 * I need to perhaps change this to return a list of such halfedges.
 */
Halfedge findFirstHalfedgeInSaddleLoop(const Vertex &v, const VertexData<double>& timeFunc) {
    for (const Halfedge &he : v.outgoingHalfedges()) {
        const Halfedge nextHe = he.twin().next();
        const Vertex v1 = he.tipVertex();
        const Vertex v2 = nextHe.tipVertex();
        const double saddleValue = timeFunc[he.tailVertex()];
        double min = std::min(timeFunc[v1], timeFunc[v2]); 
        double max = std::max(timeFunc[v1], timeFunc[v2]); 
        if (saddleValue > min && saddleValue < max) {
            return nextHe.next();
        }
        continue;
    }
    std::cout << "did not find a first halfedge" << std::endl;
    return v.halfedge(); // fail case
}

/*
 * Returns the halfedge on the **adjacent triangle** from which the isoline.
 * Need to know at least one of the halfedges that the isoline with value "v"
 * goes through. 
 *
 * If this function returns the same halfedge, then there has been an error.
 */
Halfedge findNextHalfedgeInIsoline(double v, const Halfedge &he, const VertexData<double>& timeFunc) {
    const auto heNext = he.next();
    const auto hePrev = he.next().next();

    const double nextMin = std::min(timeFunc[heNext.tailVertex()], 
                                    timeFunc[heNext.tipVertex()]);
    const double nextMax = std::max(timeFunc[heNext.tailVertex()], 
                                    timeFunc[heNext.tipVertex()]);

    const double prevMin = std::min(timeFunc[hePrev.tailVertex()],
                                    timeFunc[hePrev.tipVertex()]);
    const double prevMax = std::max(timeFunc[hePrev.tailVertex()], 
                                    timeFunc[hePrev.tipVertex()]);

    if (v > nextMin && v < nextMax) {
        return heNext.twin();
    }
    else if (v > prevMin && v < prevMax) {
        return hePrev.twin();
    }
    else{
        // this is a fail case
        return he;
    }
}

std::vector<int> findSaddleLoop(IntrinsicGeometryInterface& geometry, const VertexData<double> &saddleVertices, const VertexData<double>& timeFunc) {

    std::vector<int> edgeLoop(geometry.mesh.nEdges());
    for (const Edge &edge : geometry.mesh.edges()) {
        edgeLoop[edge.getIndex()] = 0;
    }
    // loop over all saddle vertices
    for (const Vertex &v : geometry.mesh.vertices()) {
        // skip all vertices that are not saddle points
        if (saddleVertices[v] != 1.0) {continue;}
        std::cout << "found the saddle vertex" << std::endl;
        const double saddleValue = timeFunc[v];
        std::cout << "saddle vertex has a harmonic value of: " << saddleValue << std::endl;
        const auto allHalfedges = findAllHalfedgesInSaddleLoop(v, timeFunc); 
        // this for loop will create duplicate edge loops
        for (const Halfedge &rootHe : allHalfedges) {
            // in order for the search to go correctly, the first halfedge 
            // cannot be on the first face
            Halfedge currHe = rootHe.twin();
            edgeLoop[currHe.edge().getIndex()] = 1;
            std::cout << "searching for next halfedge" << std::endl;
            while (true) {
                std::cout << "current halfedge: " << currHe << std::endl;
                auto nextHe = findNextHalfedgeInIsoline(saddleValue, currHe, timeFunc);
                std::cout << "next halfedge: " << nextHe << std::endl;
                if (currHe == nextHe) { break;}
                edgeLoop[currHe.edge().getIndex()] = 1;
                currHe = nextHe;
            }
        }
    }
    return edgeLoop;
}
/*
 * Finds the triangle strip that contains the isoline of the value 'saddleValue'
 * that starts at 'rootHe'
 */
std::vector<Face> findTriangleStripSaddleLoop(IntrinsicGeometryInterface& geometry, const VertexData<double> &timeFunc, const Halfedge &rootHe, double saddleValue) {
    // list of all faces in the triangle strip
    std::vector<Face> triangleStrip;

    // NOTE: I added this line to make sure the first and last triangles match up
    triangleStrip.push_back(rootHe.face());
    // this assume that the rootHe is a halfedge across from the saddle vertex
    Halfedge currHe = rootHe.twin();
    triangleStrip.push_back(currHe.face());
    while(true) {
        auto nextHe = findNextHalfedgeInIsoline(saddleValue, currHe, timeFunc);
        if (currHe == nextHe) {
            break;
        }
        triangleStrip.push_back(nextHe.face());
        currHe = nextHe;
    }
    return triangleStrip;
}


/*
 * Find the unique edges of a triangle strip.
 *
 */
std::vector<Edge> findUniqueEdgesInTriangleStrip(const std::vector<Face> &triangleStrip) {
    std::vector<Edge> allEdges;
    std::vector<Edge> uniqueEdges;
    // go through all of the triangles in the strip
    for (const Face &f : triangleStrip) {
        // go through all of the edges of the triangle
        for (const Edge &e : f.adjacentEdges()) {
            allEdges.push_back(e);
        }
    }
    for (const Edge &e1 : allEdges) {
        int uniqueCount = 0;
        for (const Edge &e2 : allEdges) {
            if (e1 == e2) {uniqueCount++;}
        }
        if (uniqueCount == 1) {
            uniqueEdges.push_back(e1);
        }
    }
    return uniqueEdges;
}

/*
 * Given a list of edges that contains two loops, chooses one of the loops
 * and disregards the edges in the other loop.
 *
 */
std::vector<Edge> chooseEdgeLoop(const std::vector<Edge> &doubleLoop, 
                                 const std::vector<Face> &triangleStrip) {

    // construct the halfedges bounding the triangle strip
    std::vector<Halfedge> doubleLoopHalfedges;
    for (const Edge &e : doubleLoop) {
       for (const Face &f : triangleStrip) {
           if (e.halfedge().face() == f) {
               doubleLoopHalfedges.push_back(e.halfedge());
           }
           else if (e.halfedge().twin().face() == f) {
               doubleLoopHalfedges.push_back(e.halfedge().twin());
           }
       }
    }

    // visited map
    std::unordered_map<Halfedge, bool> visited;
    for (const Halfedge &he : doubleLoopHalfedges) {
        visited[he] = false;
    }
    std::unordered_map<Vertex, bool> visitedVertex;
    for (const Halfedge &he : doubleLoopHalfedges) {
        visitedVertex[he.tipVertex()] = false;
        visitedVertex[he.tailVertex()] = false;
    }

    // choose a particular halfedge to start searching for the root
    std::vector<Halfedge> loopHalfedges;
    const Halfedge rootHe = doubleLoopHalfedges[0];
    Halfedge currHe = rootHe;
    visitedVertex[rootHe.tailVertex()] = true;
    loopHalfedges.push_back(currHe);


    // find the loop
    do {
        for (const Halfedge &he : doubleLoopHalfedges) {
            if (currHe.tipVertex() == he.tailVertex() && !visited[he]) {
                loopHalfedges.push_back(he);
                visited[he] = true;
                currHe = he;
                visitedVertex[currHe.tailVertex()] = true;
                break;
            }
        }
    } while (!visitedVertex[currHe.tipVertex()]);

    std::vector<Edge> edgeLoop;
    for (const Halfedge &he : loopHalfedges) {
        edgeLoop.push_back(he.edge());
    }

    return edgeLoop;
}

std::vector<Halfedge> chooseHalfEdgeLoop(const std::vector<Edge> &doubleLoop, 
                                 const std::vector<Face> &triangleStrip) {

    // construct the halfedges bounding the triangle strip
    std::vector<Halfedge> doubleLoopHalfedges;
    for (const Edge &e : doubleLoop) {
       for (const Face &f : triangleStrip) {
           if (e.halfedge().face() == f) {
               doubleLoopHalfedges.push_back(e.halfedge());
           }
           else if (e.halfedge().twin().face() == f) {
               doubleLoopHalfedges.push_back(e.halfedge().twin());
           }
       }
    }

    // visited map
    std::unordered_map<Halfedge, bool> visited;
    for (const Halfedge &he : doubleLoopHalfedges) {
        visited[he] = false;
    }

    // choose a particular halfedge to start searching for the root
    std::vector<Halfedge> loopHalfedges;
    const Halfedge rootHe = doubleLoopHalfedges[0];
    Halfedge currHe = rootHe;
    loopHalfedges.push_back(currHe);

    // find the loop
    do {
        for (const Halfedge &he : doubleLoopHalfedges) {
            if (currHe.tipVertex() == he.tailVertex() && !visited[he]) {
                loopHalfedges.push_back(he);
                visited[he] = true;
                currHe = he;
                break;
            }
        }
    } while (currHe != rootHe);
    return loopHalfedges;
}

void registerCurveNetworkFromEdges(VertexPositionGeometry& geometry, 
        const std::vector<Edge> &edges, std::string &name) {
    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeIndices;

    size_t nodeCounter = 0;
    for (const auto &e : edges) {
        const auto p1 = geometry.vertexPositions[e.firstVertex()];
        const auto p2 = geometry.vertexPositions[e.secondVertex()];
        positions.push_back(p1);
        positions.push_back(p2);
        edgeIndices.push_back({nodeCounter, nodeCounter+1});
        nodeCounter += 2;
    }
    polyscope::registerCurveNetwork(name, positions, edgeIndices)->setEnabled(false);
}

void visualizeAllEdgeLoops(VertexPositionGeometry& geometry, std::vector<Vertex> &saddleVertices, const VertexData<double> &timeFunc) {
    //const std::vector<Vertex> saddleVertices = extractSaddlePoints(geometry, saddleVerticesData);
    std::cout << "extracted Saddle Points" << std::endl;
    // for each saddle vertex, trace the time function isolines
    int loopCounter = 0;
    // iterate over the saddle points
    for (const Vertex &v : saddleVertices) {
        std::vector<std::vector<Face>> allTriangleStripsPerVertex;
        const double saddleValue = timeFunc[v];
        // find all of the halfedges that contain the isoline
        std::vector<Halfedge> isolineHalfedges = findAllHalfedgesInSaddleLoop(v, timeFunc);
        for (const Halfedge &rootHe : isolineHalfedges) {
            const std::vector<Face> triangleStrip = findTriangleStripSaddleLoop(geometry, timeFunc, rootHe, saddleValue);
           // make sure not adding duplicate paths
            bool stripAlreadyExists = false;
            for (const auto &strip : allTriangleStripsPerVertex) {
                if (triangleStrip.front() == strip.back() || 
                    triangleStrip.back() == strip.front()) {
                    stripAlreadyExists = true;
                    break;
                }
            }
            if (!stripAlreadyExists) {
                allTriangleStripsPerVertex.push_back(triangleStrip);
            }
        }
        for (const auto &strip : allTriangleStripsPerVertex) {
            auto uniqueEdges = findUniqueEdgesInTriangleStrip(strip);
            auto edgeLoop = chooseEdgeLoop(uniqueEdges, strip);
            std::string network = "edge loop" + std::to_string(loopCounter++);
            registerCurveNetworkFromEdges(geometry, edgeLoop, network); 
        }
        
    }
}

//the only function I ever need to care about
std::vector<std::vector<double>> findAllSaddleLoops(VertexPositionGeometry& geometry, const std::vector<Vertex> &saddleVertices, const VertexData<double>& timeFunc) {

    std::vector<std::vector<double>>  allSaddleLoops;
    
    int loopCounter = 0;
    // loop over all saddle vertices
    for (const Vertex &v : saddleVertices) {
        // skip all vertices that are not saddle points
        //if (saddleVertices[v] != 1.0) {continue;}
        // find the time function value at the saddle point
        const double saddleValue = timeFunc[v];
        // find all of the neighboring halfedge through which the isoline exits
        const auto allNeighborHalfedges = findAllHalfedgesInSaddleLoop(v, timeFunc); 

        // the set of all triangle strips per saddle vertex
        std::vector<std::vector<Face>> allTriangleStripsPerVertex;
        // iterate through all of the neighboring halfedges 
        for (const Halfedge &rootHe : allNeighborHalfedges) {
            const auto triangleStrip = findTriangleStripSaddleLoop(geometry, timeFunc, rootHe, saddleValue); 

            // this code inserts a triangle strip if it is unique
            bool stripAlreadyExists = false;
            for (const auto &strip : allTriangleStripsPerVertex) {
                if (triangleStrip.front() == strip.back() || 
                        triangleStrip.back() == strip.front()) {
                    stripAlreadyExists = true;
                    break;
                }
            }
            if (!stripAlreadyExists) {
                allTriangleStripsPerVertex.push_back(triangleStrip);
            }
        }
        int i = 0;
        // by this point all of the unique triangle strips for this particular vertex will be added to the list  
        for (const auto &strip : allTriangleStripsPerVertex) {
            std::vector<double> saddleLoop(geometry.mesh.nEdges());
            std::fill(saddleLoop.begin(), saddleLoop.end(), 0.0);
            auto uniqueEdges = findUniqueEdgesInTriangleStrip(strip);
            auto edgeLoop = chooseEdgeLoop(uniqueEdges, strip);
            auto halfedgeLoop = chooseHalfEdgeLoop(uniqueEdges, strip);
            // construct the output here
            for (const Halfedge &he : halfedgeLoop) {
                if (he.edge().halfedge() == he) {
                    saddleLoop[he.edge().getIndex()] = 1.0;
                }
                else {
                    saddleLoop[he.edge().getIndex()] = -1.0;
                }
            }
            allSaddleLoops.push_back(saddleLoop);
            // visualization (uncomment if you don't care)
            std::string network = "edge loop" + std::to_string(loopCounter++);
            registerCurveNetworkFromEdges(geometry, edgeLoop, network); 
            //std::cout << "the size of loop " << i++ << "is: " << edgeLoop.size() << std::endl;
            //std::cout << "the size of halfedge loop " << i++ << "is: " << halfedgeLoop.size() << std::endl;
        }
    }

    return allSaddleLoops;
}

//the only function I ever need to care about
//overloaded to handle intrinsic geometry
std::vector<std::vector<double>> findAllSaddleLoops(IntrinsicGeometryInterface& geometry, const std::vector<Vertex> &saddleVertices, const VertexData<double>& timeFunc) {

    std::vector<std::vector<double>>  allSaddleLoops;
    
    int loopCounter = 0;
    // loop over all saddle vertices
    for (const Vertex &v : saddleVertices) {
        // skip all vertices that are not saddle points
        //if (saddleVertices[v] != 1.0) {continue;}
        // find the time function value at the saddle point
        const double saddleValue = timeFunc[v];
        // find all of the neighboring halfedge through which the isoline exits
        const auto allNeighborHalfedges = findAllHalfedgesInSaddleLoop(v, timeFunc); 

        // the set of all triangle strips per saddle vertex
        std::vector<std::vector<Face>> allTriangleStripsPerVertex;
        // iterate through all of the neighboring halfedges 
        for (const Halfedge &rootHe : allNeighborHalfedges) {
            const auto triangleStrip = findTriangleStripSaddleLoop(geometry, timeFunc, rootHe, saddleValue); 

            // this code inserts a triangle strip if it is unique
            bool stripAlreadyExists = false;
            for (const auto &strip : allTriangleStripsPerVertex) {
                if (triangleStrip.front() == strip.back() || 
                        triangleStrip.back() == strip.front()) {
                    stripAlreadyExists = true;
                    break;
                }
            }
            if (!stripAlreadyExists) {
                allTriangleStripsPerVertex.push_back(triangleStrip);
            }
        }
        int i = 0;
        // by this point all of the unique triangle strips for this particular vertex will be added to the list  
        for (const auto &strip : allTriangleStripsPerVertex) {
            std::vector<double> saddleLoop(geometry.mesh.nEdges());
            std::fill(saddleLoop.begin(), saddleLoop.end(), 0.0);
            auto uniqueEdges = findUniqueEdgesInTriangleStrip(strip);
            auto edgeLoop = chooseEdgeLoop(uniqueEdges, strip);
            auto halfedgeLoop = chooseHalfEdgeLoop(uniqueEdges, strip);
            // construct the output here
            for (const Halfedge &he : halfedgeLoop) {
                if (he.edge().halfedge() == he) {
                    saddleLoop[he.edge().getIndex()] = 1.0;
                }
                else {
                    saddleLoop[he.edge().getIndex()] = -1.0;
                }
            }
            allSaddleLoops.push_back(saddleLoop);
            // visualization (uncomment if you don't care)
            //std::string network = "edge loop" + std::to_string(loopCounter++);
            //registerCurveNetworkFromEdges(geometry, edgeLoop, network); 
            //std::cout << "the size of loop " << i++ << "is: " << edgeLoop.size() << std::endl;
            //std::cout << "the size of halfedge loop " << i++ << "is: " << halfedgeLoop.size() << std::endl;
        }
    }

    return allSaddleLoops;
}

