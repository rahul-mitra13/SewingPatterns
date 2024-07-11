#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//compute time function using a vector of pairs of vertex mappings instead of the map because we miss stitches then 
VertexData<double> computeTimeFunction(VertexPositionGeometry& geometry, std::vector<std::pair<int,int>>& vertexMappingsPairs,  globalBoundaryConditions& boundaryConditions){

    SurfaceMesh& mesh = geometry.mesh;
    int numVertices = mesh.nVertices(); 
    int numStitches = vertexMappingsPairs.size();
    //require the cotan edge weights 
    geometry.requireEdgeCotanWeights();
    //store mappings between index in the original mesh and index in the Laplacian matrix
    std::map<int, int> originalIndexToLaplacianMatrixIndex;
    //store mappings between index in the Laplacian matrix to index in the orignal mesh 
    std::map<int, int> laplacianMatrixIndexToOriginalIndex;
    //all the pairs that have been seen so far
    std::vector<std::pair<int, int>> seenPairs;
    //number of "unique vertices" i.e., only consider one vertex per stitch
    int numUniqueVertices = 0;
   
    for (Vertex v : mesh.vertices()){
        size_t iV = v.getIndex();
        //iterate over the mappings 
        for (auto p : vertexMappingsPairs){
            if (p.first == iV || p.second == iV){
                if (originalIndexToLaplacianMatrixIndex.find(p.first) == originalIndexToLaplacianMatrixIndex.end()
                && originalIndexToLaplacianMatrixIndex.find(p.second) == originalIndexToLaplacianMatrixIndex.end()){
                    originalIndexToLaplacianMatrixIndex.insert({p.first, numUniqueVertices});
                    originalIndexToLaplacianMatrixIndex.insert({p.second, numUniqueVertices});
                    numUniqueVertices++;
                }
                if (originalIndexToLaplacianMatrixIndex.find(p.first) != originalIndexToLaplacianMatrixIndex.end()&&
                    originalIndexToLaplacianMatrixIndex.find(p.second) == originalIndexToLaplacianMatrixIndex.end()){
                    originalIndexToLaplacianMatrixIndex.insert({p.second, originalIndexToLaplacianMatrixIndex.at(p.first)});
                }
                if (originalIndexToLaplacianMatrixIndex.find(p.second) != originalIndexToLaplacianMatrixIndex.end() &&
                    originalIndexToLaplacianMatrixIndex.find(p.first) == originalIndexToLaplacianMatrixIndex.end()){
                    originalIndexToLaplacianMatrixIndex.insert({p.first, originalIndexToLaplacianMatrixIndex.at(p.second)});
                }
            }
        }
        if (originalIndexToLaplacianMatrixIndex.find(iV) == originalIndexToLaplacianMatrixIndex.end()){
            originalIndexToLaplacianMatrixIndex.insert({iV, numUniqueVertices});
            numUniqueVertices++;
        }
    }

    Eigen::SparseMatrix<double> L(originalIndexToLaplacianMatrixIndex.size(), originalIndexToLaplacianMatrixIndex.size());
    std::vector<Eigen::Triplet<double>> tripletList;
    //keep a set of indices you've already populated in the Laplacian 
    std::set<int> setIndices;//set of vertex indices we've already set in the laplacian

    for (Vertex v : mesh.vertices()){
        if (std::find(setIndices.begin(), setIndices.end(), originalIndexToLaplacianMatrixIndex.at(v.getIndex())) != setIndices.end()) continue;//we've handled this already
        double L_diag = 0.0;
        //iterate over the 1-ring of the vertex 
        //iterate over the one-ring of the vertex 
        for (Halfedge he : v.outgoingHalfedges()){
            //off diagonal entries 
            tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(he.tipVertex().getIndex()),
            -geometry.edgeCotanWeights[he.edge()]);
            setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));

            for (auto p : vertexMappingsPairs){
                if (p.first == v.getIndex()){//grab the contributes from the second in the pair
                    Vertex mappedVertex = mesh.vertex(p.second);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), 
                                        originalIndexToLaplacianMatrixIndex.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                        setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
                    }
                }
                if (p.second == v.getIndex()){
                    Vertex mappedVertex = mesh.vertex(p.first);
                    //iterate over the 1-ring of the mapped vertex 
                    for (Halfedge mappedVertexHalfedge : mappedVertex.outgoingHalfedges()){
                        tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), 
                                        originalIndexToLaplacianMatrixIndex.at(mappedVertexHalfedge.tipVertex().getIndex()),
                                        -geometry.edgeCotanWeights[mappedVertexHalfedge.edge()]);
                        setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
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
        tripletList.emplace_back(originalIndexToLaplacianMatrixIndex.at(v.getIndex()), originalIndexToLaplacianMatrixIndex.at(v.getIndex()), L_diag);
        setIndices.insert(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
    }

    L.setFromTriplets(tripletList.begin(), tripletList.end());

    //force boundary conditions
    Eigen::VectorXd b = Eigen::VectorXd::Zero(numUniqueVertices);
    
    for (Vertex v : boundaryConditions.courseStartBoundaryVertices){
        int updatedIndex = originalIndexToLaplacianMatrixIndex.at(v.getIndex());
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 0.0;
    }

    for (Vertex v : boundaryConditions.courseEndBoundaryVertices){
        int updatedIndex = originalIndexToLaplacianMatrixIndex.at(v.getIndex());
        L.row(updatedIndex) *= 0.0;
        L.coeffRef(updatedIndex, updatedIndex) = 1.0;
        b(updatedIndex) = 1.0;
    }
    
    Eigen::SparseQR<SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //Eigen::SparseLU<SparseMatrix<double>> solver;
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
        timeFunction[v] = u(originalIndexToLaplacianMatrixIndex.at(v.getIndex()));
    }
    return timeFunction;
}


//compute the gradient of a function defined as a scalar over vertices
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
    Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G*U).eval().data(),F.rows(),3);

    //return as face data
    for (Face f : mesh.faces()){
        faceGradients[f] = Vector3{GU(f.getIndex(), 0), GU(f.getIndex(), 1), GU(f.getIndex(), 2)};
    }

    return faceGradients;

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
        std::complex<double> r(cos (2.0 * alpha), sin(1.0 * alpha));
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
EdgeData<double> computeOneForm(VertexPositionGeometry& geometry, Model& gbModel, polyscope::SurfaceMesh& psMesh){

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
    //std::cout << "d1 * d0: " << d_one * d_zero << std::endl;

    //some testing code for sparse matrix updates 
    // Eigen::SparseMatrix<double, Eigen::RowMajor> d0test(6, 6);
    // d0test.coeffRef(0, 0) = 1;
    // d0test.coeffRef(0, 1) = -1;
    // d0test.coeffRef(1, 1) = 1;
    // d0test.coeffRef(1, 2) = -1;
    // d0test.coeffRef(2, 0) = 1;
    // d0test.coeffRef(2, 2) = -1;
    // d0test.coeffRef(3, 3) = 1;
    // d0test.coeffRef(3, 4) = -1;
    // d0test.coeffRef(4, 4) = 1;
    // d0test.coeffRef(4, 5) = -1;
    // d0test.coeffRef(5, 3) = -1;
    // d0test.coeffRef(5, 5) = 1;
    // Eigen::SparseMatrix<double, Eigen::RowMajor> d1test(2, 6);
    // d1test.coeffRef(0, 0) = -1;
    // d1test.coeffRef(0, 1) = -1;
    // d1test.coeffRef(0, 2) = 1;
    // d1test.coeffRef(1, 3) = -1;
    // d1test.coeffRef(1, 4) = -1;
    // d1test.coeffRef(1, 5) = -1;
    // std::cout << "d0test: " << d0test << std::endl;
    // std::cout << "d1test: " << d1test << std::endl;
    // Eigen::SparseMatrix<double, Eigen::ColMajor> d1testColMajor;
    // d1testColMajor = d1test;
    // for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(d0test, 3); it; ++it){
    //     it.valueRef() = -it.value();
    // }

    // for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(d1testColMajor, 3); it; ++it){
    //     it.valueRef() = -it.value();
    // }
    // d1test = d1testColMajor;
    // std::cout << "d0test after update: " << d0test << std::endl;
    // std::cout << "d1test after update: " << d1test << std::endl;
    // std::cout << "test result: " << d1test * d0test << std::endl;
    
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);

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
    }
    catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    psMesh.addFaceScalarQuantity("Integrated 1-form value", integratedOneForm);
    return oneForm;
}

//compute \omega that is the 1-form we're trying to match over each edge 
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3> faceGradients){
    
    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireFaceNormals();
    EdgeData<double> d0_f_avg(mesh);


    //if we're computing the matching 1-form in the wale direction, rotate all the gradients 
    if (direction == 1){
        for (Face f : mesh.faces()){
            faceGradients[f] = faceGradients[f].rotateAround(geometry.faceNormals[f], PI/2);
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

//comput matching 1-form while taking into account "stitched together" edges
EdgeData<double> computeMatchingOneForm(VertexPositionGeometry& geometry, int direction, FaceData<Vector3> faceGradients, std::vector<std::pair<int, int>>& edgeMappingsPairs){

    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireFaceNormals();
    EdgeData<double> omega(mesh);
    //if we're computing the matching 1-form in the wale direction, rotate all the gradients 
    if (direction == 1){
        for (Face f : mesh.faces()){
            faceGradients[f] = faceGradients[f].rotateAround(geometry.faceNormals[f], PI/2.);
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
            faceGradients[e.halfedge().face()] = faceGradients[e.halfedge().face()].normalize();
            faceGradients[e.halfedge().twin().face()] = faceGradients[e.halfedge().twin().face()].normalize();
            omega[e] = 0.5 * dot((faceGradients[e.halfedge().face()] + faceGradients[e.halfedge().twin().face()]),
                                    geometry.vertexPositions[e.halfedge().tipVertex()] - geometry.vertexPositions[e.halfedge().tailVertex()]);
            seenEdges[e.getIndex()] = true;
        }
        else{//found a boundary halfedge
            if (edgeMap.find(e.getIndex()) != edgeMap.end()){//found a stitched together edge
                numStitchedEdges++;
                Vector3 faceGradients1 = faceGradients[e.halfedge().face()].normalize();
                Vector3 faceGradients2 = faceGradients[mesh.edge(edgeMap.at(e.getIndex())).halfedge().face()].normalize();
                //take the average direction vector? 
                Vector3 e1 = geometry.vertexPositions[e.halfedge().tipVertex()] - geometry.vertexPositions[e.halfedge().tailVertex()];
                Vector3 e2 = geometry.vertexPositions[mesh.edge(edgeMap.at(e.getIndex())).halfedge().tipVertex()] 
                                            - geometry.vertexPositions[mesh.edge(edgeMap.at(e.getIndex())).halfedge().tailVertex()];
                Vector3 avgVector = (e1 + e2)/2.;
                //just pick the original edge as the "canonical" direction in the global mesh
                //I'm not really sure the polyscope Whitney interpolation scheme is the best way to visualize these
                omega[e] = 0.5 * dot((faceGradients1 + faceGradients2), e1);
                omega[mesh.edge(edgeMap.at(e.getIndex()))] = 0.5 * dot((faceGradients1 + faceGradients2), e1);
                seenEdges[e.getIndex()] = true;
                seenEdges[edgeMap.at(e.getIndex())] = true;
            }
            else{//found a boundary edge that's not stitched to anything
                faceGradients[e.halfedge().face()] = faceGradients[e.halfedge().face()].normalize();
                omega[e] = dot(faceGradients[e.halfedge().face()], geometry.vertexPositions[e.halfedge().tipVertex()] - 
                                                                                geometry.vertexPositions[e.halfedge().tailVertex()]);
                seenEdges[e.getIndex()] = true;
            }
        }
    }
    return omega;
}

