#include "knitting_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//solve Laplace equation over an input mesh 
VertexData<double> solveLaplace(VertexPositionGeometry& geometry, std::vector<Vertex>& zeroVertices, std::vector<Vertex>& oneVertices){

    SurfaceMesh& mesh = geometry.mesh;

    Eigen::VectorXi knittingStartVertices(zeroVertices.size());
    Eigen::VectorXi knittingEndVertices(oneVertices.size());

    for (int i = 0; i < knittingStartVertices.size(); i++){
        knittingStartVertices(i) = zeroVertices[i].getIndex();
    }
    for (int i = 0; i < knittingEndVertices.size(); i++){
        knittingEndVertices(i) = oneVertices[i].getIndex();
    }

    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    std::tie(V, F) = getVertexPositionsandFaceLists(geometry);


    //set up eigen code for solving laplace equation
    Eigen::VectorXi zero_vertices(zeroVertices.size());
    Eigen::VectorXi one_vertices(oneVertices.size());

    for (int i = 0; i < zero_vertices.size(); i++){
        zero_vertices(i) = zeroVertices[i].getIndex();
    }
    for (int i = 0; i < one_vertices.size(); i++){
        one_vertices(i) = oneVertices[i].getIndex();
    }

    //join the vertices that will serve as our boundary condition
    Eigen::VectorXi joined(zero_vertices.size() + one_vertices.size());
    joined << zero_vertices, one_vertices;

    //List of all vertex indices
    Eigen::VectorXi all, in, IA;
    igl::colon<int>(0, V.rows()-1, all);
    //List of interior vertices
    igl::setdiff(all, joined, in, IA);
    // Construct and slice up Laplacian
    Eigen::SparseMatrix<double> L,L_in_in,L_in_b;
    igl::cotmatrix(V,F,L);
    igl::slice(L,in,in,L_in_in);
    igl::slice(L,in,joined,L_in_b);
    //set boundary conditions
    Eigen::VectorXd bc;
    Eigen::VectorXd Z(V.rows());
    
    //set boundary values
    //setting the zero vertices
    for ( int i = 0 ; i < zero_vertices.size(); i++){
        Z(zero_vertices(i)) = 0.0;
    }
    //setting the one vertices
    for ( int i = 0 ; i < one_vertices.size(); i++){
      	Z(one_vertices(i)) = 1.0;
    }

    igl::slice(Z,joined,bc); //makes bc a #b  x 1 matrix
    // Solve PDE for any given boundary condition
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(-L_in_in);
    Eigen::VectorXd Z_in = solver.solve(L_in_b*bc);
    // slice into solution
    igl::slice_into(Z_in,in,Z);
    VertexData<double> timeFunction(mesh);

    for (Vertex v : mesh.vertices()){
        timeFunction[v] = Z(v.getIndex());
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
EdgeData<double> computeOneForm(VertexPositionGeometry& geometry, Model& gbModel){

    SurfaceMesh& mesh = geometry.mesh;
    geometry.requireDECOperators();

    Eigen::SparseMatrix<double, Eigen::RowMajor> d_one;
    d_one = geometry.d1;

    EdgeData<double> oneForm(mesh);

    std::vector<double> omega = gbModel.getMatchingTerms();
    std::vector<int> bdyEdges = gbModel.getBdyEdges();
    double period = gbModel.getPeriod();

    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //set the timeout
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 120);

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
      		if (edge.halfedge().orientation()){
                d0_f_avg[edge] = dot(faceGradients[edge.halfedge().face()], geometry.vertexPositions[edge.halfedge().tipVertex()] - 
                                                                                geometry.vertexPositions[edge.halfedge().tailVertex()]);
            }
      		else{//should never actually hit this case
                d0_f_avg[edge] = dot(faceGradients[edge.halfedge().face()], geometry.vertexPositions[edge.halfedge().tailVertex()] - 
                                                                                geometry.vertexPositions[edge.halfedge().tipVertex()]);
      		}
      	}
    }

    return d0_f_avg;
}

