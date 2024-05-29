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

    for (Vertex v : mesh.vertices()){
        Vector3 pos = geometry.vertexPositions[v];
        V(v.getIndex(), 0) = pos.x;
        V(v.getIndex(), 1) = pos.y;
        V(v.getIndex(), 2) = pos.z;
    }

    std::vector<std::vector<size_t>> faceVertexIndices = mesh.getFaceVertexList();

    for (int i = 0; i < faceVertexIndices.size(); i++){
        F(i, 0) = faceVertexIndices[i][0];
        F(i, 1) = faceVertexIndices[i][1];
        F(i, 2) = faceVertexIndices[i][2];
    }


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

    for (Vertex v : mesh.vertices()){
        Vector3 pos = geometry.vertexPositions[v];
        V(v.getIndex(), 0) = pos.x;
        V(v.getIndex(), 1) = pos.y;
        V(v.getIndex(), 2) = pos.z;
    }

    std::vector<std::vector<size_t>> faceVertexIndices = mesh.getFaceVertexList();

    for (int i = 0; i < faceVertexIndices.size(); i++){
        F(i, 0) = faceVertexIndices[i][0];
        F(i, 1) = faceVertexIndices[i][1];
        F(i, 2) = faceVertexIndices[i][2];
    }

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
VertexData<Vector3> computeVertexValuedField(VertexPositionGeometry& geometry, VertexData<double>& vertexScalarFunction){

    SurfaceMesh& mesh = geometry.mesh;
    VertexData<Vector3> vertexValuedField(mesh);

    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    for (Vertex v : mesh.vertices()){
        Vector3 pos = geometry.vertexPositions[v];
        V(v.getIndex(), 0) = pos.x;
        V(v.getIndex(), 1) = pos.y;
        V(v.getIndex(), 2) = pos.z;
    }

    std::vector<std::vector<size_t>> faceVertexIndices = mesh.getFaceVertexList();

    for (int i = 0; i < faceVertexIndices.size(); i++){
        F(i, 0) = faceVertexIndices[i][0];
        F(i, 1) = faceVertexIndices[i][1];
        F(i, 2) = faceVertexIndices[i][2];
    }

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
    }

    return vertexValuedField;

}