#include "experiments.h"

std::tuple<HalfedgeData<double>, VertexData<double>> strategy1Impl(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& gluedTimeFunction,
                                                                    FaceData<Vector3>& globalFaceGradients, std::map<int, std::vector<Halfedge>>& gluedOneRingMap, polyscope::SurfaceMesh& psMesh,
                                                                    std::map<int, int>& globalToGluedVertexMap, double period){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    HalfedgeData<double> sigmaTilde(globalMesh);
    VertexData<double> vertexSingularities(globalMesh, 0.0);
    //singular vertex indices 
    std::vector<std::pair<int, int>> singularVertices;

    int maxPairs = 1;
    int numPairs = 0;

    //compute curl per vertex of gradient field
    VertexData<double> curl = computeVertexCurl(globalGeometry, gluedGeometry, globalFaceGradients, gluedOneRingMap);

    //find max/min vertex over entire mesh
    std::pair<int, int> singVertexPair = findVertexSingularityPair(globalGeometry, gluedGeometry, curl,
                                        gluedTimeFunction, psMesh, globalToGluedVertexMap, 0.0, true);
    vertexSingularities[globalMesh.vertex(singVertexPair.first)] = 1.0;
    vertexSingularities[globalMesh.vertex(singVertexPair.second)] = -1.0;
    psMesh.addVertexScalarQuantity("vertex singularities", vertexSingularities);
    
    //make the model
    //we will solve it in the glued mesh setting 
    Model model;
    std::vector<std::array<double, 3>> faceGradients(gluedMesh.nFaces());
    for (Face f : globalMesh.faces()){
        faceGradients[f.getIndex()] = std::array{globalFaceGradients[f][0], globalFaceGradients[f][1], globalFaceGradients[f][2]};
    }
    model.setFaceGradients(faceGradients);
    model.setPeriod(period);

    return std::tie(sigmaTilde, vertexSingularities);

}

//compute the per-face gradient of a 1-form (in global setting)
FaceData<Vector3> computeOneFormFaceGrad(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& sigmaTilde){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    FaceData<Vector3> gradients(globalMesh);

    globalGeometry.requireFaceNormals();
    for (Face f : globalMesh.faces()){
        Eigen::MatrixXd faceSystem(3, 3);
        Eigen::VectorXd faceSigmas(3);
        Vector3 hijVec = globalGeometry.vertexPositions[f.halfedge().tipVertex()] - globalGeometry.vertexPositions[f.halfedge().tailVertex()];
        Vector3 hjkVec = globalGeometry.vertexPositions[f.halfedge().next().tipVertex()] - 
                            globalGeometry.vertexPositions[f.halfedge().next().tailVertex()];
        faceSigmas(0) = sigmaTilde[f.halfedge()];
        faceSigmas(1) = sigmaTilde[f.halfedge().next()];
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

//find max/min curl vertex for a given isoline (on global setting)
std::pair<int, int> findVertexSingularityPair(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, VertexData<double>& curl,
                                            VertexData<double>& gluedTimeFunction, polyscope::SurfaceMesh& psMesh, std::map<int, int>& globalToGluedVertexMap,
                                            double isoVal, bool useAllVertices){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh; 
    VertexData<double> globalTimeFunction = convertGluedToGlobalVertexFunction(globalGeometry, gluedGeometry, gluedTimeFunction, globalToGluedVertexMap);
    int maxVertex, minVertex = -1;
    double maxCurl = -DBL_MAX;
    double minCurl = DBL_MAX;
    
    if (useAllVertices){
        //find max curl over entire mesh 
        //find vertex with max  curl
        for (Vertex v : globalMesh.vertices()){
            if (curl[v] > maxCurl){
                maxCurl = curl[v];
                maxVertex = v.getIndex();
            }
        }
        double isoVal = globalTimeFunction[globalMesh.vertex(maxVertex)];
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);

        polyscope::registerCurveNetwork("vertex isoline", iV, iE);

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
    }
    else{
        Eigen::MatrixXd iV;
        Eigen::MatrixXd iE;
        std::vector<int> f;
        std::tie(iV, iE, f) = getTimeFunctionIsoLine(globalGeometry, globalTimeFunction, isoVal);
        //polyscope::registerCurveNetwork("vertex isoline", iV, iE);
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
    }

    return std::make_pair(maxVertex, minVertex);
}
