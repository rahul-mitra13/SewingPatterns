#include "path_constraints.h"

//rotate the face gradients clockwise 
FaceData<Vector3> clockWiseRotatedGradients(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    FaceData<Vector3> rotatedFaceGradients(globalMesh);
    globalGeometry.requireFaceNormals();
    //rotate the first gradient by 90 degrees wrt to the face normal
    for (Face f : globalMesh.faces()){
        rotatedFaceGradients[f] = globalFaceGradients[f].rotateAround(globalGeometry.faceNormals[f], -PI/2.);
        rotatedFaceGradients[f] = rotatedFaceGradients[f].normalize();
    }
    return rotatedFaceGradients;
}

//find maximum dot product halfedge
double maximumDotProduct(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& rotatedFaceGradients){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    double maxDotProd = -DBL_MAX;
    for (Halfedge he : globalMesh.halfedges()){
        Vector3 vector = (globalGeometry.vertexPositions[he.tipVertex()] - globalGeometry.vertexPositions[he.tailVertex()]).normalize();
        double dotProd = dot(vector, rotatedFaceGradients[he.face()]);
        if (dotProd > maxDotProd){
            maxDotProd = dotProd;
        }
    }
    return maxDotProd;
}

//construct glued halfedge weights 
HalfedgeData<double> constructGluedHalfedgeWeights(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, FaceData<Vector3>& rotatedFaceGradients,
                                            double maxDotProd){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;
    HalfedgeData<double> heWeights(gluedMesh, 0.);
    for (Face f : globalMesh.faces()){

        Halfedge hij = f.halfedge();
        double dotIJ = dot((globalGeometry.vertexPositions[hij.tipVertex()] - globalGeometry.vertexPositions[hij.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[gluedMesh.face(f.getIndex()).halfedge()] = std::fabs(dotIJ - maxDotProd);

        Halfedge hjk = hij.next();
        double dotJK = dot((globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[gluedMesh.face(f.getIndex()).halfedge().next()] = std::fabs(dotJK - maxDotProd);

        Halfedge hki = hjk.next();
        double dotKI = dot((globalGeometry.vertexPositions[hki.tipVertex()] - globalGeometry.vertexPositions[hki.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[gluedMesh.face(f.getIndex()).halfedge().next().next()] = std::fabs(dotKI - maxDotProd);
    }

    return heWeights;
}

//e1, e2 are in the global setting 
std::tuple<std::vector<double>, std::vector<double>> constructEdgePath(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Edge e1, Edge e2, 
                                        std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, FaceData<Vector3>& globalFaceGradients, HalfedgeData<double>& gluedHeWeights){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    std::vector<double> toReturnGlued(gluedMesh.nEdges(), 0.0);
    std::vector<double> toReturnGlobal(globalMesh.nEdges(), 0.0);

    //startVert, endVert are in the glued mesh setting 
    Vertex startVert, endVert;

    if (dot((globalGeometry.vertexPositions[e1.halfedge().tipVertex()] - globalGeometry.vertexPositions[e1.halfedge().tailVertex()]).normalize(), 
            globalFaceGradients[e1.halfedge().face()]) > dot((globalGeometry.vertexPositions[e1.halfedge().twin().tipVertex()] - globalGeometry.vertexPositions[e1.halfedge().twin().tailVertex()]).normalize(), 
            globalFaceGradients[e1.halfedge().twin().face()])){
        startVert = gluedMesh.vertex(vertexMap[e1.halfedge().tipVertex().getIndex()]);
    }
    else{
        startVert = gluedMesh.vertex(vertexMap[e1.halfedge().twin().tipVertex().getIndex()]);
    }

    if (dot((globalGeometry.vertexPositions[e2.halfedge().tipVertex()] - globalGeometry.vertexPositions[e2.halfedge().tailVertex()]).normalize(), 
            globalFaceGradients[e2.halfedge().face()]) > dot((globalGeometry.vertexPositions[e2.halfedge().twin().tipVertex()] - globalGeometry.vertexPositions[e2.halfedge().twin().tailVertex()]).normalize(),
            globalFaceGradients[e2.halfedge().twin().face()])){
        endVert = gluedMesh.vertex(vertexMap[e2.halfedge().tipVertex().getIndex()]);
    }
    else{
        endVert = gluedMesh.vertex(vertexMap[e2.halfedge().twin().tipVertex().getIndex()]);
    }

    // Early out for empty case
    if (startVert == endVert) {
        return std::tie(toReturnGlobal, toReturnGlued);
    }


    // Search state: incoming halfedges to each vertex, once discovered
    std::unordered_map<Vertex, Halfedge> incomingHalfedge;

    // Search state: visible neighbors eligible to expand to
    using WeightedHalfedge = std::tuple<double, Halfedge>;
    std::priority_queue<WeightedHalfedge, std::vector<WeightedHalfedge>, std::greater<WeightedHalfedge>> pq;

    // Helper to add a vertex's
    auto vertexDiscovered = [&](Vertex v) {
        return v == startVert || incomingHalfedge.find(v) != incomingHalfedge.end();
    };
    auto enqueueVertexNeighbors = [&](Vertex v, double dist) {
        for (Halfedge he : v.outgoingHalfedges()) {
            if (!vertexDiscovered(he.twin().vertex())) {
                double len = gluedHeWeights[he];
                double targetDist = dist + len;
                pq.emplace(targetDist, he);
            }
        }
    };

    // Add initial halfedges
    enqueueVertexNeighbors(startVert, 0.);
    std::vector<Halfedge> path;

    while (!pq.empty()) {

        // Get the next closest neighbor off the queue
        double currDist = std::get<0>(pq.top());
        Halfedge currIncomingHalfedge = std::get<1>(pq.top());
        pq.pop();

        Vertex currVert = currIncomingHalfedge.twin().vertex();
        if (vertexDiscovered(currVert)) continue;

        // Accept the neighbor
        incomingHalfedge[currVert] = currIncomingHalfedge;

        // Found path! Walk backwards to reconstruct it and return
        if (currVert == endVert) {
            Vertex walkV = currVert;
            while (walkV != startVert) {
                Halfedge prevHe = incomingHalfedge[walkV];
                path.push_back(prevHe);
                walkV = prevHe.vertex();
            }
            std::reverse(std::begin(path), std::end(path));
            for (Halfedge he : path){
            toReturnGlued[he.edge().getIndex()] = he.orientation() ? 1.0 : -1.0;
            }
            //toReturnGlobal = convertGluedToGlobalEdgeFunction(globalGeometry, gluedGeometry, toReturnGlued, edgeMap);
            for (Edge e : globalMesh.edges()){
                toReturnGlobal[e.getIndex()] = toReturnGlued[gluedMesh.edge(edgeMap[e.getIndex()]).getIndex()];
            }
            return std::tie(toReturnGlobal, toReturnGlued);
        }
        // Enqueue neighbors
        enqueueVertexNeighbors(currVert, currDist);
    }

    return std::tie(toReturnGlobal, toReturnGlued);
}
