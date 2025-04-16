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
    double penaltyFactor = 100.;//penalize anti aligned edges, ugh
    HalfedgeData<double> heWeights(gluedMesh, 0.);
    for (Face f : globalMesh.faces()){

        Halfedge hij = f.halfedge();
        double dotIJ = dot((globalGeometry.vertexPositions[hij.tipVertex()] - globalGeometry.vertexPositions[hij.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        //heWeights[hij] = std::fabs(dotIJ - maxDotProd);
        heWeights[hij] = dotIJ < 0. ? penaltyFactor * ((1 - dotIJ) / 2.) : (1 - dotIJ) / 2;
        //also set weights for boundary halfedges 
        if (!hij.twin().isInterior()){
            double dotIJTwin = dot((globalGeometry.vertexPositions[hij.tailVertex()] - globalGeometry.vertexPositions[hij.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            //heWeights[hij.twin()] = std::fabs(dotIJTwin - maxDotProd);
            heWeights[hij.twin()] = dotIJTwin < 0. ? penaltyFactor * (1 - dotIJTwin)/2. : (1 - dotIJTwin)/2.;
        }

        Halfedge hjk = hij.next();
        double dotJK = dot((globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        //heWeights[hjk] = std::fabs(dotJK - maxDotProd);
        heWeights[hjk] = dotJK < 0. ? penaltyFactor * (1 - dotJK) / 2. : (1 - dotJK) / 2;
        //also set weights for boundary halfedges 
        if (!hjk.twin().isInterior()){
            double dotJKTwin = dot((globalGeometry.vertexPositions[hjk.tailVertex()] - globalGeometry.vertexPositions[hjk.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            //heWeights[hjk.twin()] = std::fabs(dotJKTwin - maxDotProd);
            heWeights[hjk.twin()] = dotJKTwin < 0. ? penaltyFactor * (1 - dotJKTwin) / 2. : (1 - dotJKTwin) / 2.;

        }

        Halfedge hki = hjk.next();
        double dotKI = dot((globalGeometry.vertexPositions[hki.tipVertex()] - globalGeometry.vertexPositions[hki.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        //heWeights[hki] = std::fabs(dotKI - maxDotProd);
        heWeights[hki] = dotKI < 0. ? penaltyFactor * (1 - dotKI) / 2. : (1 - dotKI) / 2.;
        //also set weights for boundary halfedges 
        if (!hki.twin().isInterior()){
            double dotKITwin = dot((globalGeometry.vertexPositions[hki.tailVertex()] - globalGeometry.vertexPositions[hki.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            //heWeights[hki.twin()] = std::fabs(dotKITwin - maxDotProd);
            heWeights[hki.twin()] = dotKITwin < 0. ? penaltyFactor * (1 - dotKITwin) / 2. : (1 - dotKITwin) / 2.;
        }
    }

    return heWeights;
}

//construct an edge path between two vertices
//e1, e2 are in the global setting
std::tuple<std::vector<double>, std::vector<double>> constructEdgePath(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Edge e1, Edge e2,
                                        std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, FaceData<Vector3>& globalFaceGradients, HalfedgeData<double>& gluedHeWeights,
                                        HalfedgeData<double>& gluedSigmaTilde){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    std::vector<double> toReturnGlued(gluedMesh.nEdges(), 0.0);
    std::vector<double> toReturnGlobal(globalMesh.nEdges(), 0.0);

    //startVert, endVert are in the glued mesh setting
    Vertex startVert, endVert;
    Halfedge additionalStartHe, additionalEndHe;

    //don't take paths through singular edges
    // gluedHeWeights[e1.halfedge()] = DBL_MAX;
    // gluedHeWeights[e1.halfedge().twin()] = DBL_MAX;
    // gluedHeWeights[e2.halfedge()] = DBL_MAX;
    // gluedHeWeights[e2.halfedge().twin()] = DBL_MAX;

    if (dot((globalGeometry.vertexPositions[e1.halfedge().tipVertex()] - globalGeometry.vertexPositions[e1.halfedge().tailVertex()]).normalize(),
            globalFaceGradients[e1.halfedge().face()]) > dot((globalGeometry.vertexPositions[e1.halfedge().twin().tipVertex()] - globalGeometry.vertexPositions[e1.halfedge().twin().tailVertex()]).normalize(),
            globalFaceGradients[e1.halfedge().twin().face()])){

        if (gluedSigmaTilde[e1.halfedge()] < gluedSigmaTilde[e1.halfedge().twin()]){
            std::cout << "something went wrong " << std::endl;
            //exit(1);
        }
        startVert = gluedMesh.vertex(vertexMap[e1.halfedge().tipVertex().getIndex()]);

    }
    else{
        if (gluedSigmaTilde[e1.halfedge().twin()] < gluedSigmaTilde[e1.halfedge()]){
            std::cout << "something went wrong " << std::endl;
            //exit(1);
        }
        startVert = gluedMesh.vertex(vertexMap[e1.halfedge().twin().tipVertex().getIndex()]);
    }

    if (dot((globalGeometry.vertexPositions[e2.halfedge().tipVertex()] - globalGeometry.vertexPositions[e2.halfedge().tailVertex()]).normalize(),
            globalFaceGradients[e2.halfedge().face()]) > dot((globalGeometry.vertexPositions[e2.halfedge().twin().tipVertex()] - globalGeometry.vertexPositions[e2.halfedge().twin().tailVertex()]).normalize(),
            globalFaceGradients[e2.halfedge().twin().face()])){
        if (gluedSigmaTilde[e2.halfedge()] < gluedSigmaTilde[e2.halfedge().twin()]){
            std::cout << "something went wrong " << std::endl;
            //exit(1);
        }
        std::cout << "things are okay!! " << std::endl;
        endVert = gluedMesh.vertex(vertexMap[e2.halfedge().tipVertex().getIndex()]);
    }
    else{
        if (gluedSigmaTilde[e2.halfedge().twin()] < gluedSigmaTilde[e2.halfedge()]){
            std::cout << "something went wrong " << std::endl;
            //exit(1);
        }
        std::cout << "things are okay!! " << std::endl;
        endVert = gluedMesh.vertex(vertexMap[e2.halfedge().twin().tipVertex().getIndex()]);
    }

    // Early out for empty case
    if (startVert == endVert) {
        return std::tie(toReturnGlobal, toReturnGlued);
    }

    if (e1.halfedge().tipVertex() == startVert){
        additionalStartHe = e1.halfedge().next().twin();
    }
    else{
        additionalStartHe = e1.halfedge().twin().next().twin();
    }

    if (e2.halfedge().tipVertex() == endVert){
        additionalEndHe = e2.halfedge().twin().next().next().twin();
    }
    else{
        additionalEndHe = e2.halfedge().next().next().twin();
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

            //insert the additional halfedge at the start of the vector
            path.insert(path.begin(), additionalStartHe);
            //insert the additional halfedge at the end of the vector
            path.push_back(additionalEndHe);

            for (Halfedge he : path){
                toReturnGlued[he.edge().getIndex()] = he.orientation() ? 1.0 : -1.0;
            }
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

//update glued halfedge weights 
//set the halfedges that have been take by some edge path i.e., gluedPath to infinity 
void updateGluedHalfedgeWeights(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, std::vector<double>& gluedPath,
                                HalfedgeData<double>& gluedHeWeights){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    for (Edge e : gluedMesh.edges()){
        if (gluedPath[e.getIndex()] > 0){
            gluedHeWeights[e.halfedge().getIndex()] = DBL_MAX;
        }
        else if (gluedPath[e.getIndex()] < 0){
            gluedHeWeights[e.halfedge().twin().getIndex()] = DBL_MAX;
        }
    }

}