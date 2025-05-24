#include "helpers.h"
#include "path_constraints.h"

//rotate the face gradients clockwise 
FaceData<Vector3> clockWiseRotatedGradients(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& globalFaceGradients){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    FaceData<Vector3> rotatedFaceGradients(globalMesh);
    globalGeometry.requireFaceNormals();
    //rotate the first gradient by 90 degrees wrt to the face normal
    for (Face f : globalMesh.faces()){
        rotatedFaceGradients[f] = globalFaceGradients[f].rotateAround(globalGeometry.faceNormals[f], -PI/2.);
        //rotatedFaceGradients[f] = globalFaceGradients[f].rotateAround(globalGeometry.faceNormals[f], PI/2.);//debug
        rotatedFaceGradients[f] = rotatedFaceGradients[f].normalize();
    }
    return rotatedFaceGradients;
}

//find maximum dot product halfedge
double maximumDotProduct(VertexPositionGeometry& globalGeometry, FaceData<Vector3>& rotatedFaceGradients){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    double maxDotProd = -DBL_MAX;
    for (Halfedge he : globalMesh.halfedges()) if (he.isInterior()) {
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
    double penaltyFactor = DBL_MAX;//penalize anti aligned edges, ugh
    HalfedgeData<double> heWeights(gluedMesh, 0.);
    for (Face f : globalMesh.faces()){

        Face fGlued = gluedMesh.face(f.getIndex()); // corresponding face in glued mesh

        Halfedge hij = f.halfedge(), hijGlued = fGlued.halfedge();
        double dotIJ = dot((globalGeometry.vertexPositions[hij.tipVertex()] - globalGeometry.vertexPositions[hij.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[hijGlued] = dotIJ < 0. ? penaltyFactor * ((1 - dotIJ) / 2.) : (1 - dotIJ) / 2.;
        //also set weights for boundary halfedges 
        if (!hij.twin().isInterior()){
            double dotIJTwin = dot((globalGeometry.vertexPositions[hij.tailVertex()] - globalGeometry.vertexPositions[hij.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            heWeights[hijGlued.twin()] = dotIJTwin < 0. ? penaltyFactor * (1 - dotIJTwin)/2. : (1 - dotIJTwin)/2.;
        }

        Halfedge hjk = hij.next(), hjkGlued = hijGlued.next();
        double dotJK = dot((globalGeometry.vertexPositions[hjk.tipVertex()] - globalGeometry.vertexPositions[hjk.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[hjkGlued] = dotJK < 0. ? penaltyFactor * (1 - dotJK) / 2. : (1 - dotJK) / 2.;
        //also set weights for boundary halfedges 
        if (!hjk.twin().isInterior()){
            double dotJKTwin = dot((globalGeometry.vertexPositions[hjk.tailVertex()] - globalGeometry.vertexPositions[hjk.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            heWeights[hjkGlued.twin()] = dotJKTwin < 0. ? penaltyFactor * (1 - dotJKTwin) / 2. : (1 - dotJKTwin) / 2.;

        }
        Halfedge hki = hjk.next(), hkiGlued = hjkGlued.next();
        double dotKI = dot((globalGeometry.vertexPositions[hki.tipVertex()] - globalGeometry.vertexPositions[hki.tailVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
        heWeights[hkiGlued] = dotKI < 0. ? penaltyFactor * (1 - dotKI) / 2. : (1 - dotKI) / 2.;
        //also set weights for boundary halfedges 
        if (!hki.twin().isInterior()){
            double dotKITwin = dot((globalGeometry.vertexPositions[hki.tailVertex()] - globalGeometry.vertexPositions[hki.tipVertex()]).normalize(), 
                            rotatedFaceGradients[f].normalize());
            heWeights[hkiGlued.twin()] = dotKITwin < 0. ? penaltyFactor * (1 - dotKITwin) / 2. : (1 - dotKITwin) / 2.;
        }
    }
    return heWeights;
}

//construct an edge path between two vertices 
//e1, e2 are in the global setting 
std::tuple<std::vector<double>, std::vector<double>> constructEdgePath(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, Edge e1, Edge e2, 
                                        std::map<int, int>& vertexMap, std::map<int, int>& edgeMap, FaceData<Vector3>& globalFaceGradients, HalfedgeData<double>& gluedHeWeights,
                                        HalfedgeData<double>& gluedSigmaTilde, bool connectSaddles){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;
    SurfaceMesh& gluedMesh = gluedGeometry.mesh;

    std::vector<double> toReturnGlued(gluedMesh.nEdges(), 0.0);
    std::vector<double> toReturnGlobal(globalMesh.nEdges(), 0.0);

    //startVert, endVert are in the glued mesh setting 
    Vertex startVert, endVert;

    //don't take paths through singular edges 
    // gluedHeWeights[e1.halfedge()] = DBL_MAX;
    // gluedHeWeights[e1.halfedge().twin()] = DBL_MAX;
    // gluedHeWeights[e2.halfedge()] = DBL_MAX;
    // gluedHeWeights[e2.halfedge().twin()] = DBL_MAX;

    // Setup correct starting vertex: "top" of the edge.
    // TODO: What happens if e1 or e2 is a glued edge? Does it break down?
    Halfedge he1 = e1.halfedge(), he2 = e2.halfedge();
    Vector3 he1vec = (globalGeometry.vertexPositions[he1.tipVertex()] - globalGeometry.vertexPositions[he1.tailVertex()]).normalize();
    Vector3 he2vec = (globalGeometry.vertexPositions[he2.tipVertex()] - globalGeometry.vertexPositions[he2.tailVertex()]).normalize();
    Vector3 he1twinvec = (globalGeometry.vertexPositions[he1.twin().tipVertex()] - globalGeometry.vertexPositions[he1.twin().tailVertex()]).normalize();
    Vector3 he2twinvec = (globalGeometry.vertexPositions[he2.twin().tipVertex()] - globalGeometry.vertexPositions[he2.twin().tailVertex()]).normalize();
    Vector3 he1grad = globalFaceGradients[he1.face()];
    Vector3 he2grad = globalFaceGradients[he2.face()];
    Vector3 he1twingrad = globalFaceGradients[he1.twin().face()];
    Vector3 he2twingrad = globalFaceGradients[he2.twin().face()];
    if (dot(he1vec, he1grad) > dot(he1twinvec, he1twingrad)) {
        ensure(gluedSigmaTilde[he1] > gluedSigmaTilde[he1.twin()]);
        startVert = gluedMesh.vertex(vertexMap[he1.tipVertex().getIndex()]);
    } else{
        ensure(gluedSigmaTilde[he1.twin()] > gluedSigmaTilde[he1]);
        startVert = gluedMesh.vertex(vertexMap[he1.twin().tipVertex().getIndex()]);
    }
    if (dot(he2vec, he2grad) > dot(he2twinvec, he2twingrad)) {
        ensure(gluedSigmaTilde[he2] > gluedSigmaTilde[he2.twin()]);
        endVert = gluedMesh.vertex(vertexMap[he2.tipVertex().getIndex()]);
    } else{
        ensure(gluedSigmaTilde[he2.twin()] > gluedSigmaTilde[he2]);
        endVert = gluedMesh.vertex(vertexMap[he2.twin().tipVertex().getIndex()]);
    }

    // Early out for empty case
    if (startVert == endVert) {
        return std::tie(toReturnGlobal, toReturnGlued);
    }

    // Setup half-edges that connect to the saddle vertex
    Halfedge additionalStartHe, additionalEndHe;
    if (he1.tipVertex() == startVert)
        additionalStartHe = e1.halfedge().next().twin();
    else
        additionalStartHe = e1.halfedge().twin().next().twin();

    if (e2.halfedge().tipVertex() == endVert)
        additionalEndHe = e2.halfedge().twin().next().next().twin();
    else
        additionalEndHe = e2.halfedge().next().next().twin();

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

            // Insert additional half-edges at start and end to connect sadles
            if (connectSaddles) {
                path.insert(path.begin(), additionalStartHe);
                path.push_back(additionalEndHe);
            }
                
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

//compute the path cost between two singular vertices 
double computePathCost(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& heWeights, Vertex& startVert, Vertex& endVert){
    
    SurfaceMesh& globalMesh = globalGeometry.mesh;

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
                double len = heWeights[he];
                double targetDist = dist + len;
                pq.emplace(targetDist, he);
            }
        }
    };

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
        }
        // Enqueue neighbors
        enqueueVertexNeighbors(currVert, currDist);
    }

    double pathCost = 0; 

    for (Halfedge he : path){
        pathCost += heWeights[he];
    }

    return pathCost;

}

//given a set of singularities perform an optimal matching between them 
//always done in the global setting for now 
std::vector<std::pair<Vertex, Vertex>> performOptimalMatching(VertexPositionGeometry& globalGeometry, HalfedgeData<double>& heWeights, std::vector<std::pair<Vertex, int>>& singularities){

    SurfaceMesh& globalMesh = globalGeometry.mesh;
    
    int numSingVertices = singularities.size();
    std::vector<std::pair<Vertex, Vertex>> toReturn;

    //compute the cost matrix 
    Eigen::MatrixXd C(numSingVertices, numSingVertices);
    C.setZero();

    //construct the cost matrix 
    for (int i = 0; i < numSingVertices; i++){
        int singValI = singularities[i].second;
        for (int j = 0; j < numSingVertices; j++){
            int singValJ = singularities[j].second;
            if (singValI < 0 && singValJ > 0){
                C(i, j) = 10000;//put high costs on paths from negative singularities to positive singularities
            }
            else if (singValI == singValJ) C(i, j) = 10000; //put large weights on singularities of the same sign
            else{
                C(i, j) = computePathCost(globalGeometry, heWeights, singularities[i].first, singularities[j].first);
            }
        }
    }

    //solve the optimization model
    using namespace std;
    try {
        
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "1-form computation.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        //add a matrix of Gurobi variables that signify the matching
        std::vector<std::vector<GRBVar>> T;
    
        for (int i = 0; i < numSingVertices; i++){
            std::vector<GRBVar> currRow;
            for (int j = 0; j < numSingVertices; j++){
                GRBVar t = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
                currRow.push_back(t);
            }
            T.push_back(currRow);
        }

        //first constraint (T\mathbb{1} = \mathbb{1})
        for (int i = 0; i < numSingVertices; i++){
            GRBLinExpr lhs = 0;
            for (int j = 0; j < numSingVertices; j++){
                lhs += T[i][j];
            }
            model.addConstr(lhs == 1.0);
        }

        //second constraint (\mathbb{1}^T T = \mathbb{1})
        for (int i = 0; i < numSingVertices; i++){
            GRBLinExpr lhs = 0;
            for (int j = 0; j < numSingVertices; j++){
                lhs += T[j][i];
            }
            model.addConstr(lhs == 1.0);
        }

        //add the objective
        GRBLinExpr matProd = 0;
        for (int i = 0; i < numSingVertices; i++){
            for (int j = 0; j < numSingVertices; j++){
                matProd += C(i, j) * T[i][j];
            }
        }

        GRBLinExpr obj = matProd;

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();
        
        for (int i = 0; i < numSingVertices; i++){
            for (int j = 0; j < numSingVertices; j++){
                if (T[i][j].get(GRB_DoubleAttr_X) > 0 && (singularities[i].second > 0 && singularities[j].second < 0)){
                    std::cout << "match vertex " << singularities[i].first << " to vertex " << singularities[j].first << std::endl;
                    toReturn.push_back(std::make_pair(singularities[i].first, singularities[j].first));//always return pairs in (posVertex, negVertex) order
                }
            }
        }

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return toReturn;

}


//----------------New path constraint mechanisms----------------//


//trace an isoline between to vertices 
std::vector<double> findIsoPath(SurfaceMesh& mesh,
                                  VertexData<double>& timeFunction,
                                  Vertex vStart,
                                  Vertex vEnd) {

    // Target function value to stay near
    double target = 0.5 * (timeFunction[vStart] + timeFunction[vEnd]);

    // Priority queue: (cost, vertex)
    using QueueEntry = std::pair<double, Vertex>;
    auto cmp = [](const QueueEntry& a, const QueueEntry& b) { return a.first > b.first; };
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, decltype(cmp)> queue(cmp);

    // Bookkeeping
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());
    VertexData<Halfedge> prevHalfedge(mesh); // for path reconstruction

    dist[vStart] = 0.0;
    queue.push({0.0, vStart});

    std::unordered_set<Vertex> visited;

    while (!queue.empty()) {
        auto [curCost, v] = queue.top();
        queue.pop();

        if (visited.count(v)) continue;
        visited.insert(v);

        if (v == vEnd) break;

        for (Halfedge he : v.outgoingHalfedges()) {
            if (he.isInterior()) {
                Vertex neigh = he.tipVertex();

                double cost = std::abs(timeFunction[neigh] - target);
                double newDist = dist[v] + cost;

                if (newDist < dist[neigh]) {
                    dist[neigh] = newDist;
                    prevHalfedge[neigh] = he;
                    queue.push({newDist, neigh});
                }
            }
        }
    }

    // Reconstruct path
    std::vector<Halfedge> path;
    Vertex v = vEnd;

    while (v != vStart) {
        if ((prevHalfedge[v] == Halfedge())) {
            throw std::runtime_error("No path found between vStart and vEnd");
        }
        Halfedge he = prevHalfedge[v];
        path.push_back(he);
        v = he.tailVertex();
    }

    std::reverse(path.begin(), path.end());

    //convert path to edge scalar
    std::vector<double> edgePath(mesh.nEdges(), 0.0);
    for (Halfedge he : path){
        if (he.orientation())
            edgePath[he.edge().getIndex()] = 1.0;
        else 
            edgePath[he.edge().getIndex()] = -1.0;
    }

    return edgePath;
}


std::vector<std::tuple<Vector3, Halfedge>>
computeIsolineFaceIntersection(const Face face, double iso,
                               const VertexData<double>& fValues,
                               const VertexPositionGeometry& geom) {
    std::vector<std::tuple<Vector3, Halfedge>> intersections;

    for (Halfedge he : face.adjacentHalfedges()){
        Vertex v0 = he.vertex();
        Vertex v1 = he.next().vertex();

        double f0 = fValues[v0];
        double f1 = fValues[v1];

        if ((f0 - iso) * (f1 - iso) < 0.0 || f0 == iso || f1 == iso) { // scalar crosses iso
            double t = (iso - f0) / (f1 - f0);
            Vector3 p0 = geom.inputVertexPositions[v0];
            Vector3 p1 = geom.inputVertexPositions[v1];
            Vector3 p = (1.0 - t) * p0 + t * p1;
            intersections.emplace_back(p, he);
        }
    }
    return intersections;
}


//trace the faces that an isoline passes through
//first argument is the vector of faces 
//second is argument is currently being used as an error code 
std::tuple<std::vector<Face>, int> traceIsolineFaces(
    const VertexPositionGeometry& globalGeometry,
    const VertexData<double>& timeFunction,
    const FaceData<Vector3>& timeFunctionGradients,
    const FaceData<Vector3>& rotatedFaceGradients,
    double isoVal,
    Halfedge startHe,
    Halfedge endHe
){

    std::unordered_set<Face> visitedFaces;
    std::vector<Face> pathFaces;
    
    Face currentFace = startHe.face();//face to start at 
    Face endFace = endHe.face();//face to end at

    //expand face
    std::queue<Face> frontier;
    frontier.push(currentFace);
    visitedFaces.insert(currentFace);
    visitedFaces.insert(endFace);

    while(!frontier.empty()){

        Face f = frontier.front();
        frontier.pop();
        if (f == endFace) break;

        auto intersections = computeIsolineFaceIntersection(f, isoVal, timeFunction, globalGeometry);
        if (intersections.size() != 2){
            std::cout << "intersection size = " << intersections.size() << std::endl;
            break;
        }
        auto [p1, he1] = intersections[0];
        auto [p2, he2] = intersections[1];
        Vector3 v1 = (p2 - p1).normalize();
        Vector3 v2 = (p1 - p2).normalize();
        Face neighbor;
        if (dot(v1, rotatedFaceGradients[f]) > 0){
            neighbor = he2.twin().face();
        }
        else{
            neighbor = he1.twin().face();
        }
        if (neighbor == endFace){
            return std::make_pair(pathFaces, 0);
        }
        if (visitedFaces.find(neighbor) == visitedFaces.end()){
            pathFaces.emplace_back(neighbor);
            frontier.push(neighbor);
            visitedFaces.insert(neighbor);
        }
        
    }

    if (isoVal < std::min(timeFunction[startHe.tipVertex()], timeFunction[startHe.tailVertex()]) || isoVal > std::max(timeFunction[startHe.tipVertex()], timeFunction[startHe.tailVertex()])){
        std::cout << "average doesn't lie on start halfedge " << std::endl;
    }
    if (isoVal < std::min(timeFunction[endHe.tipVertex()], timeFunction[endHe.tailVertex()]) || isoVal > std::max(timeFunction[endHe.tipVertex()], timeFunction[endHe.tailVertex()])){
        std::cout << "average doesn't lie on end halfedge " << std::endl;
    }

    return std::make_pair(pathFaces, -1);

}