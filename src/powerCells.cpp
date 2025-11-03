#include "powerCells.h"
void computeCourseSingularities(powerCellOptions& options){

    EdgeLengthGeometry* gluedGeometry = options.gluedGeometry;
    SurfaceMesh& gluedMesh = gluedGeometry->mesh;
    VertexPositionGeometry* globalGeometry = options.globalGeometry;
    SurfaceMesh& globalMesh = globalGeometry->mesh;
    polyscope::SurfaceMesh* psMesh = options.psMesh;
    auto saddleLoops = options.saddleLoops;
    double period = options.period;
    // Create the Heat Method solver
    HeatMethodDistanceSolver heatSolver(*gluedGeometry);
    //first compute the measure in the course direction 
    VertexData<double> courseMeasure = computeCourseCurlMeasure(options); 

    //if we want to mask the saddle 
    if (options.maskSaddle){
        std::vector<Vertex> heatSourceVerts; 
        for (int i = 0; i < saddleLoops.size(); i++){
            std::vector<double> path = saddleLoops[i];
            for (int j = 0; j < path.size(); j++){
                if (std::fabs(path[j]) > 0){
                    Edge e = gluedGeometry->mesh.edge(j);
                    heatSourceVerts.push_back(e.halfedge().tailVertex());
                    heatSourceVerts.push_back(e.halfedge().tipVertex());
                }
            }
        }
        //mask the saddle 
        VertexData<double> allDist = heatSolver.computeDistance(heatSourceVerts);
        VertexData<double> courseWeighting(globalMesh);
        double maxVal = std::numeric_limits<double>::min();
        double maxSourceVal = std::numeric_limits<double>::min();
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
        //have a hard mask in the course direction
        for (Vertex v : globalMesh.vertices()){
            courseWeighting[v] = (allDist[v] > 1.5 * period);
            courseMeasure[v] = courseWeighting[v] * courseMeasure[v];
        }
    }//GENERAL TO-DO: \sigma should integrate to 0 around saddle loops

    //find the decompositions per face 
    FaceData<int> faceComponentsGlued =  componentsCutByLoops(*gluedGeometry, saddleLoops);

    // find k = #components
    int maxID = 0;
    for (Face f : gluedMesh.faces()) {
        maxID = std::max(maxID, faceComponentsGlued[f]);
    }
    int k = maxID + 1;
    //prepare buckets and assign vertex buckets
    std::vector<std::vector<Vertex>> buckets(k);
    //assign vertices (multi-membership)
    for (Vertex v : gluedMesh.vertices()) {
        std::vector<char> seen(k, 0);
        for (Face f : v.adjacentFaces()) {
            int cid = faceComponentsGlued[f];
            if (!seen[cid]) {
                buckets[cid].push_back(v);
                seen[cid] = 1;
            }
        }
    }
}

//compute the course curl measure in the glued setting 
VertexData<double> computeCourseCurlMeasure(powerCellOptions& options){

    SurfaceMesh& globalMesh = options.globalGeometry->mesh;
    SurfaceMesh& gluedMesh = options.gluedGeometry->mesh;
    VertexData<double> curlGlobal(globalMesh);
    VertexData<double> curlGlued(gluedMesh);
    options.globalGeometry->requireFaceAreas();
    auto vertexMap = options.vertexMap;
    auto field = options.normalizedTFGrad;
    for (Vertex vi : globalMesh.vertices()){
        double sum = 0.0;
        double area = 0.0; // area of the 1-ring of faces
        for (Halfedge he : options.gluedOneRingMap[vi.getIndex()]){
            Halfedge hjk = he.next();
            if (!hjk.isInterior()) continue;
            Vector3 hjkVec = options.globalGeometry->vertexPositions[hjk.tipVertex()] - options.globalGeometry->vertexPositions[hjk.tailVertex()];
            field[he.face()] = field[he.face()].normalize(); //always normalize the field
            area += options.globalGeometry->faceArea(he.face());
            sum += dot(hjkVec, field[he.face()]);
        }

        if (vi.isBoundary()) curlGlobal[vi] = 0.;//set curl to 0 on boundary vertices
        else curlGlobal[vi] = sum / area;
    }

    //convert curl to glued setting 
    for (Vertex v : globalMesh.vertices()){
        curlGlued[vertexMap[v.getIndex()]] = curlGlobal[v];
    }
   
    return curlGlued;
}