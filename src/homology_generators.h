#include "geometrycentral/surface/vertex_position_geometry.h"
#include <vector>

/*
 * Determines whether or not a halfedge is in the primal spanning tree
 */
bool inPrimalSpanningTree(const Halfedge &he, const std::unordered_map<Vertex, Vertex> &tree) {
    const Vertex v = he.vertex();
    const Vertex w = he.next().vertex();
    if (tree.count(v) == 0 || tree.count(w) == 0) {
        return false;
    }
    return tree.at(v) == w || tree.at(w) == v;
}


/*
 * Given a co-tree structure, determines whether a halfedge is in
 * the dual spanning tree.
 */
bool inDualSpanningTree(const Halfedge &he, const std::unordered_map<Face, Face> &cotree) {
    const Face F = he.face();
    // the adjacent face to this halfedge
    const Face adjF = he.twin().face();
    // returns true if either face is a parent of the other
    if (cotree.count(F) == 0 || cotree.count(adjF) == 0) {
        return false;
    }
    return cotree.at(F) == adjF || cotree.at(adjF) == F;
}


Halfedge sharedHalfedge(const Face &f, const Face &g) {
    for (const Halfedge &he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // should never reach this if f and g are adjacent
    return f.halfedge();
}

Halfedge sharedHalfedge(const Vertex &v, const Vertex &w) {
    for (const Halfedge &he : v.outgoingHalfedges()) {
        if (he.tipVertex() == w) {
            return he;
        }
    }
    // should never reach this if f and g are adjacent
    return v.halfedge();
}


/*
 * Constructs a dictionary that maps each vertex of a mesh to it parent.
 * This dictionary represents the tree structure of a primal spanning tree.
 */
std::unordered_map<Vertex, Vertex> buildPrimalSpanningTree(ManifoldSurfaceMesh &mesh,
        const std::unordered_map<Face, Face> &cotree ) {

    std::unordered_map<Vertex, Vertex> vertexParent;

    // choose a random vertex to be the root
    const Vertex rootVertex = mesh.vertex(0);
    std::queue<Vertex> bag;
    bag.push(rootVertex);
    vertexParent[rootVertex] = rootVertex;

    // breadth first search
    while (bag.size() != 0) {
        const Vertex v = bag.front();
        for (const Halfedge &adjHe : v.outgoingHalfedges()) {
            if (vertexParent.count(adjHe.tipVertex()) == 0 && !inDualSpanningTree(adjHe, cotree)) {
                vertexParent[adjHe.tipVertex()] = v;
                bag.push(adjHe.tipVertex());
            }
        }
        bag.pop();
    }
    return vertexParent;
}

/*
 * Constructs a dictionary that maps each face of a mesh to it parent.
 * This dictionary represents the tree structure of a dual spanning tree.
 */
std::unordered_map<Face, Face> buildDualSpanningTree(ManifoldSurfaceMesh &mesh) {

    std::unordered_map<Face, Face> cotree;

    // if the mesh has boundary, then make one of the boundaries a root face
    const Face rootFace = mesh.hasBoundary() ? mesh.boundaryLoop(0).asFace() :
                                               mesh.face(0);
    std::queue<Face> bag;
    bag.push(rootFace);
    cotree[rootFace] = rootFace;

    // first iteration of breadth first search
    // this step adds a single dual halfedge to the spanning tree 
    const Face boundary = bag.front();
    bag.pop();
    for (const Halfedge &he : boundary.adjacentHalfedges()) {
        cotree[he.twin().face()] = boundary;
        bag.push(he.twin().face());
        break;
    }
    
    // breadth first search
    while (bag.size() != 0) {
        const Face f = bag.front();
        for (const Halfedge &he : f.adjacentHalfedges()) {
            if ((cotree.count(he.twin().face()) == 0) && !he.edge().isBoundary()) {
                bag.push(he.twin().face());
                cotree[he.twin().face()] = f;
            }
        }
        bag.pop();
    }
    return cotree;
}

/*
 * When i first add a dual halfedge connected to the boundary, only want to 
 * add one. and then never added any dualedges that cross boundaries.
 *
 */

/*
 * Builds the homology generators of a mesh.
 */
std::vector<std::vector<Halfedge>> buildHomologyGenerators(ManifoldSurfaceMesh &mesh) {

    const auto cotree = buildDualSpanningTree(mesh);
    std::cout << "built dual tree" << std::endl;
    const auto tree = buildPrimalSpanningTree(mesh, cotree);
    std::cout << "built primal tree" << std::endl;

    std::vector<Halfedge> edgeGenerators;
    for (const Edge &e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), tree) and !inDualSpanningTree(e.halfedge(), cotree)) {
            edgeGenerators.push_back(e.halfedge());
        }
    }
    std::cout << "found edgeGenerators" << std::endl;
    std::cout << "found " << edgeGenerators.size() << "generators" << std::endl;

    std::vector<std::vector<Halfedge>> homologyGenerators;
    // find all of the loops corresponding to the edgeGenerators
    for (const Halfedge &he : edgeGenerators) {
        Vertex currentVertex = he.tipVertex();
        // construct the path to the root face
        std::vector<Halfedge> pathToRoot1;
        do {
            const Vertex parentVertex = tree.at(currentVertex);
            const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
            pathToRoot1.push_back(sharedHe);
            currentVertex = parentVertex;
        }
        while ( tree.at(currentVertex) != currentVertex);

        currentVertex = he.tailVertex();
        // construct the path to the root face
        std::vector<Halfedge> pathToRoot2;
        do {
            const Vertex parentVertex = tree.at(currentVertex);
            const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
            pathToRoot2.push_back(sharedHe);
            currentVertex = parentVertex;
        }
        while (tree.at(currentVertex) != currentVertex);

        std::vector<Halfedge> homologyRing;
        homologyRing.push_back(he);
        // only add the halfedges that are not in both paths to the homology ring
        for (const Halfedge &he : pathToRoot1) {
            bool inOtherPath = false;
            for (const Halfedge &pathHe : pathToRoot2) {
                if (he == pathHe || he.twin() == pathHe) {
                    inOtherPath = true;
                    break;
                }
            }
            if (!inOtherPath) {
                homologyRing.push_back(he);
            }
        }
        // only add the halfedge that are not in both paths to the homology ring
        for (const Halfedge &he : pathToRoot2) {
            bool inOtherPath = false;
            for (const Halfedge &pathHe : pathToRoot1) {
                if (he == pathHe || he.twin() == pathHe) {
                    inOtherPath = true;
                    break;
                }
            }
            if (!inOtherPath) {
                homologyRing.push_back(he);
            }
        }
        homologyGenerators.push_back(homologyRing);
    }
    return homologyGenerators;
}

void visualizeHomologyGenerators(const std::vector<std::vector<Halfedge>> &homologyGenerators, VertexPositionGeometry &geometry) {
    int numRings = 0;
    for (const auto &homologyRing : homologyGenerators) {
        std::string name = "Homology Generator" + std::to_string(numRings++);
        std::vector<Vector3> positions;
        std::vector<std::array<int, 2>> edgeIndices;
        int nodeCounter = 0;
        for (const Halfedge &he : homologyRing) {
            const auto p1 = geometry.vertexPositions[he.tailVertex()];
            const auto p2 = geometry.vertexPositions[he.tipVertex()];
            positions.push_back(p1);
            positions.push_back(p2);
            edgeIndices.push_back({nodeCounter, nodeCounter+1});
            nodeCounter += 2;
        }
        polyscope::registerCurveNetwork(name, positions, edgeIndices);
    }
}

/*
 * This is the only function you need to worry about.
 *
 */
// std::vector<std::vector<double>> findAllHomologyGenerators(VertexPositionGeometry &geometry,
//         ManifoldSurfaceMesh &mesh) {
    
//     const auto allGenerators = buildHomologyGenerators(mesh);
//     visualizeHomologyGenerators(allGenerators, geometry);
//     std::vector<std::vector<double>> allGeneratorsVector;
//     int i = 0;
//     for (const auto &generator : allGenerators) {
//         std::string ringname = "homology ring " + std::to_string(i++);
//         std::vector<double> ring(mesh.nEdges());
//         for (const Halfedge &he : generator) {
//             if (he.orientation()) {
//                 ring[he.edge().getIndex()] = 1.0;
//             }
//             else {
//                 ring[he.edge().getIndex()] = -1.0;
//             }
//         }
//         allGeneratorsVector.push_back(ring);
//         psMesh->addEdgeScalarQuantity(ringname, ring);
//     }
//     return allGeneratorsVector;

// }

std::vector<Edge> halfedgesToEdges(const std::vector<Halfedge> &halfedges) {
    std::vector<Edge> edges;
    for (const Halfedge &he: halfedges) {
        edges.push_back(he.edge());
    }
    return edges;
}


/*
 * Builds the homology generators of a mesh.
 */
std::vector<std::vector<double>> buildHomologyGeneratorsVector(VertexPositionGeometry &geometry, ManifoldSurfaceMesh& mesh) {
    
    const auto cotree = buildDualSpanningTree(mesh);
    const auto tree = buildPrimalSpanningTree(mesh, cotree);

    std::vector<Halfedge> edgeGenerators;
    for (const Edge &e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), tree) and !inDualSpanningTree(e.halfedge(), cotree)) {
            edgeGenerators.push_back(e.halfedge());
        }
    }

    std::vector<std::vector<Halfedge>> homologyGeneratorsHalfedges;
    std::vector<std::vector<double>> homologyGenerators;
    // find all of the loops corresponding to the edgeGenerators
    int i = 0;
    for (const Halfedge &he : edgeGenerators) {
        std::string ringname = "homology ring " + std::to_string(i++);
        Vertex currentVertex = he.tipVertex();
        // construct the path to the root face
        std::vector<double> pathToRoot1(mesh.nEdges());
        std::vector<Halfedge> pathToRoot1Halfedges;
        do {
            const Vertex parentVertex = tree.at(currentVertex);
            const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
            if (sharedHe.orientation()) {
                pathToRoot1[sharedHe.edge().getIndex()] = 1.0;
            }
            else {
                pathToRoot1[sharedHe.edge().getIndex()] = -1.0;
            }
            pathToRoot1Halfedges.push_back(sharedHe);
            currentVertex = parentVertex;
        }
        while ( tree.at(currentVertex) != currentVertex);

        currentVertex = he.tailVertex();
        // construct the path to the root face
        std::vector<double> pathToRoot2(mesh.nEdges());
        std::vector<Halfedge> pathToRoot2Halfedges;
        do {
            const Vertex parentVertex = tree.at(currentVertex);
            const Halfedge sharedHe = sharedHalfedge(currentVertex, parentVertex);
            if (sharedHe.orientation()) {
                pathToRoot2[sharedHe.edge().getIndex()] = -1.0;
            }
            else {
                pathToRoot2[sharedHe.edge().getIndex()] = 1.0;
            }
            pathToRoot2Halfedges.push_back(sharedHe);
            currentVertex = parentVertex;
        }
        while (tree.at(currentVertex) != currentVertex);

        std::vector<double> homologyRing(mesh.nEdges());
        if (he.orientation()) {
            homologyRing[he.edge().getIndex()] = 1.0;
        }
        else {
            homologyRing[he.edge().getIndex()] = -1.0;
        }
        // only add the halfedges that are not in both paths to the homology ring
        for (int i = 0; i < pathToRoot1.size(); i++) {
            if (pathToRoot1[i] != 0 && pathToRoot2[i] == 0) {
                homologyRing[i] = pathToRoot1[i];
            }
        }
        for (int i = 0; i < pathToRoot2.size(); i++) {
            if (pathToRoot2[i] != 0 && pathToRoot1[i] == 0) {
                homologyRing[i] = pathToRoot2[i];
            }
        }


        std::vector<Halfedge> homologyRingHalfedges;
        homologyRingHalfedges.push_back(he);
        // only add the halfedges that are not in both paths to the homology ring
        for (const Halfedge &he : pathToRoot1Halfedges) {
            bool inOtherPath = false;
            for (const Halfedge &pathHe : pathToRoot2Halfedges) {
                if (he == pathHe || he.twin() == pathHe) {
                    inOtherPath = true;
                    break;
                }
            }
            if (!inOtherPath) {
                homologyRingHalfedges.push_back(he);
            }
        }
        // only add the halfedge that are not in both paths to the homology ring
        for (const Halfedge &he : pathToRoot2Halfedges) {
            bool inOtherPath = false;
            for (const Halfedge &pathHe : pathToRoot1Halfedges) {
                if (he == pathHe || he.twin() == pathHe) {
                    inOtherPath = true;
                    break;
                }
            }
            if (!inOtherPath) {
                homologyRingHalfedges.push_back(he);
            }
        }

        homologyGenerators.push_back(homologyRing);
        homologyGeneratorsHalfedges.push_back(homologyRingHalfedges);
        //psMesh->addEdgeScalarQuantity(ringname, homologyRing);
    }
    visualizeHomologyGenerators(homologyGeneratorsHalfedges, geometry);
    return homologyGenerators;
}
