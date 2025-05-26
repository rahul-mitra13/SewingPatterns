

#include "gmsh_helpers.h"

#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/remeshing.h"

using namespace std;

void parseMsh(
    const std::filesystem::path mshPath, 
    std::unique_ptr<ManifoldSurfaceMesh>& mesh, 
    std::unique_ptr<VertexPositionGeometry>& geometry, 
    globalBoundaryConditions& globalBdyConditions) {

    gmsh::initialize();
    gmsh::open(mshPath);
  
    // Get vertices, triangles and create GC mesh
    vector<size_t> nodeTags;
    vector<double> coord, paramCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord);
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags, elementNodeTagsByType;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTagsByType, 2);
    ensure(elementTypes.size() == 1 && elementTypes[0] == 2); // check that it's only triangles
  
    vector<Vector3> vertexPositions;
    map<size_t, size_t> nodeTagMap; // Gmsh tag to compact GC index
    for (int i = 0; i < nodeTags.size(); i++) {
      vertexPositions.push_back({coord[3*i], coord[3*i+1], coord[3*i+2]});
      nodeTagMap[nodeTags[i]] = i;
    }
  
    vector<vector<size_t>> elementNodeTags(elementTags[0].size(), vector<size_t>(3));
    for (int i = 0; i < elementNodeTags.size(); i++)
      for (int j = 0; j < 3; j++)
        elementNodeTags[i][j] = nodeTagMap[elementNodeTagsByType[0][3*i+j]];
  
    tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(elementNodeTags, vertexPositions);
  
    // //make the mesh Delaunay 
    // fixDelaunay(*mesh, *geometry); // we make the mesh approximately Delaunay

    // Define map from pair of vertices to edge
    map<pair<int,int>, int> vertexPairToHalfedge;
    for (Halfedge he : mesh->halfedges())
      vertexPairToHalfedge[{he.tailVertex().getIndex(), he.tipVertex().getIndex()}] = he.edge().getIndex();
      // vertexPairToHalfedge[{e.firstVertex().getIndex(), e.secondVertex().getIndex()}] = e.getIndex();
    
    // Get all physical group names
    vector<pair<int,int>> dimTags;
    gmsh::model::getPhysicalGroups(dimTags, 1);
    set<string> physicalGroupNames;
    for (auto &[dim,tag] : dimTags) {
      string name;
      gmsh::model::getPhysicalName(dim, tag, name);
      physicalGroupNames.insert(name);
    }

    // Get physical groups to find start and end loops
    vector<pair<int,int>> startEntities, endEntities;
    gmsh::model::getEntitiesForPhysicalName("start", startEntities);
    gmsh::model::getEntitiesForPhysicalName("end", endEntities);
    set<int> uniqueStartNodes, uniqueEndNodes;
    for (auto &[entDim, entTag] : startEntities) {
      vector<int> elementTypes;
      vector<vector<size_t>> elementTags, elementNodeTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
      uniqueStartNodes.insert(elementNodeTags[0].begin(), elementNodeTags[0].end());
      for (int i = 0; i < elementTags[0].size(); i++) {
        int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
        int edge = vertexPairToHalfedge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
        globalBdyConditions.courseBdyEdges.push_back(edge);
      }
    }
    for (auto &[entDim, entTag] : endEntities) {
      vector<int> elementTypes;
      vector<vector<size_t>> elementTags, elementNodeTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
      uniqueEndNodes.insert(elementNodeTags[0].begin(), elementNodeTags[0].end());
      for (int i = 0; i < elementTags[0].size(); i++) {
        int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
        int edge = vertexPairToHalfedge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
        globalBdyConditions.courseBdyEdges.push_back(edge);
      }
    }
    for (int node : uniqueStartNodes)
      globalBdyConditions.courseStartBoundaryVertices.push_back(nodeTagMap[node]);
    for (int node : uniqueEndNodes)
      globalBdyConditions.courseEndBoundaryVertices.push_back(nodeTagMap[node]);
  
    // Get physical groups for course alignment
    if (physicalGroupNames.count("alignCourse")) {
      vector<pair<int,int>> alignCourseEntities;
      gmsh::model::getEntitiesForPhysicalName("alignCourse", alignCourseEntities);
      for (auto &[entDim, entTag] : alignCourseEntities) {
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags, elementNodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
        for (int i = 0; i < elementTags[0].size(); i++) {
          int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
          int edge = vertexPairToHalfedge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
          globalBdyConditions.courseBdyEdges.push_back(edge);
        }
      }
    }

    // Get physical groups for wale alignment (TODO: what to do with edges?)
    if (physicalGroupNames.count("alignWale")) {
      vector<pair<int,int>> alignWaleEntities;
      gmsh::model::getEntitiesForPhysicalName("alignWale", alignWaleEntities);
      for (auto &[entDim, entTag] : alignWaleEntities) {
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags, elementNodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, entDim, entTag);
        for (int i = 0; i < elementTags[0].size(); i++) {
          int nodeTag1 = elementNodeTags[0][2*i], nodeTag2 = elementNodeTags[0][2*i+1];
          int edge = vertexPairToHalfedge[{nodeTagMap[nodeTag1], nodeTagMap[nodeTag2]}];
          // globalBdyConditions.courseBdyEdges.push_back(edge);
        }
      }
    }

    H(globalBdyConditions.courseStartBoundaryVertices.size());
    H(globalBdyConditions.courseEndBoundaryVertices.size());
    H(globalBdyConditions.courseBdyEdges.size());
  
    // vector<pair<int,int>> physicalGroups;
    // gmsh::model::getPhysicalGroups(physicalGroups);
    // for (auto &[dim,tag] : physicalGroups) {
    //   string name; gmsh::model::getPhysicalName(dim, tag, name);
    //   H(name);
    //   vector<int> entities; gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);
    //   H(entities);
    // }
    // H(physicalGroups);
  
    // gmsh::fltk::run();
  
  }