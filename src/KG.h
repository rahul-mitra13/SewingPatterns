//geometry-central includes 
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/surface_point.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

//file includes 
#include "helpers.h"

#pragma once 
#include <vector> 
#include <optional>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <unordered_map>
#include <cmath>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//redefining for consistency
struct KGVertex{
    int id = -1;
    Vector3 position;//position embedded in R^3 (if position information is available)
    Vector3 baryCoords;//barycentric coordinates of the vertex if position information is not available
    // Row neighbors
    KGVertex* row_in_vertex  = nullptr;  
    KGVertex* row_out_vertex = nullptr; 
    // Column neighbors 
    std::array<KGVertex*, 2> col_in_vertex  { nullptr, nullptr };
    std::array<KGVertex*, 2> col_out_vertex { nullptr, nullptr };

    double alpha_tag = -1;
    double beta_tag = -1;
    std::optional<Halfedge> halfedge;//associated halfedge of a vertex (for virtual vertices)
    SurfacePoint surfacePoint;//surfacePoint representation of this knit graph vertex
    bool isAlphaVirtual = false;//is a virtual vertex in the course direction
    bool isBetaVirtual = false;//is a virtual vertex in the wale direction 

    void printVertexInfo(){
          
        std::cout << "vertex id = " << id << std::endl;
        if (row_in_vertex != nullptr) std::cout << "row in = " << row_in_vertex->id << std::endl;
        else std::cout << "row_in = -1" << std::endl;
        if (row_out_vertex != nullptr) std::cout << "row out = " << row_out_vertex->id << std::endl;
        else std::cout << "row_out = -1" << std::endl;
        if (col_in_vertex[0] != nullptr) std::cout << "col_in[0] = " << col_in_vertex[0]->id << std::endl;
        else std::cout << "col_in[0] = -1" << std::endl;
        if (col_in_vertex[1] != nullptr) std::cout << "col_in[1] = " << col_in_vertex[1]->id << std::endl;
        else std::cout << "col_in[1] = -1" << std::endl;
        if (col_out_vertex[0] != nullptr) std::cout << "col_out[0] = " << col_out_vertex[0]->id << std::endl;
        else std::cout << "col_out[0] = -1" << std::endl;
        if (col_out_vertex[1] != nullptr) std::cout << "col_out[1] = " << col_out_vertex[1]->id << std::endl;
        else std::cout << "col_out[1] = -1" << std::endl;
        std::cout << "isAlphaVirtual = " << isAlphaVirtual << std::endl;
        std::cout << "isBetaVirtual = " << isBetaVirtual << std::endl;
    }
};

class KG{

    private:
        
        //id counter of vertices 
        int vertexID = 0;

        //geometry on which this graph lives 
        VertexPositionGeometry *globalGeometry;

        //edge length geometry if we want to do this in the embedded setting 
        EdgeLengthGeometry *gluedGeometry;

        //polyscope object for this geometry 
        polyscope::SurfaceMesh *psMesh; 

        //period for sampling the course stripes 
        double coursePeriod; 

        //period for sampling the wale stripes 
        double walePeriod;

        // Flag of edges that were glued together
        EdgeData<bool> isGlued;

        //course stripe 1-form 
        CornerData<double> courseOneForm; 

        //edge indices in the course direction in the global setting
        EdgeData<double> courseSingularEdges; 

        //edge indices in the course direction in the glued setting
        EdgeData<double> courseSingularEdgesGlued;

        //wale stripe 1-form
        CornerData<double> waleOneForm;

        //edge indices in the wale direction in the global setting
        EdgeData<double> waleSingularEdges;

        //edge indices in the wale direction in the glued setting
        EdgeData<double> waleSingularEdgesGlued;

        //edge map from global mesh to glued mesh 
        std::map<int, int> *globalToGluedEdgeMap;

        //vertices in the graph (including virtuals)
        std::vector<std::unique_ptr<KGVertex>> allVertices;

        //vertices in the graph after adjusting the level sets
        std::vector<std::unique_ptr<KGVertex>> adjustedVertices;

        // Pairs of stitched (real) vertices
        std::vector<std::pair<int, int>> stitchedVertices;

        //course level sets per face (will become useful when we update the vertex positions of virtual vertices)
        FaceData<std::vector<double>> faceCourseLevelSets; 

        //wale level sets per face (will become useful when we update the vertex positions of virtual vertices)
        FaceData<std::vector<double>> faceWaleLevelSets; 

        //graph vertices per face before adjustment
        FaceData<std::vector<KGVertex*>> faceKGVertices;

        //graph vertices per face after adjustment
        FaceData<std::vector<KGVertex*>> adjustedFaceKGVertices;

        //vertices matched on singular edges in the course direction 
        std::vector<std::pair<KGVertex*, KGVertex*>> courseMatchings;

        //vertices matched on singular edges in the wale direction 
        std::vector<std::pair<KGVertex*, KGVertex*>> waleMatchings;


        private: 

            //hashing floating point numbers 
            int hashFloat(double val);

            //make virtual vertices in the course direction
            void makeCourseVirtualVertices();

            //make virtual vertices in the wale direction 
            void makeWaleVirtualVertices();

            //make virtual vertices on the border of all faces
            void makeVirtualVerticesOnBorder(Face& f, bool isCourseDirection);

            //make the real vertices 
            void makeRealVertices();

            //compute real vertices
            void makeRealVerticesOnInterior(Face& f);

            //make connections over a face 
            void makeFaceConnections();

            //render a set of knit graph nodes as a point cloud 
            void renderKnitGraphPoints();

            //render a set of knit graph nodes a curve network 
            void renderKnitGraphCurveNetwork();

            //intrinsic merging of virtual vertices on the triangle borders 
            void intrinsicMerge();

            //get the position of a KG vertex embedded in R^3
            Vector3 getKGPosition(const KGVertex *v){
                int fIndex = v->halfedge->face().getIndex();
                Face f = globalGeometry->mesh.face(fIndex);
                Vertex vI = f.halfedge().vertex();
                Vertex vJ = f.halfedge().next().vertex();
                Vertex vK = f.halfedge().next().next().vertex();
                Vector3 pI = globalGeometry->vertexPositions[vI];
                Vector3 pJ = globalGeometry->vertexPositions[vJ];
                Vector3 pK = globalGeometry->vertexPositions[vK];
                Vector3 pos = v->baryCoords[0] * pI + v->baryCoords[1] * pJ + v->baryCoords[2] * pK;
                return pos;
            }
            

            //now that we have the matchings, we update their barycoords and alpha/beta tags
            void updateSingularMatchings();

            //after updating the barycoords for matched vertices, update the stripe level sets over the faces 
            void adjustFaceLevelSets();

            //make virtual vertices on the triangle border in the course direction with the adjusted level sets 
            void makeAdjustedCourseVirtualVertices();

            //make virtual vertices on the triangle border in the wale direction with the adjusted level sets 
            void makeAdjustedWaleVirtualVertices();
            
            //make virtual vertices on the border of all faces with the adjusted level sets
            void makeAdjustedVirtualVerticesOnBorder(Face& f, bool isCourseDirection);

            //make adjusted real vertices
            void makeAdjustedRealVertices();

            //compute real vertices on the triangle interior with adjusted level sets
            void makeAdjustedRealVerticesOnInterior(Face& f);

            
        public: 

            //Default Constructor 
            KG(){};

            //Constructor
            KG(VertexPositionGeometry& globalGeometry,
               EdgeLengthGeometry& gluedGeometry,
               polyscope::SurfaceMesh& psMesh,
               double coursePeriod, double walePeriod,
               CornerData<double>& courseOneForm,
               EdgeData<double>& courseSingularEdges,
               CornerData<double>& waleOneForm,
               EdgeData<double>& waleSingularEdges,
               std::map<int, int>& globalToGluedEdgeMap);

            //get the vertices in this knit graph (including all virtual vertices)
            std::vector<std::unique_ptr<KGVertex>>&  getAllVertices(){
                return this->allVertices;
            }

            //get the adjusted vertices in this knit graph (including all virtual vertices)
            std::vector<std::unique_ptr<KGVertex>>&  getAdjustedVertices(){
                return this->adjustedVertices;
            }

            

            //build the knit graph 
            void buildGraph();

};
