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

using namespace geometrycentral;
using namespace geometrycentral::surface;

//redefining for consistency
struct KGVertex{
    int id = -1;
    Vector3 position;//position embedded in R^3 (if position information is available)
    Vector3 baryCoords;//barycentric coordinates of the vertex if position information is not available
    int row_in = -1;
    int row_out = -1;
    int col_in[2] = {-1, -1};
    int col_out[2] = {-1, -1};
    double alpha_tag = -1;
    double beta_tag = -1;
    Face face;//associated face of a vertex (on the underlying mesh)
    std::optional<Edge> edge;//associated edge of a vertex (for virtual vertices)
    std::optional<Halfedge> halfedge;//associated halfedge of a vertex (for virtual vertices)
    SurfacePoint surfacePoint;//surfacePoint representation of this knit graph vertex
    bool isAlphaVirtual = false;//is a virtual vertex in the course direction
    bool isBetaVirtual = false;//is a virtual vertex in the wale direction 

    void printVertexInfo(){
          
        std::cout << "vertex id = " << id << std::endl;
        std::cout << "row in = " << row_in << std::endl;
        std::cout << "row out = " << row_out << std::endl;
        std::cout << "col_in[0] = " << col_in[0] << std::endl;
        std::cout << "col_in[1] = " << col_in[1] << std::endl;
        std::cout << "col_out[0] = " << col_out[0] << std::endl;
        std::cout << "col_out[1] = " << col_out[1] << std::endl;
        std::cout << "isAlphaVirtual = " << isAlphaVirtual << std::endl;
        std::cout << "isBetaVirtual = " << isBetaVirtual << std::endl;
    }
};

class KG{

    private:
        
        //id counter of vertices 
        int vertexID = 0;

        //period for sampling the course stripes 
        double coursePeriod; 

        //period for sampling the wale stripes 
        double walePeriod;

        //geometry on which this graph lives 
        VertexPositionGeometry *globalGeometry;

        //edge length geometry if we want to do this in the embedded setting 
        EdgeLengthGeometry *gluedGeometry;

        //polyscope object for this geometry 
        polyscope::SurfaceMesh *psMesh; 

        //edge map from global mesh to glued mesh 
        std::map<int, int> *globalToGluedEdgeMap;

        // Flag of edges that were glued together
        EdgeData<bool> isGlued;

        //course stripe 1-form 
        CornerData<double> courseOneForm; 
        //wale stripe 1-form
        CornerData<double> waleOneForm;

        //edge indices in the course direction
        EdgeData<double> courseSingularEdges; 

        //edge indices in the wale direction
        EdgeData<double> waleSingularEdges;

        //vertices in the graph (including virtuals)
        std::vector<std::unique_ptr<KGVertex>> allVertices;

        std::vector<std::unique_ptr<KGVertex>> finalVertices;

        // Pairs of stitched (real) vertices
        std::vector<std::pair<int, int>> stitchedVertices;

        private: 

            //hashing floating point numbers 
            int hashFloat(double val);


            //get a knit graph vertex by id
            knitGraphVertex& get(int id);

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

            //get the vertices in this knit graph 
            std::vector<std::unique_ptr<KGVertex>>&  getAllVertices(){
                return this->allVertices;
            }

            //get the real vertices in this knit graph 
            std::vector<std::unique_ptr<KGVertex>>& getFinalVertices(){
                return this->finalVertices;
            }

            //build the knit graph 
            void buildGraph();

            //handle course singular edges  
            void handleCourseSingularEdges();

            //handle wale singular edges 
            void handleWaleSingularEdges();
};
