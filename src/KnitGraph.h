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

struct knitGraphVertex{
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
    std::optional<Edge> edge;//associated edge of a vertex
    std::optional<Halfedge> halfedge;//associated halfedge of a vertex
    SurfacePoint surfacePoint;//used for intersections on singular faces
    bool isBaryCenter = false;//if this is a vertex at the barycenter
    bool hasBeenHandled = false;//if this vertex has been handled by a merge
    bool isVirtual = false;//if this is a virtual vertex
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
        std::cout << "isVirtual = " << isVirtual << std::endl;
    }
};


class KnitGraph{

    private:
        
        //id counter of vertices 
        int vertexID = 0;

        //period for sampling the knit graph (unused for now)
        // double period; 

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

        // Pairs of stitched (real) vertices
        std::vector<std::pair<int, int>> stitchedVertices;

        //course stripe 1-form 
        CornerData<double> courseOneForm; 
        //wale stripe 1-form
        CornerData<double> waleOneForm;

        //course singular edges (in the glued setting)
        EdgeData<double> courseSingularEdges; 

        //wale singular edges (in the glued setting)
        EdgeData<double> waleSingularEdges;

        //vertices in the graph (including virtuals)
        std::vector<knitGraphVertex> vertices;

        std::vector<knitGraphVertex> realVertices;

        //vertex info map 
        //vertex id to information
        std::map<int, knitGraphVertex> vertexInfoMap;



        //vertex positions of this knit graph 
        std::vector<Vector3> vertexPositions; 
        //edges in this knit graph 
        std::vector<std::array<int, 2>> edges; 
        //faces in this knit graph 
        std::vector<std::vector<int>> faces;//or std::vector<std::vector<std::array<int, 4>> faces?

        private: 

            //hashing floating point numbers 
            int hashFloat(double val);


            //get a knit graph vertex by id
            knitGraphVertex& get(int id);

        public: 

            //Default Constructor 
            KnitGraph(){};

            //Constructor
            KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, double coursePeriod, double walePeriod,
                      CornerData<double>& courseOneForm, EdgeData<double>& courseSingularEdges, CornerData<double>& waleOneForm, EdgeData<double>& waleSingularEdges,
                      std::map<int, int>& globalToGluedEdgeMap);

            //get the vertices in this knit graph 
            std::vector<knitGraphVertex>& getVertices(){
                return this->vertices;
            }

            //get the real vertices in this knit graph 
            std::vector<knitGraphVertex>& getRealVertices(){
                return this->realVertices;
            }

            //build the knit graph 
            void buildGraph();

            //find knit graph vertices on faces that are smooth in both directions
            void handleCourseNonSingularFaceWaleNonSingularFace(Face &f);

            //update connections on the vertices of a smooth face 
            void connectOnSmoothFace(std::vector<knitGraphVertex>& faceVertices);

            //render the knit graph 
            void renderGraph();

            //handle the merge in the instrinsic setting 
            void intrinsicMerge();

            //perform epsilon merging to 
            //merge connections across faces
            void epsilonMerging();

            //reorder indices of knit graph
            void makeRealVertices();

            //find a cluster of vertices to merge
            std::vector<knitGraphVertex> findCluster(const knitGraphVertex &v, double eps);

            //merge a cluster
            void mergeCluster(std::vector<knitGraphVertex>& vCluster, double eps);

            //tag increases and decreases
            void tagIncreasesDecreases();

            //trace short rows to check for helices 
            void traceShortRows();

            //sanity check the graph 
            void sanityCheck();

            //perform real merges on misaligned edges in the knit graph 
            void realMerge();

            //write knit graph to txt file
            void writeKnitGraphToTxtFile(const std::string& file_name);

            //write knit graph as line element
            void writeKnitGraphLineElement();
};

