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

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct knitGraphVertex{
    int id;
    Vector3 position;//position embedded in R^3 (if position information is available)
    Vector3 baryCoords;//barycentric coordinates of the vertex if position information is not available
    int row_in = -1;
    int row_out = -1;
    int col_in[2] = {-1, -1};
    int col_out[2] = {-1, -1};
    double alpha_tag = -1;
    double beta_tag = -1;
    Face face;//associated face of a vertex (on the underlying mesh)
    SurfacePoint surfacePoint;//used for intersections on singular faces
    bool isBaryCenter = false;//if this is a vertex at the barycenter
    bool hasBeenHandled = false;//if this vertex has been handled by a merge
    bool isVirtual = false;//if this is a virtual vertex
    bool isAlphaVirtual = false;//is a virtual vertex in the course direction
    bool isBetaVirtual = false;//is a virtual vertex in the wale direction 
};

class KnitGraph{

    private:
        
        //id counter of vertices 
        int vertexID = 0;

        //period for sampling the knit graph 
        double period; 

        //geometry on which this graph lives 
        VertexPositionGeometry *globalGeometry;

        //edge length geometry if we want to do this in the embedded setting 
        EdgeLengthGeometry *gluedGeometry;

        //polyscope object for this geometry 
        polyscope::SurfaceMesh *psMesh; 

        //edge map from global mesh to glued mesh 
        std::map<int, int> *globalToGluedEdgeMap;

        //course stripe 1-form 
        CornerData<double> courseOneForm; 
        //wale stripe 1-form
        CornerData<double> waleOneForm;

        //course singular edges (in the glued setting)
        EdgeData<double> courseSingularEdges; 

        //wale singular edges (in the glued setting)
        EdgeData<double> waleSingularEdges;

        //vertices in the graph 
        std::vector<knitGraphVertex> vertices;

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

            //Constructor
            KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, double period, 
                      CornerData<double>& courseOneForm, EdgeData<double>& courseSingularEdges, CornerData<double>& waleOneForm, EdgeData<double>& waleSingularEdges,
                      std::map<int, int>& globalToGluedEdgeMap);

            //build the knit graph 
            void buildGraph();

            //find knit graph vertices on faces that are smooth in both directions
            void handleCourseNonSingularFaceWaleNonSingularFace(Face &f);

            //update connections on the vertices of a smooth face 
            void connectOnSmoothFace(std::vector<knitGraphVertex>& faceVertices);

            //render the knit graph 
            void renderGraph();

            //mkae obj out of knit graph 
            void makeObj();
};

