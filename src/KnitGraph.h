//geometry-central includes 
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/surface_point.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

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
        //period for sampling the knit graph 
        double period; 

        //geometry on which this graph lives 
        VertexPositionGeometry *globalGeometry;

        //edge length geometry if we want to do this in the embedded setting 
        EdgeLengthGeometry *gluedGeometry;

        //polyscope object for this geometry 
        polyscope::SurfaceMesh *psMesh; 
    
        //course stripe 1-form 
        EdgeData<double> courseOneForm; 
        //wale stripe 1-form
        EdgeData<double> waleOneForm;

        //vertices in the graph 
        std::vector<knitGraphVertex> vertices;

        //vertex positions of this knit graph 
        std::vector<Vector3> vertexPositions; 
        //edges in this knit graph 
        std::vector<std::array<int, 2>> edges; 
        //faces in this knit graph 
        std::vector<std::vector<int>> faces;//or std::vector<std::vector<std::array<int, 4>> faces?

        public: 

            //Constructor
            KnitGraph(VertexPositionGeometry& globalGeometry, EdgeLengthGeometry& gluedGeometry, polyscope::SurfaceMesh& psMesh, double period, EdgeData<double>& courseOneForm, EdgeData<double>& waleOneForm);

            //build the knit graph 
            void buildGraph();
};

