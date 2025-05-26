#include <vector>
#pragma once

//geometry central includes
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Model{

    private:

        //period for optimization
        double period;

        //timeout for optimization
        double timeOut;

        //if we want to specify the integrability constraint or not 
        bool integrabilityContraint;

        //stores a pair of (edge indices, integral value of the path that those edge indices specify) 
        //the length of the vector is mesh->nEdges()
        std::vector<std::pair<std::vector<double>, double>> edgePathConstraints;


        //stores a pair of (edge indices, integral value of the path that those edge indices specify) 
        //the length of the vector is mesh->nHalfedges()
        std::vector<std::pair<std::vector<double>, double>> halfedgePathConstraints;


        //homology generator constraints 
        std::vector<std::vector<double>> homologyGenerators;

        //smooth face indices 
        std::vector<int> smoothFaces;

        //singular face indices
        std::vector<std::pair<int, int>> singularFaces;

        //singular vertex indices 
        std::vector<std::pair<int, int>> singularVertices;

        //edge singularities 
        //vector of pairs where each entry is pair(a, b)
        //a - index of the edge 
        //b - constrained value of the edge
        std::vector<std::pair<int , int>> singularEdges; 

        //edge singularities 
        std::vector<int> edgeIndices;

        //specify all face indices (both smooth and singular)
        std::vector<int> faceIndices;

        //boundary edges where we set the value of sigma to 0
        std::vector<int> bdyEdges;

        //boundary loops such that the value of the integral along the boundary = kP (wale direction constraint) along the wale direction
        std::vector<std::vector<double>> waleBdyPathConstraints;


        //edge pairs where the value of sigma should agree 
        std::vector<std::pair<int, int>> edgeMappingsPairs;

        //these are the values we'll try to match using a quadratic minimization 
        //these are edge based matching
        std::vector<double> matchingTerms;

        //these are the gradients over the faces that we're trying to match
        std::vector<std::array<double, 3>> faceGradients; 

    public:

        //flags for debugging purposes and using a different comparison in the objective 
        bool useEdgeAveraging = false;
        bool useFaceDifferenceViz = false;
        bool useHelicingCorrection = true;
        bool forceSaddle = false;//force the saddle vertices
           
        //set optimization period 
        void setPeriod(double period);

        //set integrability constraint 
        void setIntegrabilityConstraint(bool flag);

        //set optimization timeout 
        void setTimeOut(double timeOut);

        //set edge path constraints
        void setEdgePathConstraints(std::vector<std::pair<std::vector<double>, double>>& edgePathConstraints);

        //set halfedge path constraints
        void setHalfedgePathConstraints(std::vector<std::pair<std::vector<double>, double>>& edgePathConstraints);

        //set singular face indices 
        void setSingularFaceIndices(std::vector<std::pair<int, int>>& singularFaces);

        //set singular vertex indices 
        void setSingularVertexIndices(std::vector<std::pair<int, int>>& singularVertices);

        //set singular edge indices 
        void setSingularEdges(std::vector<std::pair<int, int>>& singularEdges);

        //set edge indices 
        void setEdgeIndices(std::vector<int>& edgeIndices);

        //set smooth face indices 
        void setSmoothFaceIndices(std::vector<int>& smoothFaces);

        //set all face indices
        void setFaceIndices(std::vector<int>& faceIndices);

        //set course boundary edges
        void setBdyEdges(std::vector<int>& bdyEdges);

        //set matching terms 
        void setMatchingTerms(std::vector<double>& matchingTerms);

        //set face gradients 
        void setFaceGradients(std::vector<std::array<double, 3>>& faceGradients);

        //set face gradients using gc FaceData<Vector3>
        void setFaceGradients(FaceData<Vector3>& faceGradients);

        //set homology generators
        void setHomologyGenerators(std::vector<std::vector<double>>& homologyGenerators);

        //set edge mappings 
        void setEdgeMappingsPairs(std::vector<std::pair<int, int>>& edgeMappingsPairs);

        //set wale path constraints 
        void setWaleBdyPathConstraints(std::vector<std::vector<double>>& waleBdyPathConstraints);

        //get edge path constraints 
        std::vector<std::pair<std::vector<double>, double>> getEdgePathConstraints();

        //get halfedge path constraints 
        std::vector<std::pair<std::vector<double>, double>> getHalfedgePathConstraints();

        //get period value 
        double getPeriod();

        //get integrability constraint
        bool getIntegrabilityConstraint();

        //get timeout value
        double getTimeOut();

        //get smooth faces
        std::vector<int> getSmoothFaces();

        //get singular faces
        std::vector<std::pair<int, int>> getSingularFaces();

        //get singular vertices 
        std::vector<std::pair<int, int>> getSingularVertices();

        //get singular edges 
        std::vector<std::pair<int, int>> getSingularEdges();

        //get edge indices 
        std::vector<int> getEdgeIndices();

        //get all face indices
        std::vector<int> getFaceIndices();
        
        //get course boundary edges 
        std::vector<int> getBdyEdges();

        //get matching energy
        std::vector<double> getMatchingTerms();

        //get face gradients
        std::vector<std::array<double, 3>> getFaceGradients();

        //get pair of mapped edges 
        std::vector<std::pair<int, int>> getEdgeMappingsPairs();

        //get wale path constraints
        std::vector<std::vector<double>> getWaleBdyPathConstraints();

        //get homology generators 
        std::vector<std::vector<double>> getHomologyGenerators();

        //clear the data in the class 
        void clear();
};