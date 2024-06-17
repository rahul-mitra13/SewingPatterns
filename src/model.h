#include <vector>
#pragma once
class Model{

    private:

        //period for optimization
        double period;

        //timeout for optimization
        double timeOut;

        //stores a pair of (edge indices, integral value of the path that those edge indices specify) 
        std::vector<std::pair<std::vector<double>, double>> edgePathConstraints;

        //smooth face indices 
        std::vector<int> smoothFaces;

        //singular face indices
        std::vector<int> singularFaces;

        //boundary edges for a single patch 
        std::vector<int> bdyEdges;

        //these are the values we'll try to match using a quadratic minimization 
        std::vector<double> matchingTerms;

    public:

        //set optimization period 
        void setPeriod(double period);

        //set optimization timeout 
        void setTimeOut(double timeOut);

        //set edge path constraints
        void setEdgePathConstraints(std::vector<std::pair<std::vector<double>, double>> edgePathConstraints);

        //set singular face indices 
        void setSingularFaceIndices(std::vector<int> singularFaces);

        //set smooth face indices 
        void setSmoothFaceIndices(std::vector<int> smoothFaces);

        //set course boundary edges
        void setBdyEdges(std::vector<int> bdyEdges);

        //set matching terms 
        void setMatchingTerms(std::vector<double> matchingTerms);

        //get edge path constraints 
        std::vector<std::pair<std::vector<double>, double>> getEdgePathConstraints();

        //get period value 
        double getPeriod();

        //get timeout value
        double getTimeOut();

        //get smooth faces
        double getSmoothFaces();

        //get singular faces
        double getSingularFaces();

        //get course boundary edges 
        std::vector<int> getBdyEdges();

        //get matching energy
        std::vector<double> getMatchingTerms();

        //clear the data in the class 
        void clear();
};