#include "model.h"

void Model::setPeriod(double period){

    this->period = period;
}

void Model::setTimeOut(double timeOut){
    
    this->timeOut = timeOut;
}

void Model::setEdgePathConstraints(std::vector<std::pair<std::vector<double>, double>> edgePathConstraints){

    this->edgePathConstraints = edgePathConstraints;
}

void Model::setSingularFaceIndices(std::vector<int> singularFaces){

    this->singularFaces = singularFaces;
}

void Model::setSmoothFaceIndices(std::vector<int> smoothFaces){

    this->smoothFaces = smoothFaces;
}

void Model::setCourseBdyEdges(std::vector<int> courseBdyEdges){

    this->courseBdyEdges = courseBdyEdges;
}

void Model::setMatchingTerms(std::vector<double> matchingTerms){

    this->matchingTerms = matchingTerms;

}

double Model::getPeriod(){

    return this->period;
}

double Model::getTimeOut(){

    return this->timeOut;
}

std::vector<std::pair<std::vector<double>, double>> Model::getEdgePathConstraints(){

    return this->edgePathConstraints;
}

std::vector<int> Model::getCourseBdyEdges(){
    
    return this->courseBdyEdges;

}

std::vector<double> Model::getMatchingTerms(){

    return this->matchingTerms;

}



