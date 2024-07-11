#include "model.h"

void Model::setPeriod(double period){

    this->period = period;
}

void Model::setTimeOut(double timeOut){
    
    this->timeOut = timeOut;
}

void Model::setEdgePathConstraints(std::vector<std::pair<std::vector<double>, double>>& edgePathConstraints){

    this->edgePathConstraints = edgePathConstraints;
}

void Model::setSingularFaceIndices(std::vector<int>& singularFaces){

    this->singularFaces = singularFaces;
}

void Model::setSmoothFaceIndices(std::vector<int>& smoothFaces){

    this->smoothFaces = smoothFaces;
}

void Model::setBdyEdges(std::vector<int>& bdyEdges){

    this->bdyEdges = bdyEdges;
}

void Model::setMatchingTerms(std::vector<double>& matchingTerms){

    this->matchingTerms = matchingTerms;

}

void Model::setEdgeMappingsPairs(std::vector<std::pair<int, int>>& edgeMappingsPairs){

    this->edgeMappingsPairs = edgeMappingsPairs;

}

void Model::setWaleBdyPathConstraints(std::vector<std::vector<double>> waleBdyPathConstraints){

    this->waleBdyPathConstraints = waleBdyPathConstraints;

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

std::vector<int> Model::getBdyEdges(){
    
    return this->bdyEdges;

}

std::vector<double> Model::getMatchingTerms(){

    return this->matchingTerms;

}

std::vector<std::pair<int, int>> Model::getEdgeMappingsPairs(){

    return this->edgeMappingsPairs;
}

std::vector<std::vector<double>> Model::getWaleBdyPathConstraints(){

    return this->waleBdyPathConstraints;

}

void Model::clear(){
    period = 0.0;
    timeOut = 0.0;
    edgePathConstraints.clear();
    smoothFaces.clear();
    singularFaces.clear();
    bdyEdges.clear();
    matchingTerms.clear();
    edgeMappingsPairs.clear();
}



