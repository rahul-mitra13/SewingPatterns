#include "model.h"

void Model::setPeriod(double period){

    this->period = period;
}

void Model::setIntegrabilityConstraint(bool flag){

    this->integrabilityContraint = flag;
}

void Model::setTimeOut(double timeOut){
    
    this->timeOut = timeOut;
}

void Model::setEdgePathConstraints(std::vector<std::pair<std::vector<double>, double>>& edgePathConstraints){

    this->edgePathConstraints = edgePathConstraints;
}

void Model::setSingularFaceIndices(std::vector<std::pair<int, int>>& singularFaces){

    this->singularFaces = singularFaces;
}

void Model::setSingularVertexIndices(std::vector<std::pair<int, int>>& singularVertices){

    this->singularVertices = singularVertices;
}

void Model::setSingularEdges(std::vector<std::pair<int, int>> singularEdges){

    this->singularEdges = singularEdges;
}

void Model::setSmoothFaceIndices(std::vector<int>& smoothFaces){

    this->smoothFaces = smoothFaces;
}

void Model::setFaceIndices(std::vector<int>& faceIndices){

    this->faceIndices = faceIndices;

}

void Model::setBdyEdges(std::vector<int>& bdyEdges){

    this->bdyEdges = bdyEdges;
}

void Model::setMatchingTerms(std::vector<double>& matchingTerms){

    this->matchingTerms = matchingTerms;

}

void Model::setFaceGradients(std::vector<std::array<double, 3>>& faceGradients){
    
    this->faceGradients = faceGradients;

}

void Model::setFaceGradients(FaceData<Vector3>& faceGradients){

    for (int i = 0; i < faceGradients.size(); i++){
        this->faceGradients[i] = std::array<double, 3>{faceGradients[i][0], faceGradients[i][1], faceGradients[i][2]};
    }
}

void Model::setEdgeMappingsPairs(std::vector<std::pair<int, int>>& edgeMappingsPairs){

    this->edgeMappingsPairs = edgeMappingsPairs;

}

void Model::setWaleBdyPathConstraints(std::vector<std::vector<double>>& waleBdyPathConstraints){

    this->waleBdyPathConstraints = waleBdyPathConstraints;

}

double Model::getPeriod(){

    return this->period;
}

double Model::getTimeOut(){

    return this->timeOut;
}

bool Model::getIntegrabilityConstraint(){
    
    return this->integrabilityContraint;
}

std::vector<int> Model::getSmoothFaces(){

    return this->smoothFaces;
}

std::vector<std::pair<int, int>> Model::getSingularFaces(){

    return this->singularFaces;

}

std::vector<std::pair<int, int>> Model::getSingularVertices(){

    return this->singularVertices;

}

std::vector<std::pair<int, int>> Model::getSingularEdges(){

    return this->singularEdges;
}

std::vector<int> Model::getFaceIndices(){

    return this->faceIndices;

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

std::vector<std::array<double, 3>> Model::getFaceGradients(){

    return this->faceGradients;

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
    waleBdyPathConstraints.clear();
}



