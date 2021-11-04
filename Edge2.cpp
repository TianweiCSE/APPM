#include "Edge2.h"

Edge2::Edge2(){

}

Edge2::Edge2(const Edge* e1, const Edge* e2):_e1(e1), _e2(e2){
    edgeCenter = getHalfwayPosition();
    type = Edge::Type::Boundary;
    Vertex* mid_v = _e1->getCoincidentVertex(_e2);
    A = _e1->getOppositeVertex(mid_v);
    B = _e2->getOppositeVertex(mid_v); 
}

Edge2::~Edge2(){

}

const Eigen::MatrixXd Edge2::getDirection() const{
    Vertex* mid_v = _e1->getCoincidentVertex(_e2);
    Eigen::Matrix<double, 3, 2> directions;
    directions.col(0) = mid_v->getPosition() - A->getPosition();
    directions.col(1) = B->getPosition() - mid_v->getPosition();
    return directions;
}

const double Edge2::getLength() const{
    return _e1->getLength() + _e2->getLength();
}

bool Edge2::hasConnectedFace(const std::vector<Edge*> & faceEdges) const{
    return true;
}

const Eigen::Vector3d Edge2::getHalfwayPosition() const{
    return _e1->getCoincidentVertex(_e2)->getPosition();
}