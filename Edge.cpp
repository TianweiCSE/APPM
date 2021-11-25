#include "Edge.h"



Edge::Edge()
{
}

Edge::Edge(const int index)
	: GeometryItem(index)
{
}

Edge::Edge(Vertex * A, Vertex * B) {
	assert(A != nullptr);
	assert(B != nullptr);
	assert(vertexList.size() == 0);
	vertexList.push_back(A);
	vertexList.push_back(B);
	A->setAdjacient(this);
	B->setAdjacient(this);
}

Edge::Edge(Edge* e1, Edge* e2) {
    Vertex* mid_v = e1->getCoincidentVertex(e2);
    vertexList = {e1->getOppositeVertex(mid_v), mid_v, e2->getOppositeVertex(mid_v)};
	e1->getOppositeVertex(mid_v)->setAdjacient(this);
	e2->getOppositeVertex(mid_v)->setAdjacient(this);
}

Edge::~Edge()
{
}

// ytw: why two getVertex? How to decide which one is invoked?
Vertex * Edge::getVertexA() const
{
	assert(vertexList.size() != 0);
	return vertexList.front();
}

Vertex * Edge::getVertexB() const
{
	assert(vertexList.size() != 0);
	return vertexList.back();
}

Vertex*  Edge::getVertexMid() const {
	if (vertexList.size() == 3) {
		return vertexList[1];
	}
	else {
		return nullptr;
	} 
}

const Eigen::MatrixXd Edge::getDirection() const
{	
	assert(getVertexMid() == nullptr && "Get direction of non-straight edge!"); // assert this edge is straight.
	return (getVertexB()->getPosition() - getVertexA()->getPosition()).normalized();
}

const double Edge::getLength() const
{
	double length = (vertexList[0]->getPosition() - vertexList[1]->getPosition()).norm();
	if (vertexList.size() == 3) {
		length +=  (vertexList[1]->getPosition() - vertexList[2]->getPosition()).norm();
	}
	assert(length > 0);
	return length;
}

bool Edge::isAdjacient(Vertex * A, Vertex * B) const
{
	assert(A != nullptr);
	assert(B != nullptr);
	assert(A != B);
	return (A == getVertexA() && B == getVertexB()) || (A == getVertexB() && B == getVertexA());
}

void Edge::setAdjacient(Face * face)
{
	assert(face != nullptr);
	this->faceList.push_back(face);
}

bool Edge::hasConnectedFace(const std::vector<Edge*> &faceEdges) const {
	for (auto face : faceList) {
		if (face->hasFaceEdges(faceEdges)) {
			return true;
		}
	}
	return false;
}

Face * Edge::getConnectedFace(const std::vector<Edge*> &faceEdges) const
{
	assert(hasConnectedFace(faceEdges));
	for (auto face : faceList) {
		if (face->hasFaceEdges(faceEdges)) {
			return face;
		}
	}
	assert(false);
	return nullptr;
}

Vertex * Edge::getOppositeVertex(Vertex * v) const
{
	assert(v != nullptr);
	assert(hasVertex(v));
	return (v == getVertexA()) ? getVertexB() : getVertexA();
}

bool Edge::hasVertex(Vertex * v) const
{
	assert(v != nullptr);
	return v == getVertexA() || v == getVertexB();
}

bool Edge::isBoundary() const
{
	if (faceList.size() <= 1) {
		return true;
	}
	else {
		for (auto face : faceList) {
			if (face->isBoundary()) {
				return true;
			}
		}
		return false;
	}
}

const Eigen::Vector3d Edge::getHalfwayPosition() const
{
	assert(vertexList.size() != 0);
	if (vertexList.size() == 2) {
		return 0.5 * (getVertexA()->getPosition() + getVertexB()->getPosition());
	}
	else {
		return vertexList[1]->getPosition();
	}
}

Vertex * Edge::getCoincidentVertex(Edge * other) const
{
	assert(other != nullptr);
	assert(hasCoincidentVertex(other));
	Vertex * v = nullptr;

	if ((getVertexA() == other->getVertexA()) || (getVertexA() == other->getVertexB()))
	{
		v = getVertexA();
	}
	if ((getVertexB() == other->getVertexA()) || (getVertexB() == other->getVertexB())) {
		v = getVertexB();
	}
	assert(v != nullptr);
	return v;
}

bool Edge::hasCoincidentVertex(Edge * other) const
{
	assert(other != nullptr);
	return (getVertexA() == other->getVertexA() || 
			getVertexA() == other->getVertexB() || 
			getVertexB() == other->getVertexA() || 
			getVertexB() == other->getVertexB());
}

// Do NOT change this function, otherwise the mesh generation would fail.
int Edge::getIncidence(Vertex * v) const
{
	assert(v != nullptr);
	assert(hasVertex(v));
	return (v == getVertexB()) ? -1 : 1;
}

std::vector<Face*> Edge::getFaceList() const
{
	return faceList;
}

void Edge::setType(const Edge::Type & type) {
	this->type = type;
}

const Edge::Type Edge::getType() const {
	return type;
}

std::ostream & operator<<(std::ostream & os, const Edge & obj)
{
	os << "Edge " << obj.getIndex() << ": ";
	os << "{" << obj.getVertexA()->getIndex() << ", ";
	os << obj.getVertexB()->getIndex() << "}";
	return os;
}
