#include "Edge.h"



Edge::Edge()
{
}

Edge::Edge(const int index)
	: GeometryItem(index)
{
}

Edge::Edge(Vertex * A, Vertex * B)
{
	assert(A != nullptr);
	assert(B != nullptr);
	this->A = A;
	this->B = B;
	A->setAdjacient(this);
	B->setAdjacient(this);
}


Edge::~Edge()
{
}

// ytw: why two getVertex? How to decide which one is invoked?
Vertex * Edge::getVertexA()
{
	assert(A != nullptr);
	return A;
}

const Vertex * Edge::getVertexA() const
{
	return A;
}

Vertex * Edge::getVertexB()
{
	assert(B != nullptr);
	return B;
}

const Vertex * Edge::getVertexB() const
{
	return B;
}

const Eigen::Vector3d Edge::getDirection() const
{
	return B->getPosition() - A->getPosition();
}

const double Edge::getLength() const
{
	const double length = getDirection().norm();
	assert(length > 0);
	return length;
}

// ytw: Why two isAjacient?
bool Edge::isAdjacient(const Vertex * A, const Vertex * B) const
{
	assert(A != nullptr);
	assert(B != nullptr);
	assert(A != B);
	return (A == this->A && B == this->B) || (A == this->B && B == this->A);
}

bool Edge::isAdjacient(const Vertex * A, const Vertex * B) 
{
	assert(A != nullptr);
	assert(B != nullptr);
	assert(A != B);
	return (A == this->A && B == this->B) || (A == this->B && B == this->A);
}

void Edge::setAdjacient(Face * face)
{
	assert(face != nullptr);
	this->faceList.push_back(face);
}

bool Edge::hasConnectedFace(const std::vector<Edge*> & faceEdges) const
{
	for (auto face : faceList) {
		if (face->hasFaceEdges(faceEdges)) {
			return true;
		}
	}
	return false;
}

// ytw: obsolete?
Face * Edge::getConnectedFace(const std::vector<Edge*>& faceEdges) const
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

Vertex * Edge::getOppositeVertex(const Vertex * v) const
{
	assert(v != nullptr);
	assert(hasVertex(v));
	return (v == A) ? B : A;
}

bool Edge::hasVertex(const Vertex * v) const
{
	assert(v != nullptr);
	return v == A || v == B;
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
	assert(A != nullptr);
	assert(B != nullptr);
	return 0.5 * (A->getPosition() + B->getPosition());
}

Vertex * Edge::getCoincidentVertex(const Edge * other) const
{
	assert(other != nullptr);
	if (hasCoincidentVertex(other)) {
		Vertex * v = nullptr;
		bool cmp_AA = this->A == other->A;
		bool cmp_AB = this->A == other->B;
		bool cmp_BA = this->B == other->A;
		bool cmp_BB = this->B == other->B;

		if ((this->A == other->A) || (this->A == other->B))
		{
			v = A;
		}
		if ((this->B == other->A) || (this->B == other->B)) {
			v = B;
		}
		assert(v != nullptr);
		return v;
	}
	return nullptr;
}

bool Edge::hasCoincidentVertex(const Edge * other) const
{
	assert(other != nullptr);
	return (this->A == other->A || this->A == other->B || this->B == other->A || this->B == other->B);
}

int Edge::getIncidence(const Vertex * v) const
{
	assert(v != nullptr);
	assert(v == A || v == B);
	return (v == A) ? 1 : -1;
}

const std::vector<Face*> Edge::getFaceList() const
{
	return faceList;
}

void Edge::setType(const Edge::Type & type)
{
	this->type = type;
}

const Edge::Type Edge::getType() const
{
	return type;
}

std::ostream & operator<<(std::ostream & os, const Edge & obj)
{
	os << "Edge " << obj.getIndex() << ": ";
	os << "{" << obj.A->getIndex() << ", ";
	os << obj.B->getIndex() << "}";
	return os;
}
