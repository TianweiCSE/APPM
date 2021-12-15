#include "Vertex.h"



Vertex::Vertex()
{
}

Vertex::Vertex(const int index)
	: GeometryItem(index)
{
}

Vertex::Vertex(const Eigen::Vector3d & position)
{
	this->position = position;
}

Vertex::Vertex(const Eigen::Vector3d & position, const int index)
	: Vertex(position)
{
	setIndex(index);
}


Vertex::~Vertex()
{
}



void Vertex::setPosition(const Eigen::Vector3d & pos)
{
	this->position = pos;
}

const Eigen::Vector3d Vertex::getPosition() const
{
	return this->position;
}

void Vertex::setAdjacient(Edge * edge) 
{
	assert(edge != nullptr);
	adjacientEdges.push_back(edge);
}

bool Vertex::isAdjacientTo(Vertex * other) 
{
	assert(other != nullptr);
	for (auto edge : adjacientEdges) {
		if (edge->isAdjacient(this, other)) {
			return true;
		}
	}
	return false;
}

Edge * Vertex::getAdjacientEdge(Vertex * other)
{
	assert(other != nullptr);
	for (auto edge : adjacientEdges) {
		if (edge->isAdjacient(this, other)) {
			return edge;
		}
	}
	return nullptr;
}

const std::vector<Edge*> Vertex::getEdges() const
{
	return adjacientEdges;
}

bool Vertex::isBoundary() const
{
	for (auto edge : adjacientEdges) {
		if (edge->isBoundary()) {
			return true;
		}
	}
	return false;
}

void Vertex::setType(const Vertex::Type & type)
{
	this->type = type;
}

const Vertex::Type Vertex::getType() const
{
	return this->type;
}

Vertex* Vertex::copy() const {
	assert(this != nullptr);
	Vertex* v = new Vertex(getPosition());
	return v;
} 

std::ostream & operator<<(std::ostream & os, const Vertex & obj)
{
	os << "Vertex idx: " << obj.getIndex() << "; ";
	os << "edges = {";
	for (auto e : obj.getEdges()) {
		os << e->getIndex() << ",";
	}
	os << "}";
	return os;
}
