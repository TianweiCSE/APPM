#pragma once

#include "GeometryItem.h"
#include "Edge.h"
class Edge;   //ytw: what if no such line?

class Vertex
	: public GeometryItem
{
public:

	enum class Type {
		Undefined, Interior, Electrode, Insulating
	};

	Vertex();
	Vertex(const int index);
	Vertex(const Eigen::Vector3d & position);
	Vertex(const Eigen::Vector3d & position, const int index);
	~Vertex();

	void setPosition(const Eigen::Vector3d & pos);
	const Eigen::Vector3d getPosition() const;

	void setAdjacient(Edge * edge);
	bool isAdjacientTo(Vertex * other);
	Edge* getAdjacientEdge(Vertex * other);

	const std::vector<Edge*> getEdges() const;

	friend std::ostream & operator<<(std::ostream & os, Vertex & obj);

	bool isBoundary() const;

	void setType(const Type & type);
	Type getType() const;

	// return a new vertex that has the same location
	Vertex* copy() const; 

private:
	Eigen::Vector3d position;
	std::vector<Edge*> adjacientEdges;
	Type type = Type::Undefined;
};

