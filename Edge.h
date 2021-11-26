#pragma once
#include "GeometryItem.h"
#include "Vertex.h"
class Vertex;

#include "Face.h"
class Face;   

class Edge :
	public GeometryItem
{
public:

	enum class Type {
		Undefined, Interior, Electrode, Insulating
	};

	Edge();
	Edge(const int index);
	Edge(Vertex * A, Vertex * B);
	Edge(Edge* e1, Edge* e2);
	~Edge();

	// Get the first vertex in <vertexList>
	Vertex*  getVertexA()   const;
	// Get the last vertex in <vertexList>
	Vertex*  getVertexB()   const;
	// Get <vertexList[1]> if this edge is not straight, else return nullptr.
	// This function can be invoked from outside to check if <this> is an irregular edge.
	Vertex*  getVertexMid() const;

	std::vector<Face*>  getFaceList() const;

	// Get unit direction of vector AB. Straight edge is asserted first. 
	const Eigen::Vector3d getDirection() const;  
	// Get length
	const double getLength()    const;  

	// Return true if A and B forms an edge
	bool isAdjacient(Vertex * A, Vertex * B) const;

	// Add face to faceList
	void setAdjacient(Face * face); 

	// Return true if the face represented by <faceEdges> is in face list
	bool hasConnectedFace(const std::vector<Edge*> & faceEdges) const; 
	// Get the adjacent face represented by <faceEdges>. <hasConnectedFace> is asserted first.
	Face * getConnectedFace(const std::vector<Edge*> & faceEdges) const; 

	friend std::ostream & operator<<(std::ostream & os, const Edge & obj); 

	Vertex * getOppositeVertex(Vertex * v) const;
	bool hasVertex(Vertex * v) const;
	bool isBoundary() const;

	// - For straight edge: return the midpoint position;
	// - For non-straight edge: return the poisiton of <getVertexMid()>, namely the joint vertex.
	const Eigen::Vector3d getHalfwayPosition() const;

	// Return true if <this> has a common vertex with <other>. 
	bool hasCoincidentVertex(Edge * other) const;
	// Return the common vertex of <this> and <other>. <hasCoincidentVertex> is asserted first.
	Vertex * getCoincidentVertex(Edge * other) const; 
	
	// Return 1 if v is A, else return -1. 
	int getIncidence(Vertex * v) const; 
	
	void setType(const Edge::Type & type);
	const Edge::Type getType() const;

private:
	std::vector<Vertex *> vertexList;
	std::vector<Face*> faceList;
	Type type = Type::Undefined;
};

