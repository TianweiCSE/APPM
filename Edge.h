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
		Interior, InteriorToBoundary, Boundary
	};

	Edge();
	Edge(const int index);
	Edge(Vertex * A, Vertex * B);
	~Edge();

	Vertex * getVertexA();
	const Vertex * getVertexA() const;
	Vertex * getVertexB();
	const Vertex * getVertexB() const;
	virtual const Eigen::MatrixXd getDirection() const;  

	/// return vector BA
	virtual const double getLength() const;  

	bool isAdjacient(const Vertex * A, const Vertex * B) const;
	bool isAdjacient(const Vertex * A, const Vertex * B);

	/// Add face to faceList
	void setAdjacient(Face * face); 

	/// return true if the face represented by faceEdges is in face list
	virtual bool hasConnectedFace(const std::vector<Edge*> & faceEdges) const; 

	Face * getConnectedFace(const std::vector<Edge*> & faceEdges) const; /// incomplete

	friend std::ostream & operator<<(std::ostream & os, const Edge & obj); 

	Vertex * getOppositeVertex(const Vertex * v) const;
	bool hasVertex(const Vertex * v) const;
	bool isBoundary() const;

	/// return the mid-point
	virtual const Eigen::Vector3d getHalfwayPosition() const;

	Vertex * getCoincidentVertex(const Edge * other) const; 
	bool hasCoincidentVertex(const Edge * other) const;

	/// return 1 if v is A, else return -1
	int getIncidence(const Vertex * v) const; 
	const std::vector<Face*> getFaceList() const;

	void setType(const Edge::Type & type);
	const Edge::Type getType() const;


protected:
	Vertex * A = nullptr;
	Vertex * B = nullptr;
	Eigen::Vector3d edgeCenter;
	std::vector<Face*> faceList;
	Type type;
};

