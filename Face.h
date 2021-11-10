#pragma once
#include "GeometryItem.h"
#include "Edge.h"
class Edge;
#include "Cell.h"
class Cell;

class Face :
	public GeometryItem
{
public:

	enum class FluidType {
		Undefined, Interior, Opening, Wall
	};

	Face();
	// Construct by ordered edge List. The edges can be non-straight.
	Face(const std::vector<Edge*> & faceEdges);
	// Construct by ordered vertex List. In this way the face does not contain any non-straight edge.
	Face(const std::vector<Vertex*> & faceVertices);
	~Face();

	std::vector<Edge*>   getEdgeList() const;
	// Get <vertexList> which contains the topological vertices (namely all the end points of its edges).
	std::vector<Vertex*> getVertexList() const;
	// Get <vertexListExtended> which contains both topological and auxiliary vertices. 
	std::vector<Vertex*> getVertexListExtended() const;
	std::vector<Cell*>   getCellList() const;

	// Compare the edges of this face and the given edges. 
	// If they are exactly the same up to permutation, return True, else return False. 
	bool hasFaceEdges(std::vector<Edge*> faceEdges) const;
	
	friend std::ostream & operator<<(std::ostream & os, Face & obj);

	// If the inner orientation of the edge matches the induced orientation by face normal,
	// then retrun 1, else return -1
	const int getOrientation(Edge* edge) const;

	bool isBoundary() const;
	bool hasBoundaryEdges() const;

	bool  hasCommonCell(Face* other) const;
	Cell* getCommonCell(Face* other) const;

	void setAdjacient(Cell* cell);


	/** The area, center, and normal are defined as follows 
	 *      - for plane face, they are natrually defined.
	 *      - for non-plane face, they are defined with respect to the projected face, 
	 *            i.e. the face defined by topological vertices in <vertexList>.
	 */ 
	const double          getArea()   const;
	const Eigen::Vector3d getCenter() const;
	const Eigen::Vector3d getNormal() const;
	const double          computeArea();
	const Eigen::Vector3d computeNormal();
	const Eigen::Vector3d computeCenter();
	
	void  reverseNormal();

	void setFluidType(const FluidType & fluidType);
	const FluidType getFluidType() const;
	// Check if <edges> forms a loop without reordering. Note that it only checks a necessary condition.
	bool isContinuousLoop(std::vector<Edge*> edges) const;

private:
	std::vector<Edge*> edgeList;
	std::vector<Vertex*> vertexList;
	std::vector<Vertex*> vertexListExtended;
	std::vector<Cell*> cellList;

	Eigen::Vector3d center;  /// circumcenter for triangles; average of vertices for polygons!
	                            // !This variable is not well defined for nonregular faces
	Eigen::Vector3d faceNormal; // !This variable is not well defined for nonregular faces
	double area;            // !This variable is not well defined for nonregular faces

	FluidType fluidType = FluidType::Undefined;

	/** 
	 *  - build edgeList/vertexList. The order is determined by the constructor input.
	 * 	- compute the center which is circumcenter for triangle and arithmetic center for polygon
	 *  - compute the face normal which forms a right-handed system with the order of <vertexList>
	 *  - compute the area 
	 *  - register adjacency to its edges
	 */
	void init();  
	bool isListOfVerticesUnique() const;
	bool isListOfEdgesUnique() const;

	const Eigen::Vector3d getCircumCenter() const;
	// !This function is not tailored for nonregular faces

};

