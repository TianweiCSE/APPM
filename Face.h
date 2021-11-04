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
		DEFAULT, INTERIOR, OPENING, WALL
	};
	//enum class EmagType {
	//	DEFAULT
	//};

	Face();
	// Construct by ordered edgeList
	Face(const std::vector<Edge*> & faceEdges);
	// Construct by ORDERED vertexList
	Face(const std::vector<Vertex*> & faceVertices);
	~Face();

	std::vector<Edge*> getEdgeList() const;
	std::vector<Vertex*> getVertexList() const;
	std::vector<Cell*> getCellList() const;

	/// Compare the edges of this face and the given edges. 
	/// If they are exactly the same up to permutation, return True, else return False. 
	bool hasFaceEdges(const std::vector<Edge*> faceEdges) const;
	
	friend std::ostream & operator<<(std::ostream & os, const Face & obj);

	/// If the inner orientation of the edge matches the induced orientation by face normal,
	/// then retrun 1, else return -1
	const int getOrientation(const Edge * edge);

	bool isBoundary() const;
	bool hasBoundaryEdges() const;

	const Eigen::Vector3d getCenter() const;
	const double getArea() const;

	bool hasCommonCell(const Face * other) const;
	Cell * getCommonCell(const Face * other) const;

	void setAdjacient(Cell * cell);

	const Eigen::Vector3d getNormal() const;
	void setNormal(const Eigen::Vector3d & fn);

	void setFluidType(const FluidType & fluidType);
	const FluidType getFluidType() const;


private:
	std::vector<Edge*> edgeList;
	std::vector<Vertex*> vertexList;
	std::vector<Cell*> cellList;

	Eigen::Vector3d center;  /// circumcenter for triangles; average of vertices for polygons!
	                            // !This variable is not well defined for nonregular faces
	Eigen::Vector3d faceNormal; // !This variable is not well defined for nonregular faces
	double area = 0;            // !This variable is not well defined for nonregular faces

	FluidType fluidType = FluidType::DEFAULT;

	/** 
	 * - build edgeList/vertexList. The order is determined by the constructor input.
	 *  - register adjacency to its edges
	 *  - compute the center which is circumcenter for triangle and arithmetic center for polygon
	 *  - compute the face normal which is ONLY determined by the internal orientation of the first edge
	 */
	void init();  
	bool isListOfVerticesUnique() const;
	bool isListOfEdgesUnique() const;

	const Eigen::Vector3d getCircumCenter() const;
	// !This function is not tailored for nonregular faces
	const double computeArea() const;

};

