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
		Undefined, Interior, Opening, Wall, Mixed
	};

	Face();
	// Construct by ordered edge List. The edges can be non-straight.
	Face(const std::vector<Edge*> & faceEdges);
	// Construct by ordered vertex List. In this way the face does not contain any non-straight edge.
	Face(const std::vector<Vertex*> & faceVertices);
	~Face();

	const std::vector<Edge*>   getEdgeList() const;
	// Get <vertexList> which contains the topological vertices (namely all the end points of its edges).
	const std::vector<Vertex*> getVertexList() const;
	// Get <vertexListExtended> which contains both topological and auxiliary vertices. 
	const std::vector<Vertex*> getVertexListExtended() const;
	const std::vector<Cell*>   getCellList() const;
	const std::vector<Face*>   getSubFaceList() const;

	// Compare the edges of this face and the given edges. 
	// If they are exactly the same up to permutation, return True, else return False. 
	bool hasFaceEdges(std::vector<Edge*> faceEdges) const;
	
	friend std::ostream & operator<<(std::ostream & os, Face & obj);

	// If the inner orientation of the edge matches the induced orientation by face normal,
	// then retrun 1, else return -1
	const int getOrientation(const Edge* edge) const;

	bool isBoundary() const;
	bool hasBoundaryEdges() const;

	bool  hasCommonCell(Face* other) const;
	Cell* getCommonCell(Face* other) const;

	void setAdjacient(Cell* cell);



	const double getArea() const;
	const double getProjectedArea() const;
	const Eigen::Vector3d getCenter() const;
	const Eigen::Vector3d getNormal() const;
	/**
	 * @brief Compute the area of the face
	 * 		- For plane face, it is naturally defined;
	 * 		- For non-plane face, the area is the sum of the areas of the sub-faces
	 */
	const double computeArea();
	/**
	 * @brief Compute the center location of the face
	 * 		- For triangle, return the circumcenter;
	 * 		- For polygon, return the arithmetic average of end points. TODO: Need a better way.
	 * 		- For non-plane face, return the arithmetic average of end points. 
	 */
	const Eigen::Vector3d computeCenter();
	/**
	 * @brief Compute(Define) the normal vector. Note that <computeCenter> MUST be called in advance.
	 * 		- For plane face, it is well-defined and is determined by the order in which the vertices are stored.
	 * 		- For non-plane face, it is not well-defined but still computed in the same way as it is a plane face.
	 */
	const Eigen::Vector3d computeNormal();
	
	
	void  reverseNormal();

	void setFluidType(const FluidType & fluidType);
	const FluidType getFluidType() const;
	// Check if <edges> forms a loop without reordering. Note that it only checks a necessary condition.
	bool isContinuousLoop(std::vector<Edge*> edges) const;

private:
	std::vector<Edge*>   edgeList;
	std::vector<Vertex*> vertexList;
	std::vector<Vertex*> vertexListExtended;
	std::vector<Cell*>   cellList;
	std::vector<Face*>   subFaceList;

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
	/**
	 * @brief Add sub-faces for the non-plane face
	 * @note This function is very specialized to the domain geometry. For our case, we only have to deal with two cases:
	 * 			- two sub-faces
	 * 			- three sub-faces
	 */
	void addSubFaces();
	// Check if this face is plane. The criteria is to check if the face has only one or zero non-straight edge
	// which is case-specific!
	bool isPlane() const;
	bool isListOfVerticesUnique() const;
	bool isListOfEdgesUnique() const;

	const Eigen::Vector3d getCircumCenter() const;
	// !This function is not tailored for nonregular faces

};

