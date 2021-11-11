#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "H5Writer.h"
#include "H5Reader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EigenAuxiliaries.h"


#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "TriPrism.h"
#include "XmlElement.h"
#include "XdmfAttribute.h"
#include "XdmfDataItem.h"
#include "XdmfDomain.h"
#include "XdmfGeometry.h"
#include "XdmfGrid.h"
#include "XdmfRoot.h"
#include "XdmfTime.h"
#include "XdmfTopology.h"


class Mesh
{
public:
	struct FacetCounts {
		int nV_boundary;
		int nV_undefined;
		int nV_interior;
		int nV_electrode;
		int nV_insulating;
		int nE_boundary;
		int nE_undefined;
		int nE_interior;
		int nE_electrode;
		int nE_insulating;
		int nF_boundary;
		int nF_undefined;
		int nF_interior;
		int nF_opening;
		int nF_wall;
		int nC_undefined;
		int nC_fluid;
		int nC_solid;
	} facet_counts;

	Mesh();
	Mesh(const std::string & meshPrefix);
	~Mesh();

	/** - write .dat files: incidence map e2v, v2f, f2c, f2v
	 *  - write .h5 file:
	 * 		- about Vectex: position, position_extended, index, type
	 * 		- about Edge:   e2v(!XML form), index, type
	 *      - about Face:   f2v(!XML form), index, type, center, normal, area
	 *   	- about Cell:   c2v(!XML form), index, type, center, volumn
	 */    
	void writeToFile();

	/**  - prefix-mesh.xdmf   << root --> domain --> treeGrid --> (vertexGrid, edgeGrid, surfaceGrid)
	 *   - prefix-volume.xdmf << root --> domain --> VolumeGrid
	 */
	void writeXdmf();

	Vertex *   addVertex(const Eigen::Vector3d & position);
	Edge *     addEdge(Vertex * A, Vertex * B);
	Face *     addFace(const std::vector<Edge*> & faceEdges);
	Face *     addFace(const std::vector<Vertex*> & faceVertices);

	/// TriPrism not yet implemented
	TriPrism * addTriPrism(const std::vector<Face*> & sideFaces, Face * bottomFace, Face * topFace);
	Cell *     addCell(const std::vector<Face*> & cellFaces);
	
	Vertex * getVertex(const int index) const;
	Edge *   getEdge(const int index) const;
	Edge *   getEdge(Vertex * A, Vertex * B) const;
	Face *   getFace(const int index) const;
	Face *   getFace(const std::vector<Edge*> & faceEdges);
	Cell *   getCell(const int index) const;
	Cell *   getCell(const std::vector<Face*> & cellFaces);

	const std::string getPrefix() const;

	const int getNumberOfVertices() const;
	const int getNumberOfEdges()    const;
	const int getNumberOfFaces()    const;
	const int getNumberOfCells()    const;

	// return xdmf vector containing [2 #vertices idx1 idx2 ... ]
	const std::vector<int> getXdmfTopology_edge2vertexIndices() const;
	// return xdmf vector containing [facetype #vertices idx1 idx2 idx3 .. ]
	const std::vector<int> getXdmfTopology_face2vertexIndices() const;
	// return xdmf vector containing [16 #faces #vertices idx1 idx2 idx3 .. #vertices idx1 idx2 idx3 .. ]
	const std::vector<int> getXdmfTopology_cell2vertexIndices() const;
	

	// return cellList
	const std::vector<Cell*>   getCells()    const;
	// return faceList
	const std::vector<Face*>   getFaces()    const;
	// return edgeList
	const std::vector<Edge*>   getEdges()    const;
	// return vertexList
	const std::vector<Vertex*> getVertices() const;

	// This function is obsolete
	void check() const;

	const Eigen::SparseMatrix<int> & get_f2eMap() const;
	const Eigen::SparseMatrix<int> & get_e2vMap() const;

	const Eigen::VectorXi getVertexTypes() const;
	const Eigen::VectorXi getEdgeTypes()   const;
	const Eigen::VectorXi getFaceTypes()   const;
	const Eigen::VectorXi getCellTypes()   const;

	const Eigen::VectorXi getVertexIndices() const;
	const Eigen::VectorXi getEdgeIndices()   const;
	const Eigen::VectorXi getFaceIndices()   const;
	const Eigen::VectorXi getCellIndices()   const;

protected:
	std::vector<Vertex*> vertexList;
	std::vector<Edge*>   edgeList;
	std::vector<Face*>   faceList;
	std::vector<Cell*>   cellList;

	Eigen::Matrix3Xd vertexCoordinates;

	// Check if <edges> can be reordered into a continuous loop. Note that it only checks a necessary condition.
	bool canBeContinuousLoop(std::vector<Edge*> edges) const;
	// Reorder edges and make it a loop. <canBeContinuousLoop> is asserted first.
	std::vector<Edge*> makeContinuousLoop(std::vector<Edge*> edges) const;

	Eigen::SparseMatrix<int> edge2vertexMap;
	Eigen::SparseMatrix<int> face2edgeMap;
	Eigen::SparseMatrix<int> cell2faceMap;

	// VertexGrid  --> (topology, geometry, index, type)
	XdmfGrid getXdmfVertexGrid() const;
	// EdgeGrid    --> (topology, geometry, index, type)
	XdmfGrid getXdmfEdgeGrid() const;
	// SurfaceGrid --> (topology, geometry, index, type)
	XdmfGrid getXdmfSurfaceGrid() const;
	// VolumeGrid  --> (topology, geometry, index, type)
	XdmfGrid getXdmfVolumeGrid() const;

	// Update vertexCoordinates matrix
	void update_vertexCoordinates();
	// Create e2v, f2e, c2f mapping
	void createIncidenceMaps();

	// Count the facets of different types
	void facetCounting();

private:
	std::string meshPrefix = "mesh";

	// check if AB form an edge
	bool isConnected(Vertex * A, Vertex * B) const;
	// check if faceEdges form a face 
	bool isConnected(const std::vector<Edge*> & faceEdges) const;
	// check if cellFaces form a cell
	bool isConnected(const std::vector<Face*> & cellFaces) const;

	// create incidence e2v matrix: A (-) B(+) 
	void create_edge2vertex_map();
	// create incidence f2e matrix
	void create_face2edge_map();
	// create incidence c2f matrix 
	void create_cell2face_map();

	// Append the coordinates of auxiliary vertices at the end of vertexCoodinates.
	// This matrix is intended for outputing XDMF file.
	const Eigen::MatrixXd getVertexCoordinatesExtended() const;

};

