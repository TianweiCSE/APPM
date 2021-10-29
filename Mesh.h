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
	struct MeshInfo {
		int nVertices = 0;         // number of vertices
		int nVerticesBoundary = 0; // number of vertices on domain boundary
		int nVerticesTerminal = 0; // number of degrees of freedom with Dirichlet conditions

		int nEdges = 0;      // number of edges
		int nEdgesInner = 0; // number of edges in interior of domain

		int nFaces = 0;      // number of faces
		int nFacesInner = 0; // number of faces in interior of domain

		int nCells = 0; // number of cells
	};

	Mesh();
	Mesh(const std::string & meshPrefix);
	~Mesh();

	/// - write incidence map e2v, v2f, f2c, f2v to .dat files
	/// - write index v.s. type/length/center/normal/indices of subfacets to .h5 file   
	void writeToFile();

	/// prefix-mesh.xdmf   << root --> domain --> treeGrid --> (vertexGrid, edgeGrid, surfaceGrid)
	/// prefix-volume.xdmf << root --> domain --> VolumeGrid
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
	Face *   getFace(const int index) const;
	Face *   getFace(const std::vector<Edge*> & faceEdges);
	Cell *   getCell(const int index) const;
	Cell *   getCell(const std::vector<Face*> & cellFaces);

	const std::string getPrefix() const;
	
	void createIncidenceMaps();

	const int getNumberOfVertices() const;
	const int getNumberOfEdges() const;
	const int getNumberOfFaces() const;
	const int getNumberOfCells() const;

	/// return xdmf vector containing [16 #faces #vertices idx1 idx2 idx3 .. #vertices idx1 idx2 idx3 .. ]
	const std::vector<int> getXdmfTopology_cell2vertexIndices() const;
	/// return xdmf vector containing [facetype #vertices idx1 idx2 idx3 .. ]
	const std::vector<int> getXdmfTopology_face2vertexIndices() const;
	/// return cellList
	const std::vector<Cell*> getCells() const;
	/// return faceList
	const std::vector<Face*> getFaces() const;

	void check() const;

	const Eigen::SparseMatrix<int> & get_f2eMap() const;
	const Eigen::SparseMatrix<int> & get_e2vMap() const;

	const Eigen::VectorXi getVertexTypes() const;
	const Eigen::VectorXi getEdgeTypes() const;
	const Eigen::VectorXi getFaceTypes() const;

	MeshInfo getMeshInfo() const;


protected:
	std::vector<Vertex*> vertexList;
	std::vector<Edge*> edgeList;
	std::vector<Face*> faceList;
	std::vector<Cell*> cellList;
	Eigen::Matrix3Xd vertexCoordinates;

	/// reorder edges and make it a loop
	std::vector<Edge*> makeContinuousLoop(std::vector<Edge*> edges);

	Eigen::SparseMatrix<int> edge2vertexMap;
	Eigen::SparseMatrix<int> face2edgeMap;
	Eigen::SparseMatrix<int> cell2faceMap;

	/// SurfaceGrid --> (topology, geometry, index, area)
	virtual XdmfGrid getXdmfSurfaceGrid() const;
	/// VolumeGrid --> (topology, geometry, index, volume)
	virtual XdmfGrid getXdmfVolumeGrid() const;


private:
	std::string meshPrefix = "mesh";


	/// check if AB form an edge
	bool isConnected(const Vertex * A, const Vertex * B) const;
	/// check if faceEdges form a face 
	bool isConnected(const std::vector<Edge*> & faceEdges) const;
	/// check if cellFaces form a cell
	bool isConnected(const std::vector<Face*> & cellFaces) const;
	Edge * getEdge(Vertex * A, Vertex * B);

	/// useless
	void create_vertexCoordinates();

	/// create incidence e2v matrix: A (-) B(+) 
	void create_edge2vertex_map();
	/// create incidence f2e matrix
	void create_face2edge_map();
	/// create incidence c2f matrix 
	void create_cell2face_map();

	/// VertexGrid --> (topology, geometry, index, type)
	XdmfGrid getXdmfVertexGrid() const;
	/// EdgeGrid --> (topology, geometry, index, type)
	XdmfGrid getXdmfEdgeGrid() const;
	
	/// prefix-volume.xdmf << root --> domain --> VolumeGrid
	void writeXdmfVolumeMesh() const;

	//void writeXdmf_surface();
	//void writeXdmf_volume();

};

