#pragma once

#include "EigenAuxiliaries.h"
#include "Mesh.h"
#include "PrimalMesh.h"

class DualMesh :
	public Mesh
{
public:
	DualMesh(PrimalMesh* primal);
	DualMesh(const std::string & meshPrefix, PrimalMesh* primal);
	~DualMesh();

	const std::vector<Cell*> getFluidCells() const; 
	const std::vector<Face*> getFluidFaces() const;

	// Dual edges and faces follow the orientation of the associated primal faces and edges.
	void init();
	void check() const;
	void outputFaceInfo(const int idx) const; 

	// primal boundary face index ---> dual boundary vertex index 
	int pFace2dbVertex(const int idx) const {return primalFaceToDualBoundaryVertex.coeff(idx);};
	// primal boundary edge index ---> dual boundary edge index 
	int pEdge2dbEdge  (const int idx) const {return primalEdgeToDualBoundaryEdge.coeff(idx);};
	// primal boundary vertex index ---> dual boundary face index
	int pVertex2dbFace(const int idx) const {return primalVertexToDualBoundaryFace.coeff(idx);};

private:
	PrimalMesh* primal;
	const double fluidRadius = 1.0001; // Any cell with radial distance smaller than fluidRadius is considered as fluid cell

	mutable std::vector<Cell*> fluidCellList; 
	mutable std::vector<Face*> fluidFaceList; 

	Eigen::SparseVector<int> primalFaceToDualBoundaryVertex;
	Eigen::SparseVector<int> primalEdgeToDualBoundaryEdge;
	Eigen::SparseVector<int> primalVertexToDualBoundaryFace;

	using Mesh::addEdge;
	/// This function is for adding Edge2
	Edge* addEdge(Edge* e1, Edge* e2);

	/// Cell center < fluidRadius --> Fluid; Cell center > fluidRadius --> Solid 
	void init_cellFluidType();
	/**
	 * @brief Assign fluid type to each dual face
	 * 		- Inside fluid --> Interior
	 * 		- At terminal sides --> Opening
	 * 		- Interface of solid and fluid domain (or the transverse boundary of the whole domain) --> Wall
	 * 		- Multiple types in one face --> Mixed
	 * 		- else (interior faces in the solid domain) --> Undefined 
	 */
	void init_faceFluidType();

	enum class ElectrodeGeometry{
		Round, Square
	} electrodeGeo;

};

