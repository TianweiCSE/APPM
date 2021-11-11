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

	// Dual edges and faces follow the orientation of the associated primal faces and edges.
	void init_dualMesh();
	void check() const;

	// primal boundary face index ---> dual boundary vertex index 
	const int pFace2dbVertex(const int idx) {return primalFaceToDualBoundaryVertex.coeffRef(idx);};
	// primal boundary edge index ---> dual boundary edge index 
	const int pEdge2dbEdge  (const int idx) {return primalEdgeToDualBoundaryEdge.coeffRef(idx);};
	// primal boundary vertex index ---> dual boundary face index
	const int pVertex2dbFace(const int idx) {return primalVertexToDualBoundaryFace.coeffRef(idx);};

private:
	PrimalMesh* primal;

	Eigen::SparseVector<int> primalFaceToDualBoundaryVertex;
	Eigen::SparseVector<int> primalEdgeToDualBoundaryEdge;
	Eigen::SparseVector<int> primalVertexToDualBoundaryFace;

	using Mesh::addEdge;
	/// This function is for adding Edge2
	Edge* addEdge(Edge* e1, Edge* e2);

	/// Cell center < 1 --> FLUID; Cell center > 1 --> SOLID 
	void init_cellFluidType();
	/// The terminal sides of plasma --> OPENING; The interface of solid and fluid --> WALL; 
	/// Inside fluid --> INTERIOR; Else --> DEFAULT
	void init_faceFluidType();

};

