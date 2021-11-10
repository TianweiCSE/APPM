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

	/// Dual edges and faces follow the orientation of the associated primal faces and edges.
	void init_dualMesh();
	void check() const;

private:
	PrimalMesh* primal;

	using Mesh::addEdge;
	/// This function is for adding Edge2
	Edge* addEdge(Edge* e1, Edge* e2);

	/// Cell center < 1 --> FLUID; Cell center > 1 --> SOLID 
	void init_cellFluidType();
	/// The terminal sides of plasma --> OPENING; The interface of solid and fluid --> WALL; 
	/// Inside fluid --> INTERIOR; Else --> DEFAULT
	void init_faceFluidType();

};

