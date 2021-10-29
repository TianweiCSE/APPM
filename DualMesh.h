#pragma once

#include "EigenAuxiliaries.h"
#include "Mesh.h"
#include "PrimalMesh.h"

class DualMesh :
	public Mesh
{
public:
	DualMesh();
	DualMesh(const std::string & meshPrefix);
	~DualMesh();

	/// Dual edges and faces follow the orientation of the associated primal faces and edges.
	void init_dualMesh(const PrimalMesh & primal);

protected:
	/// SurfaceGrid --> attribute(topology, geometry, index, area, type)
	XdmfGrid getXdmfSurfaceGrid() const;
	/// VolumeGrid -->  attribute(topology, geometry, index, volume, type)
	XdmfGrid getXdmfVolumeGrid() const;


private:
	/// Cell center < 1 --> FLUID; Cell center > 1 --> SOLID 
	void init_cellFluidType();
	/// The terminal sides of plasma --> OPENING; The interface of solid and fluid --> WALL; 
	/// Inside fluid --> INTERIOR; Else --> DEFAULT
	void init_faceFluidType();
};

