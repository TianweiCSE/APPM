#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "Numerics.h"
#include "FluidSolver.h"
#include "SingleFluidSolver.h"
#include "TwoFluidSolver.h"
#include "MultiFluidSolver.h"
#include "MaxwellSolver.h"
#include "MaxwellSolverImplicitEuler.h"

#include <Eigen/SparseLU>


#define _RT_ONECELL 
#undef _RT_ONECELL


class AppmSolver
{

public:
	AppmSolver();
	/// Init primal mesh, dual mesh, MaxwellSolver, FluidSolver, RT_Interpolation
	AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams);
	~AppmSolver();

	/// updateFluidState --> updateMaxwellState --> interpolateMagneticFluxToPrimalVertices
	void run();

protected:
	PrimalMesh* primalMesh;
	DualMesh* dualMesh;

	FluidSolver * fluidSolver = nullptr;
	MaxwellSolver * maxwellSolver = nullptr;

	Eigen::Matrix3Xd B_vertex;  /// The interpolant of B-field on each vertex

	bool isMaxwellCurrentSource = false;

private:
	bool isMaxwellEnabled = false;
	bool isFluidEnabled = true;
	double timestepSize = 1.0;
	int maxIterations = 0;
	double maxTime = 0;
	double lambdaSquare = 1.0;


	/// For each primal triangle prism, the 1st order face element is used to reconstruct the B-field in the prism.
	/// For each primal vertex, the contributions from each face that it belongs to are summed up and averaged. 
	void interpolateMagneticFluxToPrimalVertices();

	// void test_raviartThomas();
	const Eigen::Matrix3Xd getPrismReferenceCoords(const int nSamples);

	std::vector<double> timeStamps;

	void init_meshes(const PrimalMesh::PrimalMeshParams & primalParams);

	void writeXdmf();
	void writeXdmfDualVolume();

	void writeOutput(const int iteration, const double time);

	XdmfGrid getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	
	XdmfGrid getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualVolumeGrid(const int iteration, const double time, const std::string & dataFilename);

	/// Assign Piola transformation to each triangle prism. (For interpolating B-field)
	void init_RaviartThomasInterpolation();

	std::vector<Eigen::Matrix3d> rt_piolaMatrix;
	std::vector<Eigen::Vector3d> rt_piolaVector;

	void readParameters(const std::string & filename);

	const std::string xdmf_GridPrimalEdges(const int iteration) const;
	const std::string xdmf_GridPrimalFaces(const int iteration) const;
	const std::string xdmf_GridDualEdges(const int iteration) const;
	const std::string xdmf_GridDualFaces(const int iteration) const;
	const std::string xdmf_GridDualCells(const int iteration) const;
};

