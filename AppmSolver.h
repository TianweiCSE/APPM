#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "Numerics.h"
#include "Interpolator.h"
#include "FluidSolver.h"
#include "TwoFluidSolver.h"
#include "MaxwellSolver.h"

#include <Eigen/SparseLU>

/**
 * @brief A class that wraps Maxwell solver and fluid solver. 
 * 
 */
class AppmSolver
{

public:
	AppmSolver();
	AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams);
	~AppmSolver();

	void run();

protected:
	PrimalMesh* primalMesh;
	DualMesh*   dualMesh;

	TwoFluidSolver* twofluidSolver = nullptr;
	MaxwellSolver* maxwellSolver = nullptr;

	Interpolator* interpolator;

private:
	int maxIterations = 100;
	double maxTime = 1;
	double lambda = 1.0;
	double alpha = 6.4388;  // 643.88 in 1d case
	int itersPerWrite = 1;
	bool with_friction = true;

	int iteration;
	double time;
	std::vector<std::pair<int,double>> timeStamps; //< store the (iteration, time) at which the snapshot is recorded.

	void init_meshes(const PrimalMesh::PrimalMeshParams & primalParams);

	// Collect the transient solutions on primal vertices into "solutions_primal_vertex.xdmf"
	void writeSolutionPrimalVertex() const;
	// Collect the transient solutions on primal edges into "solutions_primal_edge.xdmf"
	void writeSolutionPrimalEdge() const;
	// Collect the transient solutions on primal faces into "solutions_primal_face.xdmf"
	void writeSolutionPrimalFace() const;
	// Collect the transient solutions on dual cells into "solutions_dual_cell.xdmf"
	void writeSolutionDualCell() const;
	// Collect the transient solutions on dual cells into "solutions_dual_face.xdmf"
	void writeSolutionDualFace() const;
	
	// Compute the norm of solutions
	void writeSolutionNorms() const;

	XdmfGrid getSnapshotPrimalVertex(const int iteration) const;
	XdmfGrid getSnapshotPrimalEdge  (const int iteration) const;
	XdmfGrid getSnapshotPrimalFace  (const int iteration) const;
	XdmfGrid getSnapshotDualEdge    (const int iteration) const;
	XdmfGrid getSnapshotDualFace    (const int iteration) const;
	XdmfGrid getSnapshotDualCell    (const int iteration) const;

	void applyInitialConditions();
	void applyInitialConditions(const std::string h5_file, const double t);

	/**
	 * @brief Write snapshot of solutions to "appm-<iteration>.h5"
	 * 
	 * This function is intended to be called at each iteration to store the transient solutions.
	 * At last by calling writeSolution<...> functions, the snapshots will be collected into .xdmf file
	 * which can be read by Paraview in XdmfReaderS mode.
	 * 
	 * @param iteration iteration number
	 * @param time present time
	 */
	void writeSnapshot(const int iteration, const double time);

	void verboseDiagnosis() const;

	void readParameters(const std::string filename);
};

