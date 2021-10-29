#pragma once
#include "MaxwellSolver.h"
class MaxwellSolverCrankNicholson :
	public MaxwellSolver
{
public:
	MaxwellSolverCrankNicholson();
	MaxwellSolverCrankNicholson(const PrimalMesh * primal, const DualMesh * dual);

	~MaxwellSolverCrankNicholson();

	void updateMaxwellState(const double dt, const double time) override;
};

