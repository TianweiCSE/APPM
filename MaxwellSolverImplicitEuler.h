#pragma once
#include "MaxwellSolver.h"
class MaxwellSolverImplicitEuler :
	public MaxwellSolver
{
public:
	MaxwellSolverImplicitEuler();
	MaxwellSolverImplicitEuler(const PrimalMesh * primal, const DualMesh * dual);
	MaxwellSolverImplicitEuler(const PrimalMesh * primal, const DualMesh * dual, MaxwellParams & param);

	~MaxwellSolverImplicitEuler();

	void updateMaxwellState(const double dt, const double time) override;

private:
	Eigen::VectorXd x_m, x_mm1;

};

