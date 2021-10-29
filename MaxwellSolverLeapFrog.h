#pragma once
#include "MaxwellSolver.h"
class MaxwellSolverLeapFrog :
	public MaxwellSolver
{
public:
	MaxwellSolverLeapFrog();
	MaxwellSolverLeapFrog(const PrimalMesh * primal, const DualMesh * dual);
	~MaxwellSolverLeapFrog();

	void updateMaxwellState(const double dt, const double time);
	
private:
	int nPvt = 0;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseMatrix<double> B;
	Eigen::SparseMatrix<double> C;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> M_d, M_f;
	Eigen::SparseMatrix<double> Q;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver;
	Eigen::VectorXd x_m, x_mm1;

};

