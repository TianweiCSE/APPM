#include "MaxwellSolverImplicitEuler.h"



MaxwellSolverImplicitEuler::MaxwellSolverImplicitEuler()
{
}


MaxwellSolverImplicitEuler::MaxwellSolverImplicitEuler(const PrimalMesh * primal, const DualMesh * dual)
	: MaxwellSolver(primal, dual)
{
	const Eigen::VectorXi vertexTypes = primal->getVertexTypes();
	const int insulatingVertexType = static_cast<int>(Vertex::Type::Insulating);
	const int electrodeVertexType = static_cast<int>(Vertex::Type::Electrode);
	const int nVerticesTerminal = (vertexTypes.array() == electrodeVertexType).count();
	const int nVerticesBoundary = (vertexTypes.array() == insulatingVertexType).count() + nVerticesTerminal;

	const Eigen::VectorXi edgeTypes = primal->getEdgeTypes();
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);
	const int nEdges = primal->getNumberOfEdges();
	const int nEdgesInner = (edgeTypes.array() == interiorEdgeType).count();

	const int ndof =
		nEdgesInner +
		nVerticesBoundary;

	maxwellState = Eigen::VectorXd::Zero(ndof);
	// vectors of previous states
	x_m = maxwellState;
	x_mm1 = x_m;
}

MaxwellSolverImplicitEuler::MaxwellSolverImplicitEuler(const PrimalMesh * primal, const DualMesh * dual, MaxwellParams & param)
	: MaxwellSolverImplicitEuler(primal, dual)
{
	this->parameters = param;
}

MaxwellSolverImplicitEuler::~MaxwellSolverImplicitEuler()
{
}


void MaxwellSolverImplicitEuler::updateMaxwellState(const double dt, const double time)
{
	// update state vectors for previous timesteps
	x_mm1 = x_m;
	x_m = maxwellState;

	const double lambda2 = this->parameters.lambdaSquare;
	const double dt2 = pow(dt, 2);

	const int nTerminalVertices = primal->count_electrode_vertices();
	assert(nTerminalVertices > 0);

	// Setup of matrices for discrete operators
	Eigen::SparseMatrix<double> Q = setupOperatorQ().cast<double>();
	Eigen::SparseMatrix<double> Meps = setupOperatorMeps();
	Eigen::SparseMatrix<double> Mnu = setupOperatorMnu();
	Eigen::SparseMatrix<double> C = primal->get_f2eMap().cast<double>();

	// System of equations
	Eigen::SparseMatrix<double> temp_A1 = lambda2 * Q.transpose() * Meps * Q;
	Eigen::SparseMatrix<double> temp_A2 = Q.transpose() * C.transpose() * Mnu * C * Q;
	Eigen::SparseMatrix<double> A = temp_A1 + dt2 * temp_A2;
	A.makeCompressed();
	//Eigen::sparseMatrixToFile(A, "A.dat");
	
	const int nDirichlet = nTerminalVertices;
	const int nFree = maxwellState.size() - nDirichlet;
	//std::cout << "nFree: " << nFree << std::endl;

	Eigen::SparseMatrix<double> A_d = A.rightCols(nDirichlet);


	// Dirichlet conditions
	Eigen::VectorXd x_d(nDirichlet);
	x_d.setZero();
	const int nHalf = x_d.size() / 2;
	x_d.bottomRows(nHalf).array() = 1;   /// should be time-dependent voltage?
	//std::cout << "max(x_d): " << x_d.cwiseAbs().maxCoeff() << std::endl;
	
	Eigen::VectorXd src = Eigen::VectorXd::Zero(x_m.size());
	Eigen::VectorXd rhs = temp_A1 * (2 * x_m - x_mm1) + dt2 * src;
	rhs -= A_d * x_d;

	// Right hand side for system of equations
	Eigen::VectorXd rhs_free = rhs.topRows(nFree);
	//std::cout << "max(rhsFree): " << rhs_free.cwiseAbs().maxCoeff() << std::endl;

	// Solve for free DoF
	Eigen::VectorXd x_free(nFree);
	Eigen::SparseMatrix<double> A_free = A.topLeftCorner(nFree, nFree);
	A_free.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver(A_free);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}
	assert(maxwellSolver.info() == Eigen::Success);
	x_free = maxwellSolver.solve(rhs_free);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solve failed" << std::endl;
	}
	assert(maxwellSolver.info() == Eigen::Success);

	// Set new state vector
	maxwellState.topRows(nFree) = x_free;
	maxwellState.bottomRows(nDirichlet) = x_d.bottomRows(nDirichlet);

	// Set state vector for discrete data
	E_h = Q * maxwellState;
	B_h -= dt * C * E_h;

	//std::cout << "max(E): " << E_h.cwiseAbs().maxCoeff() << std::endl;

}
