#include "MaxwellSolverLeapFrog.h"



MaxwellSolverLeapFrog::MaxwellSolverLeapFrog()
{
}

MaxwellSolverLeapFrog::MaxwellSolverLeapFrog(const PrimalMesh * primal, const DualMesh * dual)
	: MaxwellSolver(primal, dual)
{
	const Eigen::VectorXi vertexType = primal->getVertexTypes();
	const Eigen::VectorXi edgeType = primal->getEdgeTypes();
	const Eigen::VectorXi faceType = primal->getFaceTypes();

	// Number of primal vertices
	const int nPv = primal->getNumberOfVertices();
	// Number of primal vertices on boundary
	const int nPvb = (vertexType.array() == static_cast<int>(Vertex::Type::Boundary)).count()
		+ (vertexType.array() == static_cast<int>(Vertex::Type::Terminal)).count();
	// Number of primal vertices in interior
	const int nPvi = (vertexType.array() == static_cast<int>(Vertex::Type::Inner)).count();
	// Number of primal vertrices at terminals
	
	this->nPvt = (vertexType.array() == static_cast<int>(Vertex::Type::Terminal)).count();
	assert((nPvb + nPvi) == nPv);

	// Number of primal edges
	const int nPe = primal->getNumberOfEdges();

	// Number of primal edges in interior
	const int nPei = (edgeType.array() == static_cast<int>(Edge::Type::Interior)).count()
		+ (edgeType.array() == static_cast<int>(Edge::Type::InteriorToBoundary)).count();

	// Number of primal faces
	const int nPf = primal->getNumberOfFaces();
	const int nPfi = (faceType.array() == 0).count();

	const int nDof = nPei + nPvb;
	x_m = Eigen::VectorXd::Zero(nDof);
	x_mm1 = Eigen::VectorXd::Zero(nDof);

	// Setup of operator Q = [id, G_b]
	Q = setupOperatorQ().cast<double>();

	Eigen::SparseMatrix<double> Meps;
	Eigen::SparseMatrix<double> Mnu;
	Meps = setupOperatorMeps();
	Mnu = setupOperatorMnu().topLeftCorner(nPfi, nPfi);

	// Curl operator on primal mesh
	Eigen::SparseMatrix<int> C;
	C = primal->get_f2eMap();
	this->C = C.cast<double>();

	// Curl operator for inner-faces and inner-edges
	Eigen::SparseMatrix<int> Cii = C.topLeftCorner(nPfi, nPei);

	// Setup of operator P = [Cii 0]
	Eigen::SparseMatrix<double> P;
	P = Eigen::SparseMatrix<double>(nPfi, nPei + nPvb);
	P.leftCols(Cii.cols()) = Cii.cast<double>();

	A = Q.transpose() * Meps * Q;
	B = P.transpose() * Mnu * P;
	assert(A.rows() == B.rows());
	assert(A.cols() == B.cols());


}


MaxwellSolverLeapFrog::~MaxwellSolverLeapFrog()
{
}

void MaxwellSolverLeapFrog::updateMaxwellState(const double dt, const double time)
{
	// System of equations
	assert(dt > 0);
	assert(std::isfinite(dt));
	this->M = A + pow(dt, 2) * B;
	const int nFreeIdx = M.rows() - nPvt;
	assert(nFreeIdx > 0);
	this->M_d = M.rightCols(nPvt);
	M_f = M.topLeftCorner(nFreeIdx, nFreeIdx);

	maxwellSolver.compute(M_f);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Solver initialization failed" << std::endl;
		exit(-1);
	}


	const int nDof = x_m.size();

	// Vector of degrees of freedom (dof)
	Eigen::VectorXd x(nDof);
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(x.size());

	rhs = A * (2 * x_m - x_mm1) + pow(dt, 2) * rhs;

	const int nPvt = M_d.cols();
	
	// Electric potential at terminal vertices
	Eigen::VectorXd phi_t(nPvt);
	double phi1 = 0;
	double phi2 = 0.5 * (1 + tanh((time - 5)));
	phi_t.topRows(nPvt / 2).array() = phi1;
	phi_t.bottomRows(nPvt / 2).array() = phi2;

	Eigen::VectorXd x_d;
	x_d = phi_t;
	rhs -= M_d * x_d;

	Eigen::VectorXd rhs_f = rhs.topRows(nFreeIdx);
	Eigen::VectorXd x_f(nFreeIdx);
	x_f = maxwellSolver.solve(rhs_f);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Solver solution failed" << std::endl;
		exit(-1);
	}
	x.topRows(x_f.size()) = x_f;
	x.bottomRows(x_d.size()) = x_d;

	// Electric field (from electrostatic boundary conditions)
	E_h = Q * x;

	// Update of magnetic field
	B_h = B_h - dt * C * E_h;

	// Update state vectors
	x_mm1 = x_m;
	x_m = x;




}
