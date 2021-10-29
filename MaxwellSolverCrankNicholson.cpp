#include "MaxwellSolverCrankNicholson.h"



MaxwellSolverCrankNicholson::MaxwellSolverCrankNicholson()
{
}

MaxwellSolverCrankNicholson::MaxwellSolverCrankNicholson(const PrimalMesh * primal, const DualMesh * dual)
	: MaxwellSolver(primal, dual)
{
	const Eigen::VectorXi vertexTypes = primal->getVertexTypes();
	const int boundaryVertexType = static_cast<int>(Vertex::Type::Boundary);
	const int terminalVertexType = static_cast<int>(Vertex::Type::Terminal);
	const int nVerticesTerminal = (vertexTypes.array() == terminalVertexType).count();
	const int nVerticesBoundary = (vertexTypes.array() == boundaryVertexType).count() + nVerticesTerminal;

	const Eigen::VectorXi edgeTypes = primal->getEdgeTypes();
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);
	const int interiorToBoundaryEdgeType = static_cast<int>(Edge::Type::InteriorToBoundary);
	const int boundaryEdgeType = static_cast<int>(Edge::Type::Boundary);
	const int nEdges = primal->getNumberOfEdges();
	const int nEdgesInner =
		(edgeTypes.array() == interiorEdgeType).count() +
		(edgeTypes.array() == interiorToBoundaryEdgeType).count();

	const Eigen::VectorXi faceTypes = primal->getFaceTypes();
	const int isBoundaryFace = 1;
	const int nFaces = primal->getNumberOfFaces();
	const int nFacesInner = (faceTypes.array() != isBoundaryFace).count();

	const int ndof =
		nFacesInner + 
		nEdgesInner +
		nVerticesBoundary;
	maxwellState = Eigen::VectorXd::Zero(ndof);

}


MaxwellSolverCrankNicholson::~MaxwellSolverCrankNicholson()
{
}

void MaxwellSolverCrankNicholson::updateMaxwellState(const double dt, const double time)
{
	const Eigen::VectorXd prevState = maxwellState;

	// TODO move variables to initialization that are repeatedly used
	const int n = maxwellState.size();
	const int nFacesInterior = primalMeshInfo.nFacesInner;
	const int nEdgesInterior = primalMeshInfo.nEdgesInner;
	const int nEdges = primalMeshInfo.nEdges;
	const int nVerticesTerminal = primalMeshInfo.nVerticesTerminal;
	const int nVerticesBoundary = primalMeshInfo.nVerticesBoundary;
	const int nVertices = primalMeshInfo.nVertices;

	// Discrete Operators
	const Eigen::SparseMatrix<double> M_mu = get_Mnu();
	const Eigen::SparseMatrix<double> M_eps = get_Meps();
	const Eigen::SparseMatrix<double> M_sigma = get_Msigma();

	const Eigen::SparseMatrix<double> G = primal->get_e2vMap().cast<double>();
	const Eigen::SparseMatrix<double> C = primal->get_f2eMap().cast<double>();

	const Eigen::SparseMatrix<double> X = get_bndryInclOp();

	// Boundary gradient operator
	const Eigen::SparseMatrix<double> GX = G * X;

	// Operator from inner edges and boundary vertices to all edges
	const int nRowsQ = nEdgesInterior + nVerticesBoundary;
	const int nQ = std::max(nEdges, nRowsQ);
	Eigen::SparseMatrix<double> id(nQ, nQ);
	id.setIdentity();
	Eigen::SparseMatrix<double> Q = id.topLeftCorner(nEdges, nRowsQ); 
	Q.rightCols(GX.cols()) = -GX;

	// Curl operator in interior (inner faces, inner edges)
	const Eigen::SparseMatrix<double> C_inner = C.topLeftCorner(nFacesInterior, nEdgesInterior);
	assert(C_inner.rows() == nFacesInterior);
	assert(C_inner.cols() == nEdgesInterior);
	Eigen::SparseMatrix<double> P(nFacesInterior, nEdgesInterior + nVerticesBoundary);
	P.leftCols(C_inner.cols()) = C_inner;

	// Solve system of equations: A*u(k+1) = B*u(k) for u(k+1)
	Eigen::SparseMatrix<double> A(n, n);
	Eigen::SparseMatrix<double> B(n, n);
	Eigen::SparseMatrix<double> A_11, A_12, A_21, A_22;
	A_11 = M_mu;
	A_12 = +0.5 * dt * P;
	A_21 = -0.5 * dt * P.transpose();
	if (isMaxwellCurrentSource) {
		A_22 = Q.transpose() * M_eps * Q;
	}
	else {
		A_22 = Q.transpose() * (M_eps + 0.5 * dt * M_sigma) * Q;
	}

	const Eigen::SparseMatrix<double> Aleft = Eigen::vertcat(A_11, A_21);
	const Eigen::SparseMatrix<double> Aright = Eigen::vertcat(A_12, A_22);
	A.leftCols(Aleft.cols()) = Aleft;
	A.rightCols(Aright.cols()) = Aright;

	Eigen::SparseMatrix<double> B_11, B_12, B_21, B_22;
	B_11 = M_mu;
	B_12 = -0.5 * dt * P;
	B_21 = +0.5 * dt * P.transpose();
	if (isMaxwellCurrentSource) {
		B_22 = Q.transpose() * M_eps * Q;
	}
	else {
		B_22 = Q.transpose() * (M_eps - 0.5 * dt * M_sigma) * Q;
	}

	const Eigen::SparseMatrix<double> Bleft = Eigen::vertcat(B_11, B_21);
	const Eigen::SparseMatrix<double> Bright = Eigen::vertcat(B_12, B_22);
	B.leftCols(Bleft.cols()) = Bleft;
	B.rightCols(Bright.cols()) = Bright;

	int nDirichlet = 0;
	if (isMaxwellCurrentSource) {
		nDirichlet = nVerticesTerminal / 2; // number of fixed values by boundary condition
	}
	else {
		nDirichlet = nVerticesTerminal; // number of fixed values by boundary condition
	}
	const int nFree = n - nDirichlet; // number of free unknowns 

	const Eigen::SparseMatrix<double> M_f = A.topLeftCorner(nFree, nFree);
	const Eigen::SparseMatrix<double> M_d = A.rightCols(nDirichlet);

	// Solve system of equations for free unknowns
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver(M_f);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver failed to setup the system of equations" << std::endl;
	}
	const Eigen::VectorXd x_d = electricPotentialTerminals(time); // voltage boundary condition on terminals at new timestep
	Eigen::VectorXd rhs;  // right hand side of equation

	if (isMaxwellCurrentSource) {
		assert(M_d.cols() == nDirichlet);
		rhs = B * prevState;
		rhs -= M_d * x_d.bottomRows(nDirichlet);
	}
	else {
		rhs = B * prevState - M_d * x_d;
	}

	const int nx = Q.cols();
	std::cout << "Q.size: " << Q.rows() << " x " << Q.cols() << std::endl;
	std::cout << "nFacesInner: " << primalMeshInfo.nEdges << std::endl;
	if (isMaxwellCurrentSource) {
		const Eigen::VectorXd temp = -1 * dt * Q.transpose() * J_h.topRows(primalMeshInfo.nEdges); // TODO this should be ... + 1/2 * dt * Q^t * (J^(k+1) + J(k))
		rhs.bottomRows(nx) += temp;
	}
	else {
		rhs = B * prevState - M_d * x_d;
	}


	const Eigen::VectorXd rhs_f = rhs.topRows(M_f.rows());          // discard equations to terminal vertices
	Eigen::VectorXd x_f(nFree);
	x_f = maxwellSolver.solve(rhs_f);       // solve system of equations
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver failed to solve for vector of unknowns" << std::endl;
	}

	// assign solution to state vector
	maxwellState.topRows(nFree) = x_f;
	maxwellState.bottomRows(nDirichlet) = x_d.bottomRows(nDirichlet);

	// map data to vectors
	const int nH = primalMeshInfo.nFacesInner;
	H_h.topRows(nH) = maxwellState.topRows(nH);
	E_h = Q * maxwellState.bottomRows(nx);

	B_h.topRows(M_mu.rows()) = M_mu * H_h.topRows(M_mu.cols());
}
