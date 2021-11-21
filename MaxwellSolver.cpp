#include "MaxwellSolver.h"



MaxwellSolver::MaxwellSolver()
{
}

MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator)
: primal(primal), dual(dual), interpolator(interpolator)
{	
	//std::cout << "Dual mesh has " << dual->getNumberOfVertices() << " vertices" << std::endl;
	// init();
}


MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator, MaxwellParams & param) 
: MaxwellSolver(primal, dual, interpolator)
{
	parameters = param;
}

MaxwellSolver::~MaxwellSolver()
{

}

//void MaxwellSolver::updateMaxwellState(const double dt, const double time)
//{
//	std::cout << "You should call the inherited function, not this one" << std::endl;
//}

void MaxwellSolver::writeStates(H5Writer & writer) const
{
	const int nPrimalFaces = primal->getNumberOfFaces();
	const int nDualEdges = dual->getNumberOfEdges();
	const int nDualFaces = dual->getNumberOfFaces();

	writer.writeDoubleVector(B_h, "/bvec");
	writer.writeDoubleVector(E_h, "/evec");
	writer.writeDoubleVector(H_h, "/hvec");
	writer.writeDoubleVector(J_h, "/jvec");

	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd H(3, nDualEdges);
	Eigen::Matrix3Xd J(3, nDualFaces);

	for (int i = 0; i < nPrimalFaces; i++) {
		const Eigen::Vector3d fn = primal->getFace(i)->getNormal();
		const double fA = primal->getFace(i)->getArea();
		B.col(i) = (B_h(i) / fA) * fn;
	}
	for (int i = 0; i < nDualFaces; i++) {
		const Eigen::Vector3d fn = dual->getFace(i)->getNormal();
		const double fA = dual->getFace(i)->getArea();
		J.col(i) = J_h(i) / fA * fn;
	}
	writer.writeDoubleMatrix(B, "/B");
	writer.writeDoubleMatrix(J, "/J");

	const int nPrimalEdges = primal->getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primal->getEdge(i);
		E.col(i) = E_h(i) / edge->getLength() * edge->getDirection();
	}
	writer.writeDoubleMatrix(E, "/E");

	for (int i = 0; i < nDualEdges; i++) {
		const Edge * edge = dual->getEdge(i);
		H.col(i) = H_h(i) / edge->getLength() * edge->getDirection();
	}
	writer.writeDoubleMatrix(H, "/H");
}

const Eigen::VectorXd & MaxwellSolver::getBstate() const
{
	return B_h;
}

void MaxwellSolver::setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector)
{
	const int nFaces = primal->getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		const Face * face = primal->getFace(i);
		const double area = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		B_h(i) = area * fn.dot(fieldVector);
	}
}

void MaxwellSolver::setAzimuthalMagneticFluxField()
{
	const int nFaces = primal->getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		const Face * face = primal->getFace(i);
		const double fa = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		if (fn.cross(Eigen::Vector3d::UnitZ()).norm() > 0.99) {
			// azimuthal field
			Eigen::Vector3d fc = face->getCenter();
			fc(2) = 0;
			fc.normalize();
			B_h(i) = fn.dot(fc.cross(Eigen::Vector3d::UnitZ()));
		}
		else {
			B_h(i) = 0;
		}
	}
}

void MaxwellSolver::setTorusCurrent(const double x1, const double x2, const double z1, const double z2)
{
	std::cout << "Set torus current" << std::endl;
	const Eigen::Vector3d center = Eigen::Vector3d(0.0, 0.0, 0.5);

	const double tol = 8 * std::numeric_limits<double>::epsilon();

	const int nCells = dual->getNumberOfCells();

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dual->getCell(i);
		const Eigen::Vector3d center = cell->getCenter();
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const Eigen::MatrixXd cellVertexCoords = cell->getVertexCoordinates();
		const Eigen::ArrayXd xv = cellVertexCoords.col(0).array(); // x-coordinates of cell vertices
		const Eigen::ArrayXd yv = cellVertexCoords.col(1).array(); // 1-coordinates of cell vertices
		const Eigen::ArrayXd zv = cellVertexCoords.col(2).array(); // z-coordinates of cell vertices
		const bool isXZplane = (yv > 0).any() && (yv < 0).any();

		// Skip cells that are not across xz-plane
		if (!isXZplane || (xv < x1).all() || (xv > x2).all() || (zv < z1).all() || (zv > z2).all()) { continue; }

		// For each face of this cell
		for (auto face : cellFaces) {
			const int faceIdx = face->getIndex();
			const Eigen::Vector3d fn = face->getNormal();

			// Skip faces that have a normal vector out of xz plane
			if (abs(Eigen::Vector3d::UnitY().dot(fn)) > tol) { continue; }

			// Get vertex coordinates of this face
			const std::vector<Vertex*> faceVertices = face->getVertexList();
			Eigen::Matrix3Xd faceVertexPos(3, faceVertices.size());
			for (int i = 0; i < faceVertices.size(); i++) {
				faceVertexPos.col(i) = faceVertices[i]->getPosition();
			}
			const Eigen::ArrayXd x = faceVertexPos.row(0).array();
			const Eigen::ArrayXd z = faceVertexPos.row(2).array();
			const bool isXnormal = fn.cross(Eigen::Vector3d::UnitX()).norm() < tol;
			const bool isZnormal = fn.cross(Eigen::Vector3d::UnitZ()).norm() < tol;

			// Determine if face is pierced by one of the loop edges:
			// edge 1: z = z1, x = [x1,x2]
			bool isPiercedByEdge1 = isXnormal && (z <= (z1 + tol)).any() && (x >= x1).all() && (x <= x2).all();
			if (isPiercedByEdge1) {
				J_h(faceIdx) = fn.dot(Eigen::Vector3d::UnitX());
			}

			// edge 2: x = x1, z = [z1,z2]
			bool isPiercedByEdge2 = isZnormal && (x <= (x1 + tol)).any() && (z >= z1).all() && (z <= z2).all();
			if (isPiercedByEdge2) {
				J_h(faceIdx) = -fn.dot(Eigen::Vector3d::UnitZ());
			}

			// edge 3: z = z2, x = [x1,x2]
			bool isPiercedByEdge3 = isXnormal && (z >= (z2 - tol)).any() && (x >= x1).all() && (x <= x2).all();
			if (isPiercedByEdge3) {
				J_h(faceIdx) = -fn.dot(Eigen::Vector3d::UnitX());
			}

			// edge 4: x = x2, z = [z1,z2]
			bool isPiercedByEdge4 = isZnormal && (x >= (x2 - tol)).any() && (z >= z1).all() && (z <= z2).all();
			if (isPiercedByEdge4) {
				J_h(faceIdx) = fn.dot(Eigen::Vector3d::UnitZ());
			}
		}
	}
}

Eigen::VectorXd MaxwellSolver::electricPotentialTerminals(const double time)
{
	const int n = primal->count_electrode_vertices();
	assert(n % 2 == 0);
	const double t0 = 3;
	const double sigma_t = 1;
	const double phiA = 0;
	const double phiB = 0;
	const double phi1 = phiA * 0.5 * (1 + tanh((time - t0) / sigma_t));
	const double phi2 = phiB;
	Eigen::VectorXd phi_terminals(n);
	phi_terminals.topRows(n / 2).array() = phi1;
	phi_terminals.bottomRows(n / 2).array() = phi2;
	return phi_terminals;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_M_nu()
{
	if (M_nu.size() == 0) {
		const int N_A = primal->getNumberOfFaces();
		std::vector<T> triplets;
		triplets.reserve(N_A);
		for (int i = 0; i < N_A; i++) {
			const double fA = primal->getFace(i)->getArea();
			const double eL = dual->getEdge(i)->getLength();
			assert(dual->getEdge(i)->getVertexMid() == nullptr);
			const double value = eL / fA;
			triplets.emplace_back(T(i, i, value));
		}
		M_nu.resize(N_A, N_A);
		M_nu.setFromTriplets(triplets.begin(), triplets.end());
		M_nu.makeCompressed();
		assert(M_nu.size() > 0);
		assert(M_nu.nonZeros() > 0);
	}
	return M_nu;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_M_eps() {
	if (M_eps.size() == 0) {
		const int N_L = primal->getNumberOfEdges();
		std::vector<T> triplets;
		triplets.reserve(N_L);
		for (int i = 0; i < N_L; i++) {
			const double fA = dual->getFace(i)->getArea();
			const double eL = primal->getEdge(i)->getLength();
			const double value = fA / eL;
			triplets.emplace_back(T(i, i, value));
		}
		M_eps.resize(N_L, N_L);
		M_eps.setFromTriplets(triplets.begin(), triplets.end());
		M_eps.makeCompressed();
		assert(M_eps.size() > 0);
		assert(M_eps.nonZeros() > 0);
	}
	return M_eps;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_C_L_A() {
	if (C_L_A.size() == 0) {
		const int N_L = primal->getNumberOfEdges();
		const int N_A = primal->getNumberOfFaces();
		std::vector<T> triplets;
		triplets.reserve(4*N_A);  // In our case, each primal face would have no more than 4 edges. 
		for (int i = 0; i < N_A; i++) {
			Face* f = primal->getFace(i);
			for (Edge* e : f->getEdgeList()) {
				triplets.emplace_back(T(i, e->getIndex(), f->getOrientation(e)));
			}
		}
		C_L_A.resize(N_A, N_L);
		C_L_A.setFromTriplets(triplets.begin(), triplets.end());
		C_L_A.makeCompressed();
		//const Eigen::SparseMatrix<int>& p_f2e = primal->get_f2eMap();
		//Eigen::SparseMatrix<double> delta = (C_L_A - p_f2e).pruned();
		//assert(delta.nonZeros() == 0); 
	}
	return C_L_A;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_tC_Lo_Ao() {
	if (tC_Lo_Ao.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(6*dual->getNumberOfFaces()); // In our case, each dual face would have no more than 6 edges.
		for (Face* f : dual->getFaces()) {
			if (!f->isBoundary()) {
				for (Edge* e : f->getEdgeList()) {
					if (!e->isBoundary()) {
						triplets.emplace_back(T(f->getIndex(), e->getIndex(), f->getOrientation(e)));
					}
				}
			}
		}
		tC_Lo_Ao.resize(primal->getNumberOfEdges(), primal->getNumberOfFaces()); // Notice the amount matching
		tC_Lo_Ao.setFromTriplets(triplets.begin(), triplets.end());
		tC_Lo_Ao.makeCompressed();
		// const Eigen::SparseMatrix<double>&  d_f2e = dual->get_f2eMap();
		// Eigen::SparseMatrix<double> delta = 
		//		(tC_Lo_Ao - d_f2e.topLeftCorner(primal->getNumberOfEdges(), primal->getNumberOfFaces())).pruned();
		// assert(delta.nonZeros() == 0);
	}
	return tC_Lo_Ao; 
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_tC_pL_A() {
	if (tC_pL_A.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(3*dual->facet_counts.nE_boundary); // In our case, each dual boundary edge is shared by no more than 3 dual faces.
		for (Edge* e : dual->getEdges()) {
			if (e->isBoundary()) {
				for (Face* f : e->getFaceList()) {
					triplets.emplace_back(T(f->getIndex(), dpL2ph(e->getIndex()), f->getOrientation(e)));
				}
			}
		}
		tC_pL_A.resize(dual->getNumberOfFaces(), dual->facet_counts.nE_boundary);
		tC_pL_A.setFromTriplets(triplets.begin(), triplets.end());
		tC_pL_A.makeCompressed();
	}
	return tC_pL_A;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_pP_pL() {
	if (G_pP_pL.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(2*primal->facet_counts.nV_boundary);
		for (Edge* e : primal->getEdges()) {
			if (e->isBoundary()) {
				triplets.emplace_back(T(ppL2pe(e->getIndex()), ppP2phi(e->getVertexA()->getIndex()),  1.0));
				triplets.emplace_back(T(ppL2pe(e->getIndex()), ppP2phi(e->getVertexB()->getIndex()), -1.0));
			}
		}
		G_pP_pL.resize(primal->facet_counts.nE_boundary, primal->facet_counts.nV_boundary);
		G_pP_pL.setFromTriplets(triplets.begin(), triplets.end());
		G_pP_pL.makeCompressed();
	}
	return G_pP_pL;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_S_pP_pm() {
	if (S_pP_pm.size() == 0) {
		int count = 0;
		std::vector<T> triplets;
		triplets.reserve(primal->facet_counts.nV_electrode);
		for (Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Electrode) {
				if (v->getPosition()[2] > 0) {  // Anode
					triplets.emplace_back(T(count++, ppP2phi(v->getIndex()), 1.0));
				}
			}
		}
		for (Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Electrode) {
				if (v->getPosition()[2] < 0) {  // Cathode
					triplets.emplace_back(T(count++, ppP2phi(v->getIndex()), 1.0));
				}
			}
		}
		assert(count == primal->facet_counts.nV_electrode);
		S_pP_pm.resize(primal->facet_counts.nV_electrode, primal->facet_counts.nV_boundary);
		S_pP_pm.setFromTriplets(triplets.begin(), triplets.end());
		S_pP_pm.makeCompressed();
	}
	return S_pP_pm;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_tC_pL_AI() {
	if (tC_pL_AI.size() == 0) {
		int count = 0;
		std::vector<T> triplets;
		triplets.reserve(primal->facet_counts.nV_insulating * 6); // In our case, each dual face has no more than 6 edges
		for (Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Insulating) {
				Face* f = dual->getFace(dual->pVertex2dbFace(v->getIndex()));
				for (Edge* e : f->getEdgeList()) {
					assert(e->isBoundary());
					triplets.emplace_back(T(count, dpL2ph(e->getIndex()), f->getOrientation(e)));
				}
				count++;
			}
		}
		assert(count == primal->facet_counts.nV_insulating);
		tC_pL_AI.resize(count, dual->facet_counts.nE_boundary);
		tC_pL_AI.setFromTriplets(triplets.begin(), triplets.end());
		tC_pL_AI.makeCompressed();
	} 
	return tC_pL_AI;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_Q_LopP_L() {
	if (Q_LopP_L.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(primal->facet_counts.nE_interior + 2*primal->facet_counts.nV_boundary);
		
		for (int i = 0; i < primal->facet_counts.nE_interior; i++) {
			triplets.emplace_back(T(i, i, 1.0));
		}
		for (Edge* e : primal->getEdges()) {
			if (e->isBoundary()) {
				triplets.emplace_back(
					T(e->getIndex(), ppP2phi(e->getVertexA()->getIndex()) + primal->facet_counts.nE_interior,  1.0));
				triplets.emplace_back(
					T(e->getIndex(), ppP2phi(e->getVertexB()->getIndex()) + primal->facet_counts.nE_interior, -1.0));
			}
		}
		Q_LopP_L.resize(primal->getNumberOfEdges(), primal->facet_counts.nE_interior + primal->facet_counts.nV_boundary);
		Q_LopP_L.setFromTriplets(triplets.begin(), triplets.end());
		Q_LopP_L.makeCompressed();
	}
	return Q_LopP_L;
}

// TODO: Using dense matrix seems not efficient
Eigen::MatrixX3d&& MaxwellSolver::getInterpolated_E() const {
	Eigen::VectorXd temp(e.size() + dp.size());
	temp << e, dp; // concatenate [e, dp]^T
	return std::move(interpolator->get_E_interpolator().oneContract(temp));
}

Eigen::MatrixX3d&& MaxwellSolver::getInterpolated_B() const {
	Eigen::VectorXd temp(b.size() + hp.size());
	temp << b, hp;  // concatenate [b, hp]^T
	return std::move(interpolator->get_B_interpolator().oneContract(temp));
}

void MaxwellSolver::solveLinearSystem(const double time, 
									  const double dt, 
									  Eigen::SparseMatrix<double> &&M_sigma, 
									  Eigen::VectorXd &&j_aux) {
	// For the time being, linear system based on dense matrix is implemented 
	// for lack of built-in blocking functions for sparse matrix.  
	const int N_Lo    = primal->getNumberOfEdges() - primal->facet_counts.nE_boundary;
	const int N_pP    = primal->facet_counts.nV_boundary;
	const int N_pL    = primal->facet_counts.nE_boundary;
	const int tN_pA   = dual->facet_counts.nF_boundary;
	const int N_L     = primal->getNumberOfEdges();
	const int N_pP_pm = primal->facet_counts.nV_electrode;
	const int tN_AI   = primal->facet_counts.nV_insulating;
	assert(N_Lo + N_pP + N_pL + tN_pA == N_L + tN_pA + N_pP_pm + tN_AI);
	Eigen::MatrixXd mat(N_Lo + N_pP + N_pL + tN_pA, N_Lo + N_pP + N_pL + tN_pA);
	Eigen::VectorXd vec(N_Lo + N_pP + N_pL + tN_pA);
	const double lambda2 = parameters.lambdaSquare;

	// Assemble the square matrix
	mat.setZero();
	mat.block(0, 0, N_L, N_Lo + N_pP) = 
		(lambda2 / dt * get_M_eps() + get_tC_Lo_Ao() * get_M_nu() * dt * get_C_L_A() + M_sigma.block(0, 0, N_L, N_L))
		 * get_Q_LopP_L();
	mat.block(N_L, 0, tN_pA, N_Lo + N_pP) = M_sigma.block(N_L, 0, tN_pA, N_L) * get_Q_LopP_L();
	mat.block(0, N_Lo + N_pP, N_L + tN_pA, N_pL) = - get_tC_pL_A();
	mat.block(0, N_Lo + N_pP + N_pL, N_L + tN_pA, tN_pA) = M_sigma.block(0, N_L, N_L, tN_pA);
	mat.block(N_L, N_Lo + N_pP + N_pL, tN_pA, tN_pA) = 
		lambda2 / dt * Eigen::MatrixXd::Identity(tN_pA, tN_pA) + M_sigma.block(N_L, N_L, tN_pA, tN_pA);
	mat.block(N_L + tN_pA, N_Lo, N_pP_pm, N_pP_pm + N_pP) = get_S_pP_pm();
	mat.block(N_L + tN_pA + N_pP_pm, N_Lo + N_pP, tN_AI, N_pL) = get_tC_pL_AI();
	Eigen::SparseMatrix<double> sparse_mat = mat.sparseView();

	// Assemble the right vector
	vec.segment(0, N_L) = lambda2 / dt * get_M_eps() * e + get_tC_Lo_Ao() * get_M_nu() * b - j_aux.segment(0, N_L));
	vec.segment(N_L, tN_pA) = lambda2 / dt * dp - j_aux.segment(N_L, tN_pA);
	vec.segment(N_L + tN_pA, N_pP_pm / 2).setConstant(getPotential(time));
	vec.segment(N_L + tN_pA + N_pP_pm / 2, N_pP_pm / 2).setZero();
	vec.segment(N_L + tN_pA + N_pP_pm, tN_AI).setZero();

	// Solve
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(sparse_mat);
	if (solver.info() != Eigen::Success) {
		assert(false && "linear system not valid");
	}
	Eigen::VectorXd sol = solver.solve(vec);

	
}

void MaxwellSolver::timeStepping(const double dt) {

}







