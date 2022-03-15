#include "MaxwellSolver.h"



MaxwellSolver::MaxwellSolver()
{
}

MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator)
: primal(primal), dual(dual), interpolator(interpolator)
{	
	eo.resize(primal->facet_counts.nE_interior); eo.setZero();
	phi.resize(primal->facet_counts.nV_boundary); phi.setZero();
	hp.resize(primal->facet_counts.nE_boundary); hp.setZero();
	dp.resize(dual->facet_counts.nF_boundary); dp.setZero();
	e.resize(primal->getNumberOfEdges()); e.setZero();
	b.resize(primal->getNumberOfFaces()); b.setZero();
	j.resize(dual->getNumberOfFaces()); j.setZero();
}


MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator, MaxwellParams & param) 
: MaxwellSolver(primal, dual, interpolator)
{
	parameters = param;
}

MaxwellSolver::~MaxwellSolver()
{

}

void MaxwellSolver::applyInitialCondition() {
	e.setZero();
	b.setZero();
	phi.setZero();
	hp.setZero();
	dp.setZero();
	updateInterpolated_E();
	updateInterpolated_B();
}

void MaxwellSolver::applyInitialCondition(const std::string h5_file) {
	H5Reader reader(h5_file);
	e = reader.readVectorData("/e");
	b = reader.readVectorData("/b");
	phi = reader.readVectorData("/phi");
	hp = reader.readVectorData("/hp");
	dp = reader.readVectorData("/dp");
	updateInterpolated_E();
	updateInterpolated_B();
}

void MaxwellSolver::writeSnapshot(H5Writer & writer) const
{
	const int nPrimalFaces = primal->getNumberOfFaces();
	const int nDualEdges = dual->getNumberOfEdges();
	const int nDualFaces = dual->getNumberOfFaces();

	writer.writeDoubleVector(e, "/e");
	writer.writeDoubleVector(b, "/b");
	writer.writeDoubleVector(j, "/j");
	writer.writeDoubleVector(phi, "/phi");
	writer.writeDoubleVector(hp, "/hp");
	writer.writeDoubleVector(dp, "/dp");
	writer.writeDoubleMatrix(getInterpolated_E(), "/E");
	writer.writeDoubleMatrix(getInterpolated_B(), "/B");

	Eigen::VectorXd phi_extended(primal->getNumberOfVertices());
	phi_extended.setConstant(-1.0);
	for (int i = 0; i < primal->getNumberOfVertices(); i++) {
		if (primal->getVertex(i)->isBoundary()) {
			phi_extended[i] = phi[ppP2phi(i)];
		}
	}
	Eigen::VectorXd h_extended(dual->getNumberOfEdges());
	h_extended.setZero();
	for (int i = 0; i < dual->getNumberOfEdges(); i++) {
		if (dual->getEdge(i)->isBoundary()) {
			h_extended[i] = hp[dpL2ph(i)];
		}
	}
	writer.writeDoubleVector(phi_extended, "/phi_extended");
	writer.writeDoubleVector(h_extended, "/h_extended");
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
			const Face* f = primal->getFace(i);
			for (const Edge* e : f->getEdgeList()) {
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
		for (const Face* f : dual->getFaces()) {
			if (!f->isBoundary()) {
				for (const Edge* e : f->getEdgeList()) {
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
		for (const Edge* e : dual->getEdges()) {
			if (e->isBoundary()) {
				for (const Face* f : e->getFaceList()) {
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
		for (const Edge* e : primal->getEdges()) {
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
		for (const Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Electrode) {
				if (v->getPosition()[2] > 0) {  // Anode
					triplets.emplace_back(T(count++, ppP2phi(v->getIndex()), 1.0));
				}
			}
		}
		for (const Vertex* v : primal->getVertices()) {
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
		for (const Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Insulating) {
				const Face* f = dual->getFace(dual->pVertex2dbFace(v->getIndex()));
				for (const Edge* e : f->getEdgeList()) {
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
			assert(!primal->getEdge(i)->isBoundary());
			triplets.emplace_back(T(i, i, 1.0));
		}
		for (const Edge* e : primal->getEdges()) {
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

const Eigen::SparseMatrix<double>& MaxwellSolver::get_tS_pA_AI() {
	if (tS_pA_AI.size() == 0) {
		int count = 0;
		std::vector<T> triplets;
		triplets.reserve(primal->facet_counts.nV_insulating);
		for (const Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Insulating) {
				const Face* f = dual->getFace(dual->pVertex2dbFace(v->getIndex()));
				triplets.emplace_back(T(count, dpA2pd(f->getIndex()), 1.0));
				count++;
			}
		}
		assert(count == primal->facet_counts.nV_insulating);
		tS_pA_AI.resize(count, dual->facet_counts.nF_boundary);
		tS_pA_AI.setFromTriplets(triplets.begin(), triplets.end());
		tS_pA_AI.makeCompressed();
	}
	return tS_pA_AI;
}

Eigen::MatrixXd MaxwellSolver::updateInterpolated_E() {
	Eigen::VectorXd temp(e.size() + dp.size());
	temp << e, dp; // concatenate [e, dp]^T
	E = interpolator->get_E_interpolator().oneContract(temp);
	return E;
}

Eigen::MatrixXd MaxwellSolver::updateInterpolated_B() {
	Eigen::VectorXd temp(b.size() + hp.size());
	temp << b, hp;  // concatenate [b, hp]^T
	B = interpolator->get_B_interpolator().oneContract(temp);
	return B;
}

Eigen::MatrixXd MaxwellSolver::getInterpolated_E() const {
	return E;
}

Eigen::MatrixXd MaxwellSolver::getInterpolated_B() const {
	return B;
}

void MaxwellSolver::solveLinearSystem(const double time, 
									  const double dt, 
									  Eigen::SparseMatrix<double>&& M_sigma, 
									  Eigen::VectorXd&& j_aux) {
	// For the time being, linear system based on dense matrix is implemented 
	// for lack of built-in blocking functions for sparse matrix.  
	std::cout << "-- Start assembling linear system ..." << std::endl;
	const int N_Lo    = primal->getNumberOfEdges() - primal->facet_counts.nE_boundary;
	const int N_pP    = primal->facet_counts.nV_boundary;
	const int N_pL    = primal->facet_counts.nE_boundary;
	const int tN_pA   = dual->facet_counts.nF_boundary;
	const int N_L     = primal->getNumberOfEdges();
	const int N_pP_pm = primal->facet_counts.nV_electrode;
	const int tN_AI   = primal->facet_counts.nV_insulating;
	assert(N_Lo + N_pP + N_pL + tN_pA == N_L + tN_pA + N_pP_pm + tN_AI);
	Eigen::SparseMatrix<double, Eigen::ColMajor> mat(N_Lo + N_pP + N_pL + tN_pA, N_Lo + N_pP + N_pL + tN_pA);
	std::vector<Eigen::Triplet<double>> triplets;
	Eigen::VectorXd vec(N_Lo + N_pP + N_pL + tN_pA);
	vec.setZero();
	const double lambda2 = parameters.lambdaSquare;

	// Assemble the square matrix
	Eigen::blockFill<double>(triplets, 0, 0, 
		((lambda2 / dt * get_M_eps() + get_tC_Lo_Ao() * get_M_nu() * dt * get_C_L_A() + M_sigma.block(0, 0, N_L, N_L))
		* get_Q_LopP_L()).pruned());
	Eigen::blockFill<double>(triplets, N_L, 0, 
		(M_sigma.block(N_L, 0, tN_pA, N_L) * get_Q_LopP_L()).pruned());
	Eigen::blockFill<double>(triplets, 0, N_Lo + N_pP,
		(- get_tC_pL_A()).pruned());
	Eigen::blockFill<double>(triplets, 0, N_Lo + N_pP + N_pL, 
		(M_sigma.block(0, N_L, N_L, tN_pA)).pruned());
	Eigen::blockFill<double>(triplets, N_L, N_Lo + N_pP + N_pL,
		(lambda2 / dt * Eigen::MatrixXd::Identity(tN_pA, tN_pA) + M_sigma.block(N_L, N_L, tN_pA, tN_pA)).pruned());
	Eigen::blockFill<double>(triplets, N_L + tN_pA, N_Lo,
		get_S_pP_pm());
	// Eigen::blockFill<double>(triplets, N_L + tN_pA + N_pP_pm, N_Lo + N_pP,
	//	get_tC_pL_AI());
	Eigen::blockFill<double>(triplets, N_L + tN_pA + N_pP_pm, N_Lo + N_pP + N_pL,
		get_tS_pA_AI());

	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();
	// Eigen::countRowNNZ(mat);

	// Assemble the right vector
	vec.segment(0, N_L) = lambda2 / dt * get_M_eps() * e + get_tC_Lo_Ao() * get_M_nu() * b - j_aux.segment(0, N_L);
	vec.segment(N_L, tN_pA) = lambda2 / dt * dp - j_aux.segment(N_L, tN_pA);
	vec.segment(N_L + tN_pA, N_pP_pm / 2).setConstant(getPotential(time + dt));
	vec.segment(N_L + tN_pA + N_pP_pm / 2, N_pP_pm / 2).setConstant(0.0);
	vec.segment(N_L + tN_pA + N_pP_pm, tN_AI).setZero();

	// Solve
	std::cout << "-- Linear system assembled. Size = "<< mat.rows();
	std::cout << " nonZero = " << mat.nonZeros() << std::endl;
	static Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	// static Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

	Eigen::VectorXd sol;
	solver.compute(mat);
	if (solver.info() != Eigen::Success) {
		std::cout << "********************************" << std::endl;
		std::cout << "*   Linear system not valid!   *" << std::endl;
		std::cout << "********************************" << std::endl; 
		exit(EXIT_FAILURE);
	}
	std::cout << "-- Linear system fractorized. Start solving ... " << std::endl;
	sol = solver.solve(vec);
	if (((mat * sol - vec).cwiseAbs().array() < 1e-10).all()) {
		std::cout << "solution confirmed" << std::endl;
	} 
	else { // switch to solver_base
		std::cout << "failed to solve accurately." << std::endl;
	}
    std::cout << "max error = " <<  (mat * sol - vec).cwiseAbs().maxCoeff() << std::endl;	
	std::cout << "-- Linear system solved" << std::endl;

	eo  = sol.segment(0, N_Lo);
	phi = sol.segment(N_Lo, N_pP);
	hp  = sol.segment(N_Lo + N_pP, N_pL);
	dp  = sol.segment(N_Lo + N_pP + N_pL, tN_pA);

	// Update edge voltage at each primal edge
	Eigen::VectorXd temp(eo.size() + phi.size());
	temp << eo, phi;
	e = get_Q_LopP_L() * temp;

	// Update j
	temp.resize(e.size() + dp.size());
	temp << e, dp;
	j = M_sigma * temp + j_aux; 

	// Update E
	updateInterpolated_E();
}

void MaxwellSolver::evolveMagneticFlux(const double dt) {  
	b += - dt * get_C_L_A() * e;
	updateInterpolated_B();
}

Eigen::VectorXd MaxwellSolver::getNorms() const {
	double E_1norm = 0, E_2norm = 0;
	double B_1norm = 0, B_2norm = 0;
	for (int i = 0; i < dual->getNumberOfCells(); i++) {
		const double vol = dual->getCell(i)->getVolume();
		E_1norm += E.row(i).norm() * vol;
		E_2norm += E.row(i).squaredNorm() * vol;
		B_1norm += B.row(i).norm() * vol;
		B_2norm += B.row(i).squaredNorm() * vol;
	}
	E_2norm = std::sqrt(E_2norm);
	B_2norm = std::sqrt(B_2norm);
	Eigen::VectorXd norms(4);
	norms << E_1norm, E_2norm, B_1norm, B_2norm;
	return norms;
}







