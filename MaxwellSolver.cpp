#include "MaxwellSolver.h"
#include <chrono>


MaxwellSolver::MaxwellSolver()
{
}

MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator)
: primal(primal), dual(dual), interpolator(interpolator)
{	
	const int N_P      = primal->getNumberOfVertices();
	const int N_L      = primal->getNumberOfEdges();
	const int N_Lo     = primal->facet_counts.nE_interior;
	const int N_PI     = dual->facet_counts.nC_solid; 
	const int N_pPI    = primal->facet_counts.nV_insulating;
	const int N_A      = primal->getNumberOfFaces();
	const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;

	eo.resize(primal->facet_counts.nE_interior); eo.setZero();
	phi.resize(primal->facet_counts.nV_boundary); phi.setZero();
	hp.resize(primal->facet_counts.nE_boundary); hp.setZero();
	ho.resize(N_Ao); ho.setZero();
	dp.resize(dual->facet_counts.nF_boundary); dp.setZero();
	sol.resize(N_Lo + N_pPI + N_PI); sol.setZero();
	
	e.resize(primal->getNumberOfEdges()); e.setZero();
	b.resize(primal->getNumberOfFaces()); b.setZero();
	j.resize(dual->getNumberOfFaces()); j.setZero();

	// Define index mapping: phi_ins ---> primal boundary vertex
	phi_ins2pP.resize(primal->facet_counts.nV_insulating);
	{
		int idx = 0;
		for (const Vertex* v : primal->getVertices()) {
			if (v->getType() == Vertex::Type::Insulating) {
				assert(v->isBoundary());
				phi_ins2pP(idx++) = v->getIndex();
			}
		}
	}
	assert(idx == primal->facet_counts.nV_insulating);

	// Define index mapping: larange multiplier p ---> primal vertex at insulating boundary
	p2PI.resize(dual->facet_counts.nC_solid);
	{
		int idx = 0;
		for (const Vertex* v : primal->getVertices()) {
			const int v_idx = v->getIndex();
			if (dual->getCell(v_idx)->getFluidType() == Cell::FluidType::Solid) {
				p2PI(idx++) = v_idx;
			}
		}
	}
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

const Eigen::SparseMatrix<double>& MaxwellSolver::get_M_mu()
{
	if (M_mu.size() == 0) {
		const int N_A = primal->getNumberOfFaces();
		std::vector<T> triplets;
		triplets.reserve(N_A);
		for (int i = 0; i < N_A; i++) {
			const double fA = primal->getFace(i)->getArea();
			const double eL = dual->getEdge(i)->getLength();
			assert(dual->getEdge(i)->getVertexMid() == nullptr);
			const double value = fA / eL;
			triplets.emplace_back(T(i, i, value));
		}
		M_mu.resize(N_A, N_A);
		M_mu.setFromTriplets(triplets.begin(), triplets.end());
		M_mu.makeCompressed();
		assert(M_mu.size() > 0);
		assert(M_mu.nonZeros() > 0);
	}
	assert((M_mu * get_M_nu() - Eigen::MatrixXd::Identity(M_mu.rows(), M_mu.rows())).norm() < 1e-7);
	return M_mu;
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

const Eigen::SparseMatrix<double>& MaxwellSolver::get_solidDiv() {
	if (solidDiv.size() == 0) {
		int count = 0;
		std::vector<T> triplets;
		triplets.reserve(dual->facet_counts.nC_solid * 8);
		for (const Cell* c : dual->getCells()) {
			if (c->getFluidType() == Cell::FluidType::Solid) {
				for (const Face* f : c->getFaceList()) {
					triplets.emplace_back(T(count, f->getIndex(), c->getOrientation(f)));
				}
				count++;
			}
		}
		assert(count == dual->facet_counts.nC_solid);
		solidDiv.resize(count, dual->getNumberOfFaces());
		solidDiv.setFromTriplets(triplets.begin(), triplets.end());
		solidDiv.makeCompressed();
	}
	return solidDiv;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_solidGrad() {
	if (solidGrad.size() == 0) {
		int count = 0;
		std::vector<T> triplets;
		triplets.reserve(2 * (dual->facet_counts.nF_wall + dual->facet_counts.nF_undefined));
		for (const Vertex* v : primal->getVertices()) {
			const Cell* c = dual->getCell(v->getIndex());
			if (c->getFluidType() == Cell::FluidType::Solid) {
				for (const Face* f : c->getFaceList()) {
					if (!f->isBoundary()) {
						triplets.emplace_back(T(f->getIndex(), count, c->getOrientation(f)));
					}
				}
				count++;
			}
		}
		assert(count == dual->facet_counts.nC_solid);
		solidGrad.resize(primal->getNumberOfEdges(), count);
		solidGrad.setFromTriplets(triplets.begin(), triplets.end());
		solidGrad.makeCompressed();
	}
	return solidGrad;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_harmonicE() {
	if (harmonicE.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(dual->facet_counts.nF_undefined);
		for (const Edge* e : primal->getEdges()) {
			if (dual->getFace(e->getIndex())->getFluidType() == Face::FluidType::Undefined) {
				if (std::abs(e->getDirection()(2)) < 1e-10) {
					double a1, b1, a2, b2, A, B, C, integral, temp;
					a1 = e->getVertexA()->getPosition()(0);
					b1 = e->getVertexA()->getPosition()(1);
					a2 = e->getVertexB()->getPosition()(0);
					b2 = e->getVertexB()->getPosition()(1);
					temp = (a2 - a1) * (a2 - a1) + (b2 - b1) * (b2 - b1);
					A = (- a2*b1 + a1*b2) / temp;
					B = (a1*(a2 - a1) + b1*(b2 - b1)) / temp;
					C = (a1*a1 + b1*b1) / temp - B*B; // In theory, C >= 0
					if (std::abs(C) < 1e-12) {
						integral = 0;
					}
					else {
						integral = A / std::sqrt(C) * (std::atan((1 + B) / std::sqrt(C)) - std::atan(B / std::sqrt(C)));
					}
					triplets.emplace_back(T(0, e->getIndex(), integral));
				}
			}
		}
		harmonicE.resize(1, primal->getNumberOfEdges());
		harmonicE.setFromTriplets(triplets.begin(), triplets.end());
		harmonicE.makeCompressed();
	}
	return harmonicE;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_harmonicD() {
	if (harmonicD.size() == 0) {
		harmonicD = get_harmonicE() * get_M_eps();
	}
	return harmonicD;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_DirichletHarmonicD() {
	if (DirichletHarmonicD.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(dual->facet_counts.nF_undefined);
		for (const Face* f : dual->getFaces()) {
			if (f->getFluidType() != Face::FluidType::Interior && f->getFluidType() != Face::FluidType::Opening) {
				double integral = 0;
				if (f->getSubFaceList().size() > 1) {
					for (const Face* subf : f->getSubFaceList()) {
						if (std::abs(subf->getNormal()(2)) < 1e-10
						&& std::abs(subf->getNormal().segment(0,2).dot(subf->getCenter().segment(0,2).normalized())) > 0.8) {
							if (subf->getNormal().segment(0,2).dot(subf->getCenter().segment(0,2)) > 0) {
								integral += subf->getArea() / subf->getCenter().segment(0,2).norm();
							}
							else {
								integral -= subf->getArea() / subf->getCenter().segment(0,2).norm();
							}
						}
					}
				}
				else {
					if (std::abs(f->getNormal()(2)) < 1e-10 
						&& std::abs(f->getNormal().segment(0,2).dot(f->getCenter().segment(0,2).normalized())) > 0.8) {
						if (f->getNormal().segment(0,2).dot(f->getCenter().segment(0,2)) > 0) {
							integral = f->getArea() / f->getCenter().segment(0,2).norm();
						}
						else {
							integral = - f->getArea() / f->getCenter().segment(0,2).norm();
						}
					}
				}
				triplets.emplace_back(T(0, f->getIndex(), integral));
			}
		}
		DirichletHarmonicD.resize(1, dual->getNumberOfFaces());
		DirichletHarmonicD.setFromTriplets(triplets.begin(), triplets.end());
		DirichletHarmonicD.makeCompressed();
	}
	return DirichletHarmonicD;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_DirichletHarmonicE() {
	if (DirichletHarmonicE.size() == 0) {
		DirichletHarmonicE = Eigen::MatrixXd(get_M_eps()).inverse() * get_DirichletHarmonicD().row(0).segment(0, primal->getNumberOfEdges());
	}
	return DirichletHarmonicE;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_M_sigma_const() {
	if (M_sigma_const.size() == 0) {
		std::vector<T> triplets;
		triplets.reserve(dual->getNumberOfFaces());
		Eigen::blockFill<double>(triplets, 0, 0, get_M_eps());
		Eigen::blockFill<double>(triplets, primal->getNumberOfEdges(), primal->getNumberOfEdges(), 
			Eigen::MatrixXd::Identity(dual->facet_counts.nF_boundary, dual->facet_counts.nF_boundary).sparseView());
		Eigen::SparseMatrix<double> temp(dual->getNumberOfFaces(), dual->getNumberOfFaces());
		temp.setFromTriplets(triplets.begin(), triplets.end());

		std::vector<T> triplets2;
		triplets.reserve(dual->getNumberOfFaces());
		for (const Face* f : dual->getFaces()) {
			if (f->getFluidType() == Face::FluidType::Interior || f->getFluidType() == Face::FluidType::Opening) {
				triplets2.emplace_back(T(f->getIndex(), f->getIndex(), 0.0));
			}
		}
		Eigen::SparseMatrix<double> temp2(dual->getNumberOfFaces(), dual->getNumberOfFaces());
		temp2.setFromTriplets(triplets2.begin(), triplets2.end());
		M_sigma_const = temp2 * temp;
	}
	return M_sigma_const;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_glb2lcl() {
	if (glb2lcl.size() == 0) {
		std::vector<T> triplets;
		for (int i = 0; i < dual->getNumberOfFaces(); i++) {
				const Face * f = dual->getFace(i);
				if (f->isBoundary()) {
					triplets.emplace_back(T(i, i, 1.0 / f->getArea()));
				}
				else {
					triplets.emplace_back(T(i, i, 1.0 / primal->getEdge(i)->getLength()));
				}
		}
		glb2lcl.resize(dual->getNumberOfFaces(), dual->getNumberOfFaces());
		glb2lcl.setFromTriplets(triplets.begin(), triplets.end());
	}
	return glb2lcl;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_pPI_Lo() {  // This is actually "-grad".
	if (G_pPI_Lo.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_pPI     = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		std::vector<T> triplets;
		for (int phi_ins_idx = 0; phi_ins_idx < N_pPI; phi_ins_idx++) {
			const Vertex* v = primal->getVertex(phi_ins2pP(phi_ins_idx));
			assert(v->getType() == Vertex::Type::Insulating);
			for (const Edge* e_v: v->getEdges()) {
				if (!e_v->isBoundary()) {
					assert(e_v->getIndex() < N_Lo);
					triplets.emplace_back(e_v->getIndex(), phi_ins_idx, e_v->getIncidence(v));
				}
			}
		}
		G_pPI_Lo.resize(N_Lo, N_pPI);
		G_pPI_Lo.setFromTriplets(triplets.begin(), triplets.end());
		G_pPI_Lo.makeCompressed();
	}
	return G_pPI_Lo;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_D_Lo_pPI() {
	if (D_Lo_pPI.rows() == 0) {
		D_Lo_pPI = get_G_pPI_Lo().transpose();
	}
	return D_Lo_pPI;
}


const Eigen::SparseMatrix<double>& MaxwellSolver::get_C_Ao_Lo() {
	if (C_Ao_Lo.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_pPI     = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		std::vector<T> triplets;
		for (const Face* f : primal->getFaces()) {
			if (!f->isBoundary()) {
				for (const Edge* e_f : f->getEdgeList()) {
					if(!e_f->isBoundary()) {
						triplets.emplace_back(e_f->getIndex(), f->getIndex(), f->getOrientation(e_f));
						assert(e_f->getIndex() < N_Lo);
						assert(f->getIndex() < N_Ao);
					}
				}
			}
		}
		C_Ao_Lo.resize(N_Lo, N_Ao);
		C_Ao_Lo.setFromTriplets(triplets.begin(), triplets.end());
		C_Ao_Lo.makeCompressed();
	}
	return C_Ao_Lo;
}
const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_pPI_L() {
	if (G_pPI_L.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_pPI     = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		
		std::vector<T> triplets;
		for (int phi_ins_idx = 0; phi_ins_idx < N_pPI; phi_ins_idx++) {
			const Vertex* v = primal->getVertex(phi_ins2pP(phi_ins_idx));
			assert(v->getType() == Vertex::Type::Insulating);
			for (const Edge* e_v: v->getEdges()) {
				triplets.emplace_back(e_v->getIndex(), phi_ins_idx, e_v->getIncidence(v));
			}
		}
		G_pPI_L.resize(N_L, N_pPI);
		G_pPI_L.setFromTriplets(triplets.begin(), triplets.end());
		G_pPI_L.makeCompressed();
	}
	return G_pPI_L;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_P_L() {
	if (G_P_L.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_pPI     = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		
		std::vector<T> triplets;
		for (int P_idx = 0; P_idx < N_P; P_idx++) {
			const Vertex* v = primal->getVertex(P_idx);
			for (const Edge* e_v: v->getEdges()) {
				triplets.emplace_back(e_v->getIndex(), P_idx, e_v->getIncidence(v));
			}
		}
		G_P_L.resize(N_L, N_P);
		G_P_L.setFromTriplets(triplets.begin(), triplets.end());
		G_P_L.makeCompressed();
	}
	return G_P_L;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_D_L_pPI() {
	if (D_L_pPI.rows() == 0) {
		D_L_pPI = G_pPI_L.transpose();
	}
	return D_L_pPI;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_PI_L() {
	if (G_PI_L.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_PI     = dual->facet_counts.nC_solid; 
		const int N_pPI    = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		std::vector<T> triplets;
		for (int i = 0; i < N_PI; i++) {
			const Vertex* v = primal->getVertex(p2PI(i));
			for (const Edge* e_v : v->getEdges()) {
				triplets.emplace_back(T(e_v->getIndex(), i, e_v->getIncidence(v)));
			}
		}
		G_PI_L.resize(N_L, N_PI);
		G_PI_L.setFromTriplets(triplets.begin(), triplets.end());
		G_PI_L.makeCompressed();
	}
	return G_PI_L;
}

const Eigen::SparseMatrix<double>& MaxwellSolver::get_G_PI_Lo() {
	if (G_PI_Lo.rows() == 0) {
		const int N_P      = primal->getNumberOfVertices();
		const int N_L      = primal->getNumberOfEdges();
		const int N_Lo     = primal->facet_counts.nE_interior;
		const int N_PI     = dual->facet_counts.nC_solid; 
		const int N_pPI    = primal->facet_counts.nV_insulating;
		const int N_A      = primal->getNumberOfFaces();
		const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
		G_PI_Lo = get_G_PI_L().topRows(N_Lo);
	}
	return G_PI_Lo;
}

const Eigen::MatrixXd& MaxwellSolver::updateInterpolated_E() {
	Eigen::VectorXd temp(e.size() + dp.size());
	temp << e, dp; // concatenate [e, dp]^T
	E = interpolator->get_E_interpolator().oneContract(temp);
	return E;
}

const Eigen::MatrixXd& MaxwellSolver::updateInterpolated_B() {
	Eigen::VectorXd temp(b.size() + hp.size());
	temp << b, hp;  // concatenate [b, hp]^T
	B = interpolator->get_B_interpolator().oneContract(temp);
	return B;
}

const Eigen::MatrixXd& MaxwellSolver::getInterpolated_E() const {
	return E;
}

const Eigen::MatrixXd& MaxwellSolver::getInterpolated_B() const {
	return B;
}

const Eigen::VectorXd& MaxwellSolver::get_e_vec() const {
	return e;
}

const Eigen::VectorXd& MaxwellSolver::get_dp_vec() const {
	return dp;
}

void MaxwellSolver::solveLinearSystem(const double time, 
									  const double dt, 
									  Eigen::SparseMatrix<double>&& M_sigma, 
									  Eigen::VectorXd&& j_aux) {
	//M_sigma = get_M_sigma_const();
	//j_aux.setZero();

	std::cout << "-- Start assembling linear system ..." << std::endl;
	const int N_Lo    = primal->getNumberOfEdges() - primal->facet_counts.nE_boundary;
	const int N_pP    = primal->facet_counts.nV_boundary;
	const int N_pL    = primal->facet_counts.nE_boundary;
	const int tN_pA   = dual->facet_counts.nF_boundary;
	const int N_L     = primal->getNumberOfEdges();
	const int N_pP_pm = primal->facet_counts.nV_electrode;
	const int tN_AI   = primal->facet_counts.nV_insulating;
	const int tN_sV	  = dual->facet_counts.nC_solid;
	assert(N_Lo + N_pP + N_pL + tN_pA == N_L + tN_pA + N_pP_pm + tN_AI);
	Eigen::SparseMatrix<double> mat(N_Lo + N_pP + N_pL + tN_pA + tN_sV, N_Lo + N_pP + N_pL + tN_pA + tN_sV);
	std::vector<Eigen::Triplet<double>> triplets;
	Eigen::VectorXd vec(N_Lo + N_pP + N_pL + tN_pA + tN_sV);
	vec.setZero();
	const double lambda2 = parameters.lambdaSquare;

	// M_sigma += get_M_sigma_const();
	// modifyM_sigma(M_sigma);
	checkM_sigma(M_sigma);

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
	Eigen::blockFill<double>(triplets, N_L + tN_pA + N_pP_pm, N_Lo + N_pP + N_pL,
		get_tS_pA_AI());

	Eigen::blockFill<double>(triplets, N_Lo + N_pP + N_pL + tN_pA, 0, 
		get_solidDiv().leftCols(N_L) * get_M_eps() * get_Q_LopP_L());
	Eigen::blockFill<double>(triplets, N_Lo + N_pP + N_pL + tN_pA, N_Lo + N_pP + N_pL,
		get_solidDiv().rightCols(tN_pA));
	Eigen::blockFill<double>(triplets, 0, N_Lo + N_pP + N_pL + tN_pA,
		(get_M_eps() * get_solidGrad()));
	
	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();
	std::cout << "-- Extended linear system assembled. Size = (" << mat.rows() << ", " << mat.cols() << ").";
	std::cout << " nonZero = " << mat.nonZeros() << std::endl;

	// Assemble the rhs vector
	vec.segment(0, N_L) = lambda2 / dt * get_M_eps() * e + get_tC_Lo_Ao() * get_M_nu() * b - j_aux.segment(0, N_L);
	vec.segment(N_L, tN_pA) = lambda2 / dt * dp - j_aux.segment(N_L, tN_pA);
	vec.segment(N_L + tN_pA, N_pP_pm / 2).setConstant(getPotential(time + dt));
	vec.segment(N_L + tN_pA + N_pP_pm / 2, N_pP_pm / 2).setConstant(0.0);
	vec.segment(N_L + tN_pA + N_pP_pm, tN_AI).setZero();
	vec.segment(N_L + tN_pA + N_pP_pm + tN_AI, tN_sV).setZero();

	// Solve
	static Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::VectorXd sol;
	solver.compute(mat);
	if (solver.info() != Eigen::Success) {
		std::cout << "*****************************************" << std::endl;
		std::cout << "*         Linear system not valid!      *" << std::endl;
		std::cout << "*****************************************" << std::endl; 
		exit(EXIT_FAILURE);
	}
	else {
		std::cout << "-- Linear system fractorized. Start solving ... " << std::endl;
		sol = solver.solve(vec);
		if (((mat * sol - vec).cwiseAbs().array() < 1e-10).all()) {
			std::cout << "-- Solution is confirmed." << std::endl;
		} 
		else { 
			std::cout << "-- Fail to solve accurately." << std::endl;
		}
		std::cout << "-- max error = " <<  (mat * sol - vec).cwiseAbs().maxCoeff() << std::endl;	
		std::cout << "-- Standard Linear system solved" << std::endl;
	}
	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << "Elapsed time: " << (std::chrono::duration_cast<std::chrono::seconds>(stop - start)).count() << std::endl;

	std::cout << "-- max p = " << (sol.segment(N_Lo + N_pP + N_pL + tN_pA, tN_sV).cwiseAbs().maxCoeff()) << std::endl;
	
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
	checkZeroDiv();
}


void MaxwellSolver::solveLinearSystem_sym(const double time, 
									  	  const double dt, 
									      Eigen::SparseMatrix<double>&& M_sigma, 
									      Eigen::VectorXd&& j_aux) {

	
	//M_sigma = get_M_sigma_const();
	//j_aux.setZero();

	std::cout << "-- Start assembling linear system sym ..." << std::endl;
	const int N_P      = primal->getNumberOfVertices();
	const int N_L      = primal->getNumberOfEdges();
	const int N_Lo     = primal->facet_counts.nE_interior;
	const int N_PI     = dual->facet_counts.nC_solid; 
	const int N_pPI    = primal->facet_counts.nV_insulating;
	const int N_A      = primal->getNumberOfFaces();
	const int N_Ao     = primal->getNumberOfFaces() - primal->facet_counts.nF_boundary;
	
	// Reconstruct eo, phi_ins, ho at the old time step
	std::cout << "Reconstruct eo, phi_ins, ho at the old time step ..." << std::endl;
	Eigen::VectorXd eo_old = e.head(N_Lo);
	Eigen::VectorXd e_old  = e;
	Eigen::VectorXd ho_old = (get_M_nu() * b).head(N_Ao); 

	// Some auxiliary matrices
	Eigen::SparseMatrix<double> M_eps_int   = get_M_eps().block(0,0,N_Lo,N_Lo);
	Eigen::SparseMatrix<double> M_sigma_int = M_sigma.block(0,0,N_Lo,N_Lo);
	Eigen::SparseMatrix<double> M_mu_int    = get_M_mu().block(0,0,N_Ao,N_Ao);

	// Reserve space for lse
	std::cout << "Reserve space for lse ..." << std::endl;
	Eigen::SparseMatrix<double> mat(N_Lo + N_pPI + N_PI + N_Ao, N_Lo + N_pPI + N_PI + N_Ao);
	std::vector<Eigen::Triplet<double>> triplets_lse;
	Eigen::VectorXd vec(N_Lo + N_pPI + N_PI + N_Ao);
	const double lambda2 = parameters.lambdaSquare;

	// Assemble the square matrix
	std::cout << "Assemble the square matrix ..." << std::endl;
	Eigen::blockFill<double>(triplets_lse, 0, 0, 
		(lambda2 / dt * M_eps_int + M_sigma_int).pruned());
	Eigen::blockFill<double>(triplets_lse, 0, N_Lo, 
		((lambda2 / dt * M_eps_int + M_sigma_int) * get_G_pPI_Lo()).pruned());
	Eigen::blockFill<double>(triplets_lse, 0, N_Lo + N_pPI,
		(M_eps_int * get_G_PI_Lo()).pruned());
	Eigen::blockFill<double>(triplets_lse, 0, N_Lo + N_pPI + N_PI,
		(- get_C_Ao_Lo()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo, 0, 
		(get_D_Lo_pPI() * (lambda2 / dt * M_eps_int + M_sigma_int)).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo, N_Lo,
		(get_D_L_pPI()  * (lambda2 / dt * M_eps + M_sigma.block(0,0,N_L,N_L)) * get_G_pPI_L()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo, N_Lo + N_pPI,
		(get_D_L_pPI()  * M_eps * get_G_PI_L()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo + N_pPI, 0,
		((M_eps_int * get_G_PI_Lo()).transpose()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo + N_pPI, N_Lo,
		((get_D_L_pPI()  * M_eps * get_G_PI_L()).transpose()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo + N_pPI + N_PI, 0,
		(get_C_Ao_Lo().transpose()).pruned());
	Eigen::blockFill<double>(triplets_lse, N_Lo + N_pPI + N_PI, N_Lo + N_pPI + N_PI,
		(1. / dt * M_mu_int).pruned());
	
	mat.setFromTriplets(triplets_lse.begin(), triplets_lse.end());
	mat.makeCompressed();

	// Construct Gpsi
	std::cout << "Construct Gpsi ..." << std::endl;
	Eigen::VectorXd psi_glb(N_P);
	psi_glb.setZero();
	for (const Vertex* v : primal->getVertices()) {
		if (v->getType() == Vertex::Type::Electrode && v->getPosition()[2] > 0) {
			assert(v->isBoundary());
			psi_glb(v->getIndex()) = 1.0;
		}
	}
	Eigen::VectorXd Gpsi = get_G_P_L() * psi_glb;

	// Assemble the rhs vector
	std::cout << "Assemble the rhs vector ..." << std::endl;
	vec.setZero();
	vec.segment(0, N_Lo) = lambda2 / dt * M_eps_int * eo_old
						   - j_aux.head(N_Lo) 	
						   - (lambda2 / dt * M_eps_int + M_sigma_int) * Gpsi.head(N_Lo);
	vec.segment(N_Lo, N_pPI) = get_D_L_pPI() * (lambda2 / dt * M_eps * e_old - j_aux.head(N_L)) 
	                                 - D_L_pPI * (lambda2 / dt * M_eps + M_sigma.block(0,0,N_L,N_L)) * Gpsi;  
	vec.segment(N_Lo + N_pPI, N_PI) = - G_PI_L.transpose() * M_eps * Gpsi;
	vec.segment(N_Lo + N_pPI + N_PI, N_Ao) = 1./dt * M_mu_int * ho_old;

	// Schur complement
	Eigen::SparseMatrix<double> A{mat.topLeftCorner(N_Lo + N_pPI + N_PI, N_Lo + N_pPI + N_PI)};
	Eigen::SparseMatrix<double> invB{dt * M_mu_int.cwiseInverse()};
	Eigen::SparseMatrix<double> C{-mat.block(0, N_Lo + N_pPI + N_PI, N_Lo + N_pPI + N_PI, N_Ao)};
	Eigen::SparseMatrix<double> mat_reduced{A + Eigen::SparseMatrix<double>(C*invB*C.transpose())};

	Eigen::VectorXd a{vec.segment(0, N_Lo + N_pPI + N_PI)};
	Eigen::VectorXd b{vec.segment(N_Lo + N_pPI + N_PI, N_Ao)};

	Eigen::VectorXd vec_reduced{a + C * invB * b};

	std::cout << "-- linear system assembled. Size = (" << mat_reduced.rows() << ", " << mat_reduced.cols() << ").";
	std::cout << " nonZero = " << mat_reduced.nonZeros() << std::endl;

	// Solve
	std::cout << "Solve ..." << std::endl;
	//static Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	static Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	
	auto start = std::chrono::high_resolution_clock::now();
	solver.compute(mat_reduced);
	if (solver.info() != Eigen::Success) {
		std::cout << "*****************************************" << std::endl;
		std::cout << "*         Linear system not valid!      *" << std::endl;
		std::cout << "*****************************************" << std::endl; 
		exit(EXIT_FAILURE);
	}
	else {
		std::cout << "-- Linear system fractorized. Start solving ... " << std::endl;
		sol = solver.solveWithGuess(vec_reduced, sol);
		if (((mat_reduced * sol - vec_reduced).cwiseAbs().array() < 1e-10).all()) {
			std::cout << "-- Solution is confirmed." << std::endl;
		} 
		else { 
			std::cout << "-- Fail to solve accurately." << std::endl;
		}
		std::cout << "-- max error = " <<  (mat_reduced * sol - vec_reduced).cwiseAbs().maxCoeff() << std::endl;	
		std::cout << "-- Standard Linear system solved" << std::endl;
	}
	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << "Elapsed time: " << (std::chrono::duration_cast<std::chrono::seconds>(stop - start)).count() << std::endl;
	
	Eigen::VectorXd eo_new = sol.segment(0, N_Lo);
	Eigen::VectorXd phi_ins_new = sol.segment(N_Lo, N_pPI);
	Eigen::VectorXd ho_new = invB * (b - C.transpose() * sol);

	// Reconstruct the variables [eo, phi] according to the old definition
	std::cout << "Reconstruct the variables [eo, phi] according to the old definition" << std::endl;
	eo  = eo_new + G_pPI_Lo * phi_ins_new + Gpsi.head(N_Lo);
	phi.setZero();
	for (int phi_ins_idx = 0; phi_ins_idx < N_pPI; phi_ins_idx++) {
		phi(ppP2phi(phi_ins2pP(phi_ins_idx))) = phi_ins_new(phi_ins_idx);
	}
	for (int phi_idx = 0; phi_idx < phi.size(); phi_idx++) {
		if (primal->getVertex(phi2ppP(phi_idx))->getType() == Vertex::Type::Electrode
			&& primal->getVertex(phi2ppP(phi_idx))->getPosition()[2] > 0) {
			assert(abs(phi(phi_idx)) < 1e-10);
			phi(phi_idx) = 1.0; // constant potential applied on anode
		}
	}

	// Update edge voltage at each primal edge
	std::cout << "Update edge voltage at each primal edge" << std::endl;
	Eigen::VectorXd temp(eo.size() + phi.size());
	temp << eo, phi;
	e = get_Q_LopP_L() * temp;

	// Reconstruct hp
	std::cout << " Reconstruct hp" << std::endl;
	assert(hp.size() == primal->facet_counts.nE_boundary);
	for (int hp_idx = 0; hp_idx < primal->facet_counts.nE_boundary; hp_idx++) {
		const Edge* dual_e = dual->getEdge(ph2dpL(hp_idx));
		assert(dual_e->isBoundary());
		const Face* f_e = nullptr;
		int count = 0; // just for checking
		for (const Face* subf_e : dual_e->getFaceList()) {
			if (!subf_e->isBoundary()) {
				f_e = subf_e;
				count++;
			}
		}
		assert(count == 1);
		const int idx_f_e = f_e->getIndex(); // which is also the index of the corresponding primal boundary edge
		assert((f_e->getNormal().cross(primal->getEdge(idx_f_e)->getDirection())).norm() < 1e-7);
		assert(primal->getEdge(idx_f_e)->isBoundary());
		double tmp = 0;
		tmp += lambda2 * M_eps.coeff(idx_f_e, idx_f_e) * (e(idx_f_e) - e_old(idx_f_e)) / dt;
		tmp += M_sigma.coeff(idx_f_e, idx_f_e) * e(idx_f_e);
		tmp += j_aux(idx_f_e);
		for (const Edge* sube_f : f_e->getEdgeList()) {
			if(!sube_f->isBoundary() && sube_f->getIndex() < N_Ao) { // Note that here we implicity use the knowledge that boundary normal component of B is zero
				tmp -= f_e->getOrientation(sube_f) * ho_new(sube_f->getIndex());
			}
		}
		hp(hp_idx) = tmp / f_e->getOrientation(dual_e);
	}

	// Check Ampere's law
	/*
	Eigen::VectorXd ho_ext(N_A);
	ho_ext.setZero();
	ho_ext.head(N_Ao) = ho_new;
	for (const Face* f : dual->getFaces()) {
		if (!f->isBoundary() && !primal->getEdge(f->getIndex())->isBoundary()) {
			const int e_idx = f->getIndex();
			const double lhs = lambda2 * M_eps_int.coeff(e_idx, e_idx) * (e(e_idx) - e_old(e_idx)) / dt;
			double rhs = - M_sigma_int.coeff(e_idx, e_idx) * e(e_idx) - j_aux(e_idx);
			for (const Edge* edge : f->getEdgeList()) {
				if (edge->isBoundary()) {
					rhs += hp(dpL2ph(edge->getIndex())) * f->getOrientation(edge);
				}
				else {
					rhs += ho_ext(edge->getIndex()) * f->getOrientation(edge);
					if (edge->getIndex() >= N_Ao) assert(abs(ho_ext(edge->getIndex())) < 1e-7);
				}
			}
			if (abs(lhs - rhs) > 1e-10) {
				std::cout << "residual of Ampere's law at dual face " << f->getIndex() << " is " << abs(lhs - rhs) 
				          << " of type ";
				if (f->getFluidType() == Face::FluidType::Interior) std::cout << "interior" << std::endl;  
				if (f->getFluidType() == Face::FluidType::Mixed)    std::cout << "mixed" << std::endl;
				if (f->getFluidType() == Face::FluidType::Opening)  std::cout << "opening" << std::endl;
				if (f->getFluidType() == Face::FluidType::Wall)     std::cout << "wall" << std::endl;      
			}
		}
	}*/

	// Check if cross(H)`n is zero
	/*
	for (const Face* f : dual->getFaces()) {
		if (f->isBoundary()) {
			double sum = 0;
			for (const Edge* e : f->getEdgeList()) {
				sum += hp(dpL2ph(e->getIndex())) * f->getOrientation(e);
			}
			if (abs(sum) > 1e-10 && f->getFluidType() != Face::FluidType::Opening) std::cout << "cross(H)`n = 0 violated at dual Face" << f->getIndex() << ":" << abs(sum) << std::endl;
		}
	}*/

	// Reconstruct dp
	// ATTENTION: we set d associated to dual faces lying at insulating boundary to zero 
	// such that this boundary condition is preserved at the limit.
	std::cout << " Reconstruct dp" << std::endl;
	Eigen::VectorXd dp_old = dp;
	dp.setZero();
	for (int dp_idx = 0; dp_idx < dual->facet_counts.nF_boundary; dp_idx++) {
		const Face* dual_f = dual->getFace(pd2dpA(dp_idx));
		assert(dual_f->isBoundary());
		if (dual_f->getCellList()[0]->getFluidType() != Cell::FluidType::Solid) { 
			const int idx_f = dual_f->getIndex();
			double tmp = 0;
			tmp -= lambda2 / dt * dp_old(dp_idx);
			tmp += j_aux(dual_f->getIndex());
			for (const Edge* sube_f : dual_f->getEdgeList()) {
				assert(sube_f->isBoundary());
				tmp -= dual_f->getOrientation(sube_f) * hp(dpL2ph(sube_f->getIndex()));
			}
			dp(dp_idx) = tmp / (- lambda2 / dt - M_sigma.coeff(idx_f, idx_f));
		}
	}

	// Update j
	temp.resize(e.size() + dp.size());
	temp << e, dp;
	j = M_sigma * temp + j_aux; 

	// Update E
	updateInterpolated_E();
	checkZeroDiv();

	/*
	int count = 0;
	for (int idx_dp = 0; idx_dp < dp.size(); idx_dp++) {
		const Face* f = dual->getFace(pd2dpA(idx_dp));
		if (!f->isBoundary()) std::cout << "--------------- wrong! ---------------" << std::endl;
		const Cell* c = f->getCellList()[0];
		if (c->getFluidType() == Cell::FluidType::Solid) {
			count++;
			//std::cout << "D`n at insulating boundary = " << dp(idx_dp) << std::endl;
		}
	}
	if (count != N_pPI) std::cout << "Something is wrong!!!!" << std::endl;
	*/
}

void MaxwellSolver::evolveMagneticFlux(const double dt) {  
	b += - dt * get_C_L_A() * e;
	updateInterpolated_B();
}

Eigen::VectorXd MaxwellSolver::getD() {
	Eigen::VectorXd temp(dual->getNumberOfFaces());
	temp << get_M_eps() * e, dp;
	return temp;
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

void MaxwellSolver::checkZeroDiv() {
	if (((get_solidDiv() * getD()).cwiseAbs().array() < 1e-10).all()) {
		std::cout << "========= Zero divergence at solid is checked" << std::endl;
		return;
	} 
	else {
		std::cout << "=========================================================" << std::endl;
		std::cout << "||       Zero divergence at solid id not satisfied!    ||" << std::endl;
		std::cout << "=========================================================" << std::endl;
		std::cout << "-- max divD = " << (get_solidDiv() * getD()).cwiseAbs().maxCoeff() << std::endl;
	} 
	/*
	else {
		for (const Vertex* v : primal->getVertices()) {
			const Cell* dual_cell = dual->getCell(v->getIndex()); 
			if (dual_cell->getFluidType() == Cell::FluidType::Solid) {
				double sum = 0;
				for (const Face* f : dual_cell->getFaceList()) {
					sum += getD()(f->getIndex());
				}
				if (abs(sum) > 1e-5) {
					std::cout << "DivD != 0 at vertex " << v->getIndex() << " of type ";
					if (v->isBoundary()) {
						std::cout << " boundary" << std::endl;
					}
					else {
						std::cout << " interior" << std::endl;
					}
				}
			}
		}
		std::cout << "**********************************************" << std::endl;
		std::cout << "* Zero divergence at solid is not fulfilled  *" << std::endl;
		std::cout << "**********************************************" << std::endl;
		return;
	}*/
}

void MaxwellSolver::enforceHarmonicE() {
	
	e = Eigen::MatrixXd(get_harmonicE()).row(0);
	// std::cout << e << std::endl;
	updateInterpolated_E();
}

void MaxwellSolver::enforceDirichletHarmonicE() {
	
	dp = Eigen::MatrixXd(get_DirichletHarmonicD().row(0).segment(primal->getNumberOfEdges(), dual->facet_counts.nF_boundary)); 
	for (int i = 0; i < primal->getNumberOfEdges(); i++) {
		e(i) = get_DirichletHarmonicD().coeff(0,i) / get_M_eps().coeff(i,i);
	}
	// std::cout << e << std::endl;
	updateInterpolated_E();
}

void MaxwellSolver::checkM_sigma(Eigen::SparseMatrix<double> &M_sigma) const {
	double min = 999;
	double max = 0;
	for (int i = 0; i < dual->getNumberOfFaces(); i++) {
		double temp = M_sigma.coeff(i,i);
		if (dual->getFace(i)->getFluidType() == Face::FluidType::Interior || dual->getFace(i)->getFluidType() == Face::FluidType::Opening) {
			min = (temp < min) ? temp : min;
		}
		if (dual->getFace(i)->getFluidType() == Face::FluidType::Undefined || dual->getFace(i)->getFluidType() == Face::FluidType::Wall) {
			max = (temp > max) ? temp : max;
		}
	}
	std::cout << "min diagonal coeff of M_sigma at interior  is " << min << std::endl;
	std::cout << "max diagonal coeff of M_sigma at insulator is " << max << std::endl;
	return;
}

void MaxwellSolver::modifyM_sigma(Eigen::SparseMatrix<double> &M_sigma) const {
	std::vector<T> triplets;
	for (int i = 0; i < dual->getNumberOfFaces(); i++) {
		double row_sum = M_sigma.row(i).sum();
		triplets.emplace_back(i, i, row_sum + 1e-10);
	}
	M_sigma.setFromTriplets(triplets.begin(), triplets.end());
}

void MaxwellSolver::outputM_sigma(Eigen::SparseMatrix<double> &M_sigma) const {

}









