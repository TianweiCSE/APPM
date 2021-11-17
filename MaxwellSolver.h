#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "H5Writer.h"
#include "FluidSolver.h"
#include "Interpolator.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;

class MaxwellSolver
{
public:
	struct MaxwellParams {
		double lambdaSquare = 1.0;
	} parameters;

	MaxwellSolver();
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual);
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, MaxwellParams & param);
	~MaxwellSolver();

	bool isMaxwellCurrentSource = false;

	void updateMaxwellState(const double dt, const double time);

	void writeStates(H5Writer & writer) const;

	const Eigen::VectorXd & getBstate() const;

	void setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector);
	void setAzimuthalMagneticFluxField();
	void setTorusCurrent(const double x1, const double x2, const double z1, const double z2);

	Eigen::VectorXd electricPotentialTerminals(const double time);

	void solveLinearSystem(const double dt, Eigen::SparseMatrix<double> &&M_sigma, Eigen::VectorXd &&j_aux);
	void timeStepping(const double dt);

	Eigen::MatrixXd&& getInterpolated_E() const;
	Eigen::MatrixXd&& getInterpolated_B() const;   


protected:
	const PrimalMesh* primal = nullptr;
	const DualMesh*   dual   = nullptr;
	Interpolator interpolator;

	Eigen::VectorXd maxwellState;

	Eigen::VectorXd B_h, E_h, H_h, J_h;

	// Compute M_nu := M_mu^{-1} of size (N_A, N_A) 
	const Eigen::SparseMatrix<double>& get_M_nu();
	// Compute M_eps of size (N_L, N_L)
	const Eigen::SparseMatrix<double>& get_M_eps();
	// Compute C_L^A of size (N_A, N_L)
	const Eigen::SparseMatrix<double>& get_C_L_A();
	// Compute tC_L^A of size (N_tAo, N_tLo)
	const Eigen::SparseMatrix<double>& get_tC_Lo_Ao();
	// Compute tC_pL^A of size (N_tA, N_ptL)
	const Eigen::SparseMatrix<double>& get_tC_pL_A();
	// Compute G_P^pL of size (N_pL, N_pP) 
	const Eigen::SparseMatrix<double>& get_G_pP_pL();
	// Compute S_pP^{+-} of size (N_pP, N_pP)
	const Eigen::SparseMatrix<double>& get_S_pP_pm();
	// Compute tC_pL^AI of size (N_tAI, N_tpL)
	const Eigen::SparseMatrix<double>& get_tC_pL_AI();
	// Compute Q_LopP^L of size (N_L, N_Lo + N_pP)
	const Eigen::SparseMatrix<double>& get_Q_LopP_L();

private:
	void init();

	Eigen::SparseMatrix<double> M_nu;
	Eigen::SparseMatrix<double> M_eps;
	Eigen::SparseMatrix<double> C_L_A;
	Eigen::SparseMatrix<double> tC_Lo_Ao;
	Eigen::SparseMatrix<double> tC_pL_A;
	Eigen::SparseMatrix<double> G_pP_pL;
	Eigen::SparseMatrix<double> S_pP_pm;
	Eigen::SparseMatrix<double> tC_pL_AI;
	Eigen::SparseMatrix<double> Q_LopP_L;

	// index of h component ---> index of dual boundary edge
	const int ph2dpL(const int ph_idx)  const {return ph_idx  + (dual->getNumberOfEdges() - dual->facet_counts.nE_boundary);};
	// index of dual boundary edge ---> index of h component
	const int dpL2ph(const int dpL_idx) const {return dpL_idx - (dual->getNumberOfEdges() - dual->facet_counts.nE_boundary);};

	// index of phi component ---> index of primal boundary vertex 
	const int phi2ppP(const int phi_idx) const {return phi_idx + primal->facet_counts.nV_boundary;};
	// index of primal boundary vertex ---> index of phi component
	const int ppP2phi(const int ppP_idx) const {return ppP_idx - primal->facet_counts.nV_boundary;};

	// index of boundary e component ---> index of primal boundary edge
	const int pe2ppL(const int pe_idx)  const {return pe_idx  + primal->facet_counts.nE_interior;};
	// index of primal boundary edge ---> index of boundary e component
	const int ppL2pe(const int ppL_idx) const {return ppL_idx - primal->facet_counts.nE_interior;};


};

