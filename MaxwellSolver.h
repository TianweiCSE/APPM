#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "H5Writer.h"
#include "FluidSolver.h"
#include "Interpolator.h"
#include "EigenAuxiliaries.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;
class AppmSolver;

/**
 * @brief A class that is for dealing with the Maxwell equation. 
 * 	      Information about the (generalized) Ohm's law (J = M_sigma * E + J_aux) must be provided.
 * 		  The work flow at each time step is:
 * 				- assemble the (sparse) linear system;
 * 				- solve the linear system and get the new E-field (as well as some auxiliary variables)
 *     			- evolve the B-field
 */
class MaxwellSolver
{
public:
	struct MaxwellParams {
		double lambdaSquare = 1.0;  //< nondimensionalized Debye length
	} parameters;

	MaxwellSolver();
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator);
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, const Interpolator* interpolator, MaxwellParams & param);
	~MaxwellSolver();

	// Store the data to .h5 file
	void writeSnapshot(H5Writer & writer) const;
	/**
	 * @brief Assemble the linear system and solve it. 
	 * 
	 * The variables to be solved simultaneouly are [e^o, phi, h^partial, d^\partial]
	 * 
	 * @param time current time
	 * @param dt time step size
	 * @param M_sigma maxtrix in Ohm's law
	 * @param j_aux vector in Ohm's law
	 */
	void solveLinearSystem(const double time, const double dt, Eigen::SparseMatrix<double>&& M_sigma, Eigen::VectorXd&& j_aux);
	// Evolve magnetic flux vector <b> through (4.31)
	void evolveMagneticFlux(const double dt);

	// Return interpolated E-field at dual cell center. Each row vector is indexed by dual cell index.
	// The row corresponding to non-fluid cell is not meaningful.
	Eigen::MatrixXd updateInterpolated_E();
	// Return interpolated B-field at dual cell center. Each row vector is indexed by dual cell index.
	// The row corresponding to non-fluid cell is not meaningful.
	Eigen::MatrixXd updateInterpolated_B();
	Eigen::MatrixXd getInterpolated_E() const;
	Eigen::MatrixXd getInterpolated_B() const;   

	// Get the electric potential assigned to the anode.
	double getPotential(const double t) const {return -1.0;};

	// Assign initital conditons to electromagnetic variables
	void applyInitialCondition();
	void applyInitialCondition(const std::string h5_file);

	Eigen::VectorXd getNorms() const;

protected:
	const PrimalMesh* primal = nullptr;
	const DualMesh* dual = nullptr;
	const Interpolator* interpolator = nullptr;

	/* These are variables to be solved in linear system */
	Eigen::VectorXd eo;  //< E-field integral at interior primal edge
	Eigen::VectorXd phi; //< electric potenial at boundary primal vertex
	Eigen::VectorXd hp;	 //< H-field integral at boundary dual edge  
	Eigen::VectorXd dp;	 //< D-field integral at boundary dual face

	/* These are auxiliary variables */
	Eigen::VectorXd e;   //< E-field integral at priaml edge
	Eigen::VectorXd b;   //< B-field integral at primal face
	Eigen::VectorXd j;   //< current at dual face
	Eigen::MatrixXd E;	 //< interpolated E-field defined at dual cell center
	Eigen::MatrixXd B;   //< interpolated B-field defined at dual cell center

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

	Eigen::SparseMatrix<double> M_nu;
	Eigen::SparseMatrix<double> M_eps;
	Eigen::SparseMatrix<double> C_L_A;
	Eigen::SparseMatrix<double> tC_Lo_Ao;
	Eigen::SparseMatrix<double> tC_pL_A;
	Eigen::SparseMatrix<double> G_pP_pL;
	Eigen::SparseMatrix<double> S_pP_pm;
	Eigen::SparseMatrix<double> tC_pL_AI;
	Eigen::SparseMatrix<double> Q_LopP_L;

	// index of boundary h component ---> index of dual boundary edge
	const int ph2dpL(const int ph_idx)  const {return ph_idx  + (dual->getNumberOfEdges() - dual->facet_counts.nE_boundary);};
	// index of dual boundary edge ---> index of boundary h component
	const int dpL2ph(const int dpL_idx) const {return dpL_idx - (dual->getNumberOfEdges() - dual->facet_counts.nE_boundary);};

	// index of boundary d component ---> index of dual boundary face
	const int pd2dpA(const int pd_idx) const {return pd_idx + (dual->getNumberOfFaces() - dual->facet_counts.nF_boundary);};
	// index of dual boundary face ---> index of boundary d component
	const int dpA2pd(const int dpA_idx) const {return dpA_idx - (dual->getNumberOfFaces() - dual->facet_counts.nF_boundary);};

	// index of phi component ---> index of primal boundary vertex 
	const int phi2ppP(const int phi_idx) const {return phi_idx + primal->facet_counts.nV_interior;};
	// index of primal boundary vertex ---> index of phi component
	const int ppP2phi(const int ppP_idx) const {return ppP_idx - primal->facet_counts.nV_interior;};

	// index of boundary e component ---> index of primal boundary edge
	const int pe2ppL(const int pe_idx)  const {return pe_idx  + primal->facet_counts.nE_interior;};
	// index of primal boundary edge ---> index of boundary e component
	const int ppL2pe(const int ppL_idx) const {return ppL_idx - primal->facet_counts.nE_interior;};
	
	friend class AppmSolver;
};

