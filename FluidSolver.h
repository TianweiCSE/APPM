#pragma once

#include <Eigen/Dense>
#include "DualMesh.h"
#include "Numerics.h"
#include "H5Writer.h"
#include "Interpolator.h"

// TODO:  - Leave all the basic fluid solver to class <BasicFluidSolver>. 
//        - Rename this class as <SingleFluidSolver>, and keep all the members relevent to plasma modelling. 

class TwoFluidSolver; 
class AppmSolver;

/**
 * @brief A class for solving Euler equation by Rusanov scheme and Euler time-stepping.
 * 
 * The typical work flow is:
 * 		- compute numberical flux at each face
 * 		- compute rhs(source term)
 * 		- compute the rate of change
 * 		- do time stepping by Euler scheme
 */
class FluidSolver
{
	public:
		FluidSolver() = delete;
		FluidSolver(const DualMesh* mesh);
		FluidSolver(const DualMesh* mesh, const double gamma, const double mass, const double charge, const std::string name);
		~FluidSolver();

		/**
		 * @brief Update the fluxes explicitly, as is done in normal FVM, at each fluid faces.
		 * 
		 * This function should be called at the beginning of each iteration to get the time step.
		 * 
		 * @return time step size that is restricted by CFL condition
		 */
		const double updateFluxExplicit();
		/**
		 * @brief  Update the semi-implicit fluxes as (4.23). 
		 * 
		 * Since the explicit fluxes are already computed in advance, modification only happens for the advection
		 * part of mass flux. Thus, the strategy is to extract the old and new velocities, and add the difference to 
		 * the existing mass flux, with momentum and energy fluxed unchanged.
		 * 
		 */
		void updateMassFluxImplicit();
		
		void updateMomentum(const double dt, const Eigen::MatrixXd &E);
		void updateMomentum(const double dt, const Eigen::MatrixXd &E, const double alpha, const FluidSolver* anotherSpecies);

		// Evolve one step without considering source terms.
		void timeStepping(const double dt);
		/**
		 * @brief Evolve one step with source terms
		 * 
		 * @param dt time step size
		 * @param E new E-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 * @param B old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 */
		void timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);

		// This function should NOT be used anymore, since in this way one of the fluid is updated first!
		void timeStepping(const double dt, 
						  const Eigen::MatrixXd &E, 
						  const Eigen::MatrixXd &B, 
						  const double alpha,
						  const FluidSolver* anotherSpecies);

		/**
		 * @brief Write density, velocity, pressure into .h5 file. 
		 *	- the output matrices are adjusted to be consistent to the facets indices.
		 *  - the data at solid cell is set zero.
		 * 
		 * @param writer h5writer object
		 */
		void writeSnapshot(H5Writer &writer) const;

		// Assign initial conditions to fluid variables
		void applyInitialCondition();
		void applyInitialCondition(const std::string h5_file);

		Eigen::VectorXd getNorms() const;

	private:

		const int U2cell(const int U_idx) const {return U2cell_map[U_idx];};
		const int cell2U(const int c_idx) const {return cell2U_map[c_idx];};
		const int F2face(const int F_idx) const {return F2face_map[F_idx];};
		const int face2F(const int f_idx) const {return face2F_map[f_idx];};

		// Return vector of number density whose entry is indexed by cell index. Entries of solid cells are set ZERO.
		Eigen::VectorXd getExtended_n() const;
		Eigen::VectorXd getExtended_s() const;

		const DualMesh* mesh = nullptr;
		const int nFaces;
		const int nCells;

		/**
		 * @brief add Lorentz force to <RHS> 
		 * 
		 * @param E New E-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 * @param B Old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 */
		void applyLorentzForce(const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);
		/**
		 * @brief add friction term to <RHS>
		 * 
		 * @param anotherSpecies pointer to the solver of another spieces
		 * @param alpha friction coefficient
		 */
		void applyFrictionTerm(const FluidSolver* anotherSpecies, const double alpha);
		/**
		 * @brief Compute the rate of change of fluid conservative variables.
		 * 
		 * @param with_rhs True if the source term needs to be considered. 
		 */
		void updateRateOfChange(const bool with_rhs);

		const double updateFluxInterior(const int faceIdx);
		const double updateFluxOpening (const int faceIdx);
		const double updateFluxWall    (const int faceIdx);
		const double updateFluxMixed   (const int faceIdx);

		Eigen::VectorXd Flux(const Eigen::VectorXd &q, const Eigen::Vector3d &fn) const;

		Eigen::VectorXd RusanovFlux(const Eigen::VectorXd & qL, 
									const Eigen::VectorXd & qR,
									const Eigen::Vector3d & fn) const;

		Eigen::VectorXd conservative2primitive(const Eigen::VectorXd &q_cons) const;
		Eigen::VectorXd primitive2conservative(const Eigen::VectorXd &q_prim) const;
		const double speedOfSound(const Eigen::VectorXd &q_cons) const;
		const double maxWaveSpeed(const Eigen::VectorXd &q_cons, const Eigen::Vector3d &normal) const;

		// Check the positivity of number density and pressure
		bool isValidState() const;
		// Check A and D
		void check_A_and_D(const Tensor3& A, const Eigen::MatrixXd& D) const;
		void check_eta() const;
		void check_eta2(const Eigen::MatrixXd &B, const double dt);
		void check_updatedMomentum() const;
		void check_updatedMomentum2(const double dt, const Eigen::MatrixXd &E, const FluidSolver* another) const;
		void check_updatedMomentum3(const double dt, const Eigen::MatrixXd &E, const FluidSolver* another) const;

		// Update eta matrix defined in (4.37)
		// Note: The size of rows is extended to number of dual cells, with soild entry being zero.
		//       The indexing of row is consisent with indexing of dual cells, not U.
		//       This is for the convenience of later usage.
		void update_eta(const double dt, const Eigen::MatrixXd &B) const;

		/**
		 * @brief Compute mu vector defined in (4.40)
		 * 
		 * @param dt time step size
		 * @param B Old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 * @param A Tensor of rank 3. See definition in (4.39) 
		 * @param D Matrix of rank 2. See definition in (4.39)
		 * @return mu vector. Indexed in face indexing.
		 */
		Eigen::VectorXd get_mu(const double dt, 
		                       const Eigen::MatrixXd &B, 
							   const Tensor3& A, 
							   const Eigen::SparseMatrix<double>& D) const;
		
		/**
		 * @brief Compute mu vector defined in (4.40)
		 * 
		 * @param dt time step size
		 * @param B Old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
		 * @param A Tensor of rank 3. See definition in (4.39) 
		 * @param D Matrix of rank 2. See definition in (4.39)
		 * @param alpha friction coefficient
		 * @param anotherSpecies pointer to the solver of another spieces
		 * @return mu vector. Indexed in face indexing.
		 */
		Eigen::VectorXd get_mu(const double dt, 
		                       const Eigen::MatrixXd &B, 
							   const Tensor3& A, 
							   const Eigen::SparseMatrix<double>& D,
							   const double alpha,
							   const FluidSolver* anotherSpecies) const;

		/**
		 * @brief Compute T matrix defined in (4.40)
		 * 
		 * @param dt time step size
		 * @param A Tensor of rank 3. See definition in (4.39)
		 * @param R Tensor for E-field interpolation
		 * @return T matrix
		 */
		Eigen::SparseMatrix<double> get_T(const double dt,
										  const Tensor3& A,
										  const Tensor3& R) const;
		
		/**
		 * @brief Compute T matrix defined in (5.18)
		 * 
		 * @param dt time step size
		 * @param A Tensor of rank 3. See definition in (4.39)
		 * @param R Tensor for E-field interpolation
		 * @param alpha friction coefficient
		 * @param anotherSpecies pointer to the solver of another spieces
		 * @return T matrix
		 */
		Eigen::SparseMatrix<double> get_T(const double dt,
										  const Tensor3& A,
										  const Tensor3& R,
										  const double alpha,
										  const FluidSolver* anotherSpecies) const;

		Eigen::MatrixXd U;   				//< conservative variable at each fluid cell
		Eigen::MatrixXd F;                  //< flux at each fluid face
		Eigen::MatrixXd rhs;				//< rhs at each fluid cell
		Eigen::MatrixXd rate_of_change;		//< rate of change of conservative variable
		Eigen::VectorXd S;                  //< artificial diffusive coefficient
		Eigen::MatrixXd updatedMomentum;    //< an temporay variable to store new momentum 

		mutable Eigen::MatrixXd eta;		//< eta defined in (4.37)
		
		// U component index    ---> dual cell index 
		Eigen::VectorXi U2cell_map;
		// dual cell index      ---> U component index
		Eigen::VectorXi cell2U_map;
		// Flux component index ---> dual face index
		Eigen::VectorXi F2face_map;
		// dual face index      ---> Flux component index
		Eigen::VectorXi face2F_map;

		const double gamma   = 1.4;
		const double vareps2 = 1.0; 
		const double charge  = 0.0;

		const std::string name = "";

		friend class TwoFluidSolver;
		friend class AppmSolver;
};

