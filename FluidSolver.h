#pragma once

#include <Eigen/Dense>
#include "DualMesh.h"
#include "Numerics.h"
#include "H5Writer.h"
#include "Interpolator.h"

// TODO:  - Leave all the basic fluid solver to class <BasicFluidSolver>. 
//        - Rename this class as <SingleFluidSolver>, and keep all the members relevent to plasma modelling. 

class TwoFluidSolver; 

class FluidSolver
{
	public:
		FluidSolver();
		FluidSolver(const DualMesh* mesh);
		FluidSolver(const DualMesh* mesh, const double gamma, const double mass, const double charge, const std::string name);
		~FluidSolver();

		// Update the fluxes explicitly at each fluid faces. Then return dt which is restricted by CFL condition.
		const double updateFluxExplicit();
		// Update the fluxes implicitly as (4.23). 
		void updateFluxImplicit(const Eigen::MatrixXd &E);
		// Evolve one step without considering electromagnetic field.
		void timeStepping(const double dt);
		// Evolve one step with electromagnetic force
		void timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);
		// Write density, velocity, pressure into .h5 file. 
		// Note: - the output matrices are adjusted to be consistent to the facets indices.
		//       - the data at solid cell is set zero.
		void writeSnapshot(H5Writer &writer) const;

		// virtual const std::string getXdmfOutput(const int iteration) const;

		void applyInitialCondition();

	private:

		const int U2cell(const int U_idx) const {return U2cell_map[U_idx];};
		const int cell2U(const int c_idx) const {return cell2U_map[c_idx];};
		const int F2face(const int F_idx) const {return F2face_map[F_idx];};
		const int face2F(const int f_idx) const {return face2F_map[f_idx];};

		bool isWriteStates = false;
		const DualMesh* mesh = nullptr;
		int nFaces;
		int nCells;

		void updateRHS(const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);
		void updateRateOfChange(const bool with_rhs);

		const double updateFluxInterior(const int faceIdx);
		const double updateFluxOpening (const int faceIdx);
		const double updateFluxWall    (const int faceIdx);

		Eigen::VectorXd Flux(const Eigen::VectorXd &q, const Eigen::Vector3d &fn) const;

		Eigen::VectorXd RusanovFlux(const Eigen::VectorXd & qL, 
									const Eigen::VectorXd & qR,
									const Eigen::Vector3d & fn) const;

		Eigen::VectorXd conservative2primitive(const Eigen::VectorXd &q_cons) const;
		Eigen::VectorXd primitive2conservative(const Eigen::VectorXd &q_prim) const;
		const double speedOfSound(const Eigen::VectorXd &q_cons) const;
		const double maxWaveSpeed(const Eigen::VectorXd &q_cons, const Eigen::Vector3d &normal) const;

		// Update eta matrix defined in (4.37)
		// Note: The size of rows is extended to number of dual cells, with soild entry being zero.
		//       The indexing of row is consisent with indexing of dual cells, not U.
		//       This is for the convenience of later usage.
		void update_eta(const double dt, const Eigen::MatrixXd& B) const;

		// Get mu matrix defined in (4.40)
		Eigen::VectorXd&& get_mu(const double dt, 
		                         const Eigen::MatrixXd& B, 
								 const Tensor3& A, 
								 const Eigen::SparseMatrix<double>& D) const;
		// Get T tensor defined in (4.40)
		Eigen::SparseMatrix<double>&& get_T(const double dt,
										    const Tensor3& A,
											const Tensor3& R) const;

		Eigen::MatrixXd U;   				//< conservative variable at each fluid cell
		Eigen::MatrixXd F;                  //< flux at each fluid face
		Eigen::MatrixXd rhs;				//< rhs at each fluid cell
		Eigen::MatrixXd rate_of_change;		//< rate of change of conservative variable
		mutable Eigen::MatrixXd eta;		//< eta defined in (4.37)
		

		// ----------- The index mapping might be realized by sparse vector

		// U component index    ---> dual cell index 
		Eigen::VectorXi U2cell_map;
		// dual cell index      ---> U component index
		Eigen::VectorXi cell2U_map;
		// Flux component index ---> dual face index
		Eigen::VectorXi F2face_map;
		// dual face index      ---> Flux component index
		Eigen::VectorXi face2F_map;

		double gamma   = 1.4;
		double vareps2 = 1.0; 
		double charge  = 0.0;

		std::string name = "";

		friend class TwoFluidSolver;
};

