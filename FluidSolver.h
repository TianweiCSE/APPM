#pragma once

#include <Eigen/Dense>

#include "DualMesh.h"
#include "Numerics.h"
#include "H5Writer.h"


class FluidSolver
{
	public:
		FluidSolver();
		FluidSolver(const DualMesh * mesh);
		~FluidSolver();

		// Stepping one step without considering electromagnetic field.
		const double timeStepping();

		virtual void writeSnapshot(H5Writer & writer) const;

		// virtual const std::string getXdmfOutput(const int iteration) const;

		virtual void applyInitialCondition();

	protected:
		const int U2cell(const int U_idx) const {return U2cell_map[U_idx];};
		const int cell2U(const int c_idx) const {return cell2U_map[c_idx];};
		const int F2face(const int F_idx) const {return F2face_map[F_idx];};
		const int face2F(const int f_idx) const {return face2F_map[f_idx];};

		bool isWriteStates = false;
		const DualMesh* mesh = nullptr;
		int nFaces;
		int nCells;

		const double updateFlux();
		void updateRateOfChange();

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

	private:
		Eigen::MatrixXd U;
		Eigen::MatrixXd F;
		Eigen::MatrixXd rate_of_change;

		// ----------- The index mapping might be realized by sparse vector

		// U component index    ---> dual cell index 
		Eigen::VectorXi U2cell_map;
		// dual cell index      ---> U component index
		Eigen::VectorXi cell2U_map;
		// Flux component index ---> dual face index
		Eigen::VectorXi F2face_map;
		// dual face index      ---> Flux component index
		Eigen::VectorXi face2F_map;

		const double gamma = 1.4;
		const double vareps2 = 1.0; 
};

