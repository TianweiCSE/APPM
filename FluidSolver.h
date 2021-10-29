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

	/// Stepping one step without considering electromagnetic field.
	const double updateFluidState();

	virtual void writeStates(H5Writer & writer) const = 0;

	virtual const std::string getXdmfOutput(const int iteration) const = 0;

	virtual void init() = 0;

protected:
	bool isWriteStates = false;
	const DualMesh * mesh = nullptr;

	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidFluxes;

	Eigen::MatrixXd fluidStateUpdate;
	Eigen::VectorXd dt_local;

	void updateFaceFluxInterior(const int faceIdx);
	void updateFaceFluxOpening(const int faceIdx);
	void updateFaceFluxWall(const int faceIdx);


	virtual Eigen::VectorXd getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc) = 0;



private:
};

