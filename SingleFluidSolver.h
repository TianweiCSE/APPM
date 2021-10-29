#pragma once
#include "FluidSolver.h"

class SingleFluidSolver :
	public FluidSolver
{
public:
	SingleFluidSolver();
	SingleFluidSolver(const DualMesh * dualMesh);
	~SingleFluidSolver();

	const std::string getXdmfOutput(const int iteration) const override;
	void writeStates(H5Writer & writer) const override;

	/// Set the initial condition: SOD shock tube
	void init() override;

protected:
	/** 
	* Rusanov flux, computed by the adjacient states, face normal, distance to cell centers.
	* @return local timestep size due to maximum wave speed.
	*/
	Eigen::VectorXd getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc) override;
	Eigen::Vector3d getFlux(const Eigen::Vector3d & q);
	double maxWaveSpeed(const Eigen::Vector3d & q);

private:
	const double massRatio = 1.0;
	const double gamma = 1.4;

	struct PrimitiveState {
		double p = 0;
		double rho = 0;
		Eigen::Vector3d u = Eigen::Vector3d::Zero();
	};
	PrimitiveState getPrimitiveState(const Eigen::VectorXd & q) const;
	Eigen::VectorXd getConservationState(const PrimitiveState & primitive) const;

	friend std::ostream & operator<<(std::ostream & os, const PrimitiveState & obj);

};

