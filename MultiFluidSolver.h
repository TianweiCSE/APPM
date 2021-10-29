#pragma once
#include "FluidSolver.h"
class MultiFluidSolver :
	public FluidSolver
{
public:
	MultiFluidSolver();
	MultiFluidSolver(const DualMesh * dualMesh);
	~MultiFluidSolver();

	void writeStates(H5Writer & writer) const override;
	const std::string getXdmfOutput(const int iteration) const override;
	/// Set the initial condition: SOD shock tube
	void init() override;

	Eigen::VectorXd getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc) override;

	const int getNfluids() const;

private:
	std::vector<std::string> particleName;
	std::vector<double> particleMasses;
	std::vector<double> particleCharges;

	// Ratio of specific heats 
	const double gamma = 1.4; 

	/** Get wavespeed from 1d-state vector. */
	const double getWavespeed(const Eigen::Vector3d & q) const;
	const Eigen::Vector3d getFlux(const Eigen::Vector3d & q) const;
	const Eigen::Vector3d getSingleFluidState1d(const Eigen::VectorXd & q, const Eigen::VectorXd & fn, const int fluidIdx) const;
};

