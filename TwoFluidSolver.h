#pragma once
#include "FluidSolver.h"
class TwoFluidSolver :
	public FluidSolver
{
public:
	TwoFluidSolver();
	TwoFluidSolver(const DualMesh * dualMesh);
	~TwoFluidSolver();

	void writeStates(H5Writer & writer) const override;
	const std::string getXdmfOutput(const int iteration) const override;
	/// Set the initial condition: SOD shock tube
	void init() override;

protected:
	Eigen::VectorXd getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc) override;
	Eigen::Vector3d getFlux(const Eigen::Vector3d & q);
	//const double maxWaveSpeed(const Eigen::VectorXd & q, const Eigen::Vector3d & fn);
	const double maxWaveSpeed(const Eigen::Vector3d & q) const;

private:
	const double gamma = 1.4;
	const double massRatio_a = 1.;
	const double massRatio_b = 1. / 1836.;

	const int stateVectorLength = 10;

	struct PrimitiveState {
		double p_a = 0.;
		double n_a = 0.;
		Eigen::Vector3d u_a = Eigen::Vector3d::Zero();
		double p_b = 0.;
		double n_b = 0.;
		Eigen::Vector3d u_b = Eigen::Vector3d::Zero();
	};

	PrimitiveState getPrimitiveState(const Eigen::VectorXd & q) const;
	Eigen::VectorXd getConservationState(const PrimitiveState & primitive) const;

	friend std::ostream & operator<<(std::ostream & os, const PrimitiveState & obj);

	const Eigen::VectorXd getInitStateLeft() const;
	const Eigen::VectorXd getInitStateRight() const;
	const Eigen::VectorXd getState1d(const Eigen::VectorXd & q_3d, const Eigen::Vector3d & n);
};

