#include "SingleFluidSolver.h"



SingleFluidSolver::SingleFluidSolver()
{
}

SingleFluidSolver::SingleFluidSolver(const DualMesh * dualMesh)
	: FluidSolver(dualMesh)
{
}

SingleFluidSolver::~SingleFluidSolver()
{
}

const std::string SingleFluidSolver::getXdmfOutput(const int iteration) const
{
	std::stringstream ss;

	ss << "<Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/density" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/pressure" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/velocity" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	return ss.str();
}

void SingleFluidSolver::writeStates(H5Writer & writer) const
{
	const int nDualCells = mesh->getNumberOfCells();

	if (isWriteStates) {
		Eigen::VectorXd q_density = fluidStates.row(0);
		Eigen::MatrixXd q_momentum = fluidStates.middleRows(1, 3);
		Eigen::VectorXd q_energy = fluidStates.row(4);
		writer.writeData(q_density, "/qDensity");
		writer.writeData(q_momentum, "/qMomentum");
		writer.writeData(q_energy, "/qEnergy");
	}

	// Fluid states in primitive variables
	Eigen::VectorXd pressure(nDualCells);
	Eigen::VectorXd density(nDualCells);
	Eigen::Matrix3Xd velocity(3, nDualCells);
	for (int i = 0; i < nDualCells; i++) {
		Eigen::VectorXd state = fluidStates.col(i);
		PrimitiveState primitiveState = getPrimitiveState(state);// FluidState(fluidStates.col(i)).getPrimitiveState();
		pressure(i) = primitiveState.p;
		density(i) = primitiveState.rho;
		velocity.col(i) = primitiveState.u;
	}
	writer.writeData(pressure, "/pressure");
	writer.writeData(density, "/density");
	writer.writeData(velocity, "/velocity");

}

void SingleFluidSolver::init()
{
	std::cout << "SingleFluidSolver::init()" << std::endl;

	const int nCells = mesh->getNumberOfCells();
	assert(nCells > 0);

	fluidStates = Eigen::MatrixXd::Zero(5, nCells);
	fluidFluxes = Eigen::MatrixXd::Zero(5, nCells);

	const double rho_L = 1.0;
	const double p_L = 1.0;
	const double rho_R = 0.125;
	const double p_R = 0.1;
	const Eigen::Vector3d uZero = Eigen::Vector3d::Zero();

	PrimitiveState primitiveL;
	primitiveL.rho = rho_L;
	primitiveL.p = p_L;
	primitiveL.u = uZero;

	PrimitiveState primitiveR;
	primitiveR.rho = rho_R;
	primitiveR.p = p_R;
	primitiveR.u = uZero;

	std::cout << "pL: " << std::endl << primitiveL << std::endl;
	std::cout << "pR: " << std::endl << primitiveR << std::endl;

	const Eigen::VectorXd qL = getConservationState(primitiveL);
	const Eigen::VectorXd qR = getConservationState(primitiveR);

	std::cout << "qL: " << qL.transpose() << std::endl;
	std::cout << "qR: " << qR.transpose() << std::endl;
	std::cout << std::endl;

	const PrimitiveState pStateL = getPrimitiveState(qL);
	const PrimitiveState pStateR = getPrimitiveState(qR);
	std::cout << "pStateL: " << std::endl << pStateL << std::endl;
	std::cout << "pStateR: " << std::endl << pStateR << std::endl;

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = mesh->getCell(i);
		const double idxC = cell->getIndex();
		const Eigen::Vector3d cellCenter = cell->getCenter();

		fluidStates.col(i) = (cellCenter(2) < 0) ? qL : qR;
	}

	std::cout << "Test for fluid flux: " << std::endl;
	double dt_loc = 0;
	const Eigen::Vector3d fn = Eigen::Vector3d(0, 0, 1);
	Eigen::VectorXd fluidFlux = getRusanovFlux(qL, qR, fn, 1, dt_loc);
	std::cout << "Flux: " << fluidFlux.transpose() << std::endl;
	std::cout << "dt_loc: " << dt_loc << std::endl;

}

Eigen::VectorXd SingleFluidSolver::getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc)
{
	assert(qL.size() == qR.size());
	assert(qL.size() == 5);
	assert(qR.size() == 5);
	assert(dx > 0);

	Eigen::Vector3d qL_1d;
	qL_1d(0) = qL(0);
	qL_1d(1) = qL.segment(1, 3).dot(fn);
	qL_1d(2) = qL(4);
	const double sL = maxWaveSpeed(qL_1d);
	const Eigen::Vector3d fL_1d = getFlux(qL_1d);

	Eigen::Vector3d qR_1d;
	qR_1d(0) = qR(0);
	qR_1d(1) = qR.segment(1, 3).dot(fn);
	qR_1d(2) = qR(4);
	const double sR = maxWaveSpeed(qR_1d);
	const double s = std::max(sL, sR);
	const Eigen::Vector3d fR_1d = getFlux(qR_1d);

	dt_loc = dx / s;

	const Eigen::Vector3d flux_1d = 0.5 * (fL_1d + fR_1d) - 0.5 * s * (qR_1d - qL_1d);

	Eigen::VectorXd flux(5);
	flux(0) = flux_1d(0);
	flux.segment(1, 3) = flux_1d(1) * fn;
	flux(4) = flux_1d(2);
	return flux;

	/**
	 * What I think is the right way:
	 * 
	 * 	PrimitiveState ps_L = getPrimitiveState(qL)
	 * 	int u_normal_L = ps_L.u.dot(fn)
	 * 	PrimitiveState ps_R = getPrimitiveState(qR)
	 * 	int u_normal_R = ps_R.u.dot(fn)
	 * 
	 * 	Eigen::Vector5d fL, fR;
	 * 
	 * 	fL(0) = ps_L.rho * u_normal_L;
	 * 	fL.segment(1,3) = ps_L.rho * u_normal_L * u + ps_L.p * fn;
	 * 	fL(4) = u_normal_L * (qL(4) + ps_L.p);
	 * 
	 * 	fR(0) = ps_R.rho * u_normal_R;
	 * 	fR.segment(1,3) = ps_R.rho * u_normal_R * u + ps_R.p * fn;
	 * 	fR(4) = u_normal_R * (qR(4) + ps_R.p);
	 * 
	 * 	Eigen::Vector5d flux = 0.5 * (fL + fR) - 0.5 * s * (qR - qL);
	 * 
	 *  return flux
	 */
}

Eigen::Vector3d SingleFluidSolver::getFlux(const Eigen::Vector3d & q) {
	Eigen::Vector3d flux;
	flux(0) = q(1);
	flux(1) = (gamma - 1) * q(2) + 0.5 * (3. - gamma) * pow(q(1), 2) / q(0);
	flux(2) = gamma * q(2) * q(1) / q(0) + 0.5 * (1. - gamma) * pow(q(1), 3) / pow(q(0), 2);
	return flux;
}

double SingleFluidSolver::maxWaveSpeed(const Eigen::Vector3d & q)
{
	const double u = q(1) / q(0);
	const double s2 = gamma * (gamma - 1) * ( q(2) / q(0) - 0.5 * pow(q(1) / q(0), 2.0) );
	assert(s2 > 0);
	const double s = sqrt(s2);
	return (u * Eigen::Vector3d::Ones() + s * Eigen::Vector3d(-1, 0, 1)).cwiseAbs().maxCoeff();
}

Eigen::VectorXd SingleFluidSolver::getConservationState(const PrimitiveState & primitive) const
{
	const double rho = primitive.rho;
	const Eigen::Vector3d u = primitive.u;
	const double p = primitive.p;
	const double rho_etot = p / (gamma - 1) + 0.5 * rho * u.squaredNorm();

	Eigen::VectorXd q(5);
	q(0) = primitive.rho;
	q.segment(1, 3) = rho * u;
	q(4) = rho_etot;
	return q;
}

SingleFluidSolver::PrimitiveState SingleFluidSolver::getPrimitiveState(const Eigen::VectorXd & q) const
{
	const Eigen::Vector3d u = q.segment(1, 3) / q(0);
	const double rho = q(0);
	const double rho_etot = q(4);
	const double p = (gamma - 1) * (rho_etot - 0.5 * rho * u.squaredNorm());
	
	PrimitiveState primitive;
	primitive.rho = rho;
	primitive.p = p;
	primitive.u = u;
	return primitive;
}

std::ostream & operator<<(std::ostream & os, const SingleFluidSolver::PrimitiveState & obj)
{
	os << "p: " << obj.p << std::endl;
	os << "rho: " << obj.rho << std::endl;
	os << "u: " << obj.u.transpose() << std::endl;
	return os;
}
