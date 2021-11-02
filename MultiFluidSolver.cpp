#include "MultiFluidSolver.h"



MultiFluidSolver::MultiFluidSolver()
	: FluidSolver()
{
}

MultiFluidSolver::MultiFluidSolver(const DualMesh * dualMesh)
	: FluidSolver(dualMesh)
{
}

MultiFluidSolver::~MultiFluidSolver()
{
}

void MultiFluidSolver::writeStates(H5Writer & writer) const
{
	const int nFluids = this->getNfluids();
	const int nCells = mesh->getNumberOfCells();

	for (int k = 0; k < nFluids; k++) {
		std::stringstream ss;
		ss << "/fluid" << k << "-";
		const std::string fluidTag = ss.str();
		const std::string pressureTag = fluidTag + "pressure";
		const std::string velocityTag = fluidTag + "velocity";
		const std::string densityTag  = fluidTag + "density";
		const std::string numberDensityTag = fluidTag + "numberDensity";

		Eigen::VectorXd numberDensity = fluidStates.row(5 * k);
		Eigen::VectorXd density(nCells);
		Eigen::MatrixXd velocity(3, nCells);
		Eigen::VectorXd pressure(nCells);
		const double epsilon2 = particleMasses[k];
		for (int i = 0; i < nCells; i++) {
			const Eigen::VectorXd state = fluidStates.block(5 * k, i, 5, 1);

			const double n = state(0);
			const double rho = epsilon2 * n;
			const Eigen::Vector3d u = epsilon2 * state.segment(1,3) / n;
			const double p = epsilon2 * (gamma - 1) * (state(4) - 0.5 * n * u.squaredNorm());

			density(i) = rho;
			velocity.col(i) = u;
			pressure(i) = p; 
		}

		writer.writeDoubleVector(density, densityTag);
		writer.writeDoubleVector(numberDensity, numberDensityTag);
		writer.writeDoubleVector(pressure, pressureTag);
		writer.writeDoubleMatrix(velocity, velocityTag);
	}
}

const std::string MultiFluidSolver::getXdmfOutput(const int iteration) const
{
	std::stringstream ss;
	const int nFluids = this->getNfluids();
	for (int k = 0; k < nFluids; k++) {
		std::stringstream ss;
		ss << "" << k;
		const std::string fluidName = ss.str();

		ss << "<Attribute Name=\"Density " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << "fluid" << k << "-density" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Number Density " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << "fluid" << k << "-numberDensity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Pressure " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << "fluid" << k << "-pressure" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;


		ss << "<Attribute Name=\"Velocity " << fluidName << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << mesh->getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << "fluid" << k << "-velocity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}
	return ss.str();
}

void MultiFluidSolver::init()
{
	// Read particle parameters from file
	const std::string filename = "particleParameters.txt";
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File not opened: " << filename << std::endl;
		exit(-1);
	}
	std::string line;
	while (std::getline(file, line)) {
		if (line.empty() || line.substr(0, 1) == "#") {
			continue;
		}
		std::string fluidName;
		double mass;
		int charge;
		char delimiter;
		std::stringstream ss(line);
		ss >> fluidName;
		ss >> mass >> delimiter;
		ss >> charge;
		particleName.push_back(fluidName.substr(0, fluidName.size() - 1));
		particleMasses.push_back(mass);
		particleCharges.push_back(charge);
	}

	assert(particleMasses.size() == particleCharges.size()); 
	const int nFluids = particleMasses.size();
	std::cout << "Fluid parameters:" << std::endl;
	std::cout << std::endl;
	std::cout << "Fluid #: Name,\tMass,\tCharge" << std::endl;
	std::cout << "=====================" << std::endl;
	for (int i = 0; i < nFluids; i++) {
		std::cout << "Fluid " << i << ": " << particleName[i] << "\t" << particleMasses[i] << "\t" << particleCharges[i] << std::endl;
	}

	// Define initial states
	const int nCells = mesh->getNumberOfCells();
	const int fluidStateLength = 5 * nFluids;
	fluidStates = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	fluidFluxes = Eigen::MatrixXd::Zero(fluidStateLength, nCells);

	// Sod shock tube test states
	Eigen::VectorXd leftState(fluidStateLength);
	Eigen::VectorXd rightState(fluidStateLength);

	for (int k = 0; k < nFluids; k++) {
		// TODO: edit such that number density is equal among the fluids (= quasi-neutral plasma state), not the mass density.
		const double epsilon2 = particleMasses[k];

		const double pL = 1;
		const double nL = 1;
		const double rhoL = epsilon2 * nL;
		const Eigen::Vector3d uL = Eigen::Vector3d::Zero();

		const double pR = 0.1;
		const double nR = 0.125; 
		const double rhoR = epsilon2 * nR;
		const Eigen::Vector3d uR = Eigen::Vector3d::Zero();

		Eigen::VectorXd singleFluidStateLeft(5);
		singleFluidStateLeft(0) = nL;
		singleFluidStateLeft.segment(1, 3).setZero();
		singleFluidStateLeft(4) = 1./epsilon2 * (pL / (gamma-1) + 0.5 * rhoL * uL.squaredNorm());

		Eigen::VectorXd singleFluidStateRight(5);
		singleFluidStateRight(0) = nR;
		singleFluidStateRight.segment(1, 3).setZero();
		singleFluidStateRight(4) = 1./epsilon2 * (pR / (gamma - 1) + 0.5 * rhoR * uR.squaredNorm());

		leftState.segment( 5 * k, 5) = singleFluidStateLeft;
		rightState.segment(5 * k, 5) = singleFluidStateRight;
	}

	// Set states on the left and right side
	for (int i = 0; i < nCells;  i++) {
		const Cell * cell = mesh->getCell(i);
		Eigen::VectorXd cellState(fluidStateLength);

		if (cell->getFluidType() == Cell::FluidType::FLUID) {
			const Eigen::Vector3d cellCenterPos = cell->getCenter();
			cellState = (cellCenterPos(2) < 0) ? leftState : rightState;
		}
		else {
			const double a = std::nan(""); // value of Not-A-Number
			cellState.setConstant(a);
		}

		fluidStates.col(i) = cellState;
	}

	// Check for Rusanov flux
	//const Eigen::Vector3d fn = Eigen::Vector3d::UnitZ();
	//const double dx = 1;
	//double dt_loc = 0;
	//const Eigen::VectorXd flux = getRusanovFlux(leftState, rightState, fn, dx, dt_loc);
	//std::cout << "flux: " << std::endl;
	//std::cout << flux << std::endl;
	//std::cout << "dt_loc: " << dt_loc << std::endl;
}

Eigen::VectorXd MultiFluidSolver::getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc)
{
	const int stateLength = 5;
	const int nFluids = this->getNfluids();
	Eigen::VectorXd flux(stateLength * nFluids);
	assert(flux.size() == fluidFluxes.rows());
	Eigen::VectorXd maxWavespeed(nFluids);

	for (int k = 0; k < nFluids; k++) {

		const Eigen::Vector3d qL_1d = getSingleFluidState1d(qL, fn, k);
		const Eigen::Vector3d qR_1d = getSingleFluidState1d(qR, fn, k);
		const double sL = getWavespeed(qL_1d);
		const double sR = getWavespeed(qR_1d);
		const double s = std::max(sL, sR);
		maxWavespeed(k) = s;

		const Eigen::Vector3d fL_1d = getFlux(qL_1d);
		const Eigen::Vector3d fR_1d = getFlux(qR_1d);
		const Eigen::Vector3d singleFluidFlux1d = 0.5 * (fL_1d + fR_1d) - 0.5 * s * (qR_1d - qL_1d);

		Eigen::VectorXd singleFluidFlux(5);
		assert(singleFluidFlux.size() == stateLength);
		singleFluidFlux(0) = singleFluidFlux1d(0);
		singleFluidFlux.segment(1,3) = singleFluidFlux1d(1) * fn;
		singleFluidFlux(4) = singleFluidFlux1d(2);

		flux.segment(5 * k, 5) = singleFluidFlux;
	}
	const double smax = maxWavespeed.maxCoeff();
//#ifdef _DEBUG
//	std::cout << "wavespeeds: " << maxWavespeed.transpose() << std::endl;
//	std::cout << "smax: " << smax << std::endl;
//	std::cout << "ratio: " << smax * Eigen::VectorXd::Ones(nFluids).cwiseQuotient(maxWavespeed).transpose() << std::endl;
//#endif
	dt_loc = dx / smax;
	return flux;
}

const int MultiFluidSolver::getNfluids() const
{
	return particleMasses.size();
}

/**
@return s = sqrt(gamma * p / rho) = sqrt( gamma * (gamma-1) * e )
*/
const double MultiFluidSolver::getWavespeed(const Eigen::Vector3d & q) const
{
	assert(q(0) > 0);
	assert(q(2) > 0);
	const double s2 = gamma * (gamma - 1) * (q(2) / q(0) - 0.5 * pow(q(1), 2) / q(0));
	assert(s2 > 0);
	const double s = sqrt(s2);
	const double u = q(1);
	return Eigen::Vector3d(u + s, u + 0, u - s).cwiseAbs().maxCoeff();
}

const Eigen::Vector3d MultiFluidSolver::getFlux(const Eigen::Vector3d & q) const
{
	Eigen::Vector3d flux;
	flux(0) = q(1);
	flux(1) = pow(q(1), 2) / q(0) + (gamma - 1) * (q(2) - 0.5 * pow(q(1), 2) / q(0));
	flux(2) = (q(2) + (gamma - 1) * (q(2) - 0.5 * pow(q(1),2) / q(0))) * q(1) / q(0);
	return flux;
}

const Eigen::Vector3d MultiFluidSolver::getSingleFluidState1d(const Eigen::VectorXd & q, const Eigen::VectorXd & fn, const int fluidIdx) const
{
	const Eigen::VectorXd state3d = q.segment(5 * fluidIdx, 5);
	assert(state3d.size() == 5);
	Eigen::Vector3d state;
	state(0) = state3d(0);
	state(1) = state3d.segment(1, 3).dot(fn);
	state(2) = state3d(4);
	return state;
}
