#include "FluidSolver.h"



FluidSolver::FluidSolver()
{
}


FluidSolver::FluidSolver(const DualMesh * mesh)
{
	this->mesh = mesh;
	const int nDualFaces = mesh->getNumberOfFaces();
	dt_local = Eigen::VectorXd(nDualFaces);
}

FluidSolver::~FluidSolver()
{
}

const double FluidSolver::updateFluidState()
{
	fluidStateUpdate = Eigen::MatrixXd::Zero(fluidStates.rows(), fluidStates.cols());
	const int nDualFaces = mesh->getNumberOfFaces();
	dt_local.setZero();

	// Compute fluid fluxes at faces
	for (int i = 0; i < nDualFaces; i++) {
		const Face * face = mesh->getFace(i);
		const Eigen::Vector3d faceNormal = face->getNormal();
		const double faceArea = face->getArea();

		const Face::FluidType faceFluidType =  face->getFluidType();

		switch (faceFluidType) {
		case Face::FluidType::INTERIOR:
			updateFaceFluxInterior(i);
			break;
			
		case Face::FluidType::OPENING:
			updateFaceFluxOpening(i);
			break;

		case Face::FluidType::WALL:
			updateFaceFluxWall(i);
			break;

		default:
			dt_local(i) = std::numeric_limits<double>::max();
		}
	}
	bool allPositive = (dt_local.array() > 0.).all();
	assert(allPositive);

	// get global timestep size
	const double dt = dt_local.minCoeff();
	assert(dt > 0);

	for (int i = 0; i < mesh->getNumberOfCells(); i++) {
		fluidStateUpdate.col(i) *= 1. / (mesh->getCell(i)->getVolume());
	}

	fluidStates -= dt * fluidStateUpdate;

	return dt;
}

void FluidSolver::updateFaceFluxInterior(const int faceIdx)
{
	const Face * face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::INTERIOR);

	const Eigen::Vector3d faceNormal = face->getNormal();
	const double faceArea = face->getArea();

	// Get left and right cell of this face
	std::vector<Cell*> faceCells = face->getCellList();
	assert(faceCells.size() == 2);
	const Cell * cell = faceCells[0];
	const int orientation = (face->getCenter() - cell->getCenter()).dot(faceNormal) > 0 ? 1 : -1;
	Cell * leftCell = nullptr;
	Cell * rightCell = nullptr;
	if (orientation > 0) {
		leftCell = faceCells[0];
		rightCell = faceCells[1];
	}
	else {
		leftCell = faceCells[1];
		rightCell = faceCells[0];
	}
	const int idxL = leftCell->getIndex();
	const int idxR = rightCell->getIndex();

	// distance between cell center and face center
	const Eigen::Vector3d fc = face->getCenter();
	const Eigen::Vector3d ccL = leftCell->getCenter();
	const Eigen::Vector3d ccR = rightCell->getCenter();
	const double dxL = std::abs((fc - ccL).dot(faceNormal));
	const double dxR = std::abs((fc - ccR).dot(faceNormal));
	const double dx = std::min(dxL, dxR); // minimum distance between cell center and face center
	assert(dx > 0);
	const Eigen::VectorXd qL = fluidStates.col(idxL);
	const Eigen::VectorXd qR = fluidStates.col(idxR);

	double dt_loc = 0;
	//const Eigen::VectorXd faceFlux = Numerics::fluidFlux_rusanov(qL, qR, faceNormal, dx, dt_loc);
	const Eigen::VectorXd faceFlux = getRusanovFlux(qL, qR, faceNormal, dx, dt_loc);
	assert(dt_loc > 0);
	dt_local(faceIdx) = dt_loc;

	fluidStateUpdate.col(idxL) += faceFlux * faceArea;
	fluidStateUpdate.col(idxR) -= faceFlux * faceArea;
}

/** 
* Update the hydrodynamic flux across a face with Opening boundary conditions, i.e., same state across the face.
*/
void FluidSolver::updateFaceFluxOpening(const int faceIdx)
{
	const Face * face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::OPENING);

	const Eigen::Vector3d faceNormal = face->getNormal();
	const double faceArea = face->getArea();

	const std::vector<Cell*> faceCells = face->getCellList();
	const int nFaceCells = faceCells.size();
	assert(nFaceCells == 1);
	Cell * cell = faceCells[0];
	assert(cell->getFluidType() == Cell::FluidType::FLUID);

	int idxC = cell->getIndex();
	const int orientation = (face->getCenter() - cell->getCenter()).dot(faceNormal) > 0 ? 1 : -1;

	Eigen::VectorXd qL, qR;
	if (orientation > 0) {
		qL = fluidStates.col(idxC);
		qR = qL;
	}
	else {
		qR = fluidStates.col(idxC);
		qL = qR;
	}
	const Eigen::Vector3d cc = cell->getCenter();
	const Eigen::Vector3d fc = face->getCenter();
	const double dx = std::abs((fc - cc).dot(faceNormal));
	double dt_loc = 0;
	const Eigen::VectorXd flux = getRusanovFlux(qL, qR, faceNormal, dx, dt_loc);
	assert(dt_loc > 0);
	dt_local(faceIdx) = dt_loc;
	fluidStateUpdate.col(idxC) += orientation * faceArea * flux;
}

/** 
* Update the hydrodynamic flux at a face with Wall boundary condition. Here, we consider only states that are perpendicular to wall. 
* Therefore, this is modeled with velocity in opposite flow direction.
*/
void FluidSolver::updateFaceFluxWall(const int faceIdx)
{
	const Face * face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::WALL);

	const Eigen::Vector3d faceNormal = face->getNormal();
	const double faceArea = face->getArea();

	const std::vector<Cell*> faceCells = face->getCellList();
	const int nFaceCells = faceCells.size();

	Cell * cell = nullptr; 
	switch (nFaceCells) {
	case 1:
		cell = faceCells[0];
		break;
	case 2:
		cell = (faceCells[0]->getFluidType() == Cell::FluidType::FLUID) ? faceCells[0] : faceCells[1];
		break;
	default:
		assert(nFaceCells >= 1 && nFaceCells <= 2);
	}
	assert(cell->getFluidType() == Cell::FluidType::FLUID);
	int idxC = cell->getIndex();

	const int orientation = (face->getCenter() - cell->getCenter()).dot(faceNormal) > 0 ? 1 : -1;

	Eigen::VectorXd qL, qR;
	if (orientation > 0) {
		qL = fluidStates.col(idxC);
		qR = qL;
		qR.segment(1, 3) *= -1;   /// Supposed to reverse the normal velocity?
	}
	else {
		qR = fluidStates.col(idxC);
		qL = qR;
		qL.segment(1, 3) *= -1;
	}
	const Eigen::Vector3d cc = cell->getCenter();
	const Eigen::Vector3d fc = face->getCenter();
	const double dx = std::abs((fc - cc).dot(faceNormal));
	double dt_loc = 0;
	const Eigen::VectorXd flux = getRusanovFlux(qL, qR, faceNormal, dx, dt_loc);
	assert(dt_loc > 0);
	dt_local(faceIdx) = dt_loc;
	fluidStateUpdate.col(idxC) += orientation * faceArea * flux;
}

