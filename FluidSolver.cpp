#include "FluidSolver.h"



FluidSolver::FluidSolver()
{
}

FluidSolver::FluidSolver(const DualMesh * mesh) : 
	mesh(mesh), 
 	nFaces(mesh->facet_counts.nF_interior + mesh->facet_counts.nF_opening + mesh->facet_counts.nF_wall),
	nCells(mesh->facet_counts.nC_fluid)
{
	U = Eigen::MatrixXd::Zero(nCells, 5);
	F = Eigen::MatrixXd::Zero(nFaces, 5);
	rate_of_change = Eigen::MatrixXd::Zero(nCells, 5);

	// init mapping : mesh index <---> vector component index
	U2cell_map = Eigen::VectorXi(nCells).setConstant(-1);
	cell2U_map = Eigen::VectorXi(mesh->getNumberOfCells()).setConstant(-1);
	F2face_map = Eigen::VectorXi(nFaces).setConstant(-1);
	face2F_map = Eigen::VectorXi(mesh->getNumberOfFaces()).setConstant(-1);

	int count = 0;
	for (Cell* c : mesh->getCells()) {
		if (c->getFluidType() == Cell::FluidType::Fluid) {
			U2cell_map[count] = c->getIndex();
			cell2U_map[c->getIndex()] = count;
			count++;
		}
	}
	assert(count == nCells);
	count = 0;
	for (Face* f : mesh->getFaces()) {
		if (f->getFluidType() != Face::FluidType::Undefined) {
			F2face_map[count] = f->getIndex();
			face2F_map[f->getIndex()] = count;
			count++;
		}
	}
	assert(count == nFaces);

}

FluidSolver::~FluidSolver()
{
}



const double FluidSolver::timeStepping()
{
	const int dt = updateFlux();
	updateRateOfChange();
	U += dt * rate_of_change;
	return dt;
}

void FluidSolver::applyInitialCondition() {
	Eigen::VectorXd qL(5), qR(5);
	qL << 1.0, 0.0, 0.0, 0.0, 1.0;
	qR << 0.125, 0.0, 0.0, 0.0, 0.1;
	qL = primitive2conservative(qL);
	qR = primitive2conservative(qR);
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		if (mesh->getCell(U2cell(U_idx))->getCenter()[2] < 0) {
			U.row(U_idx) = qL;
		}
		else {
			U.row(U_idx) = qR;
		}
	} 
}

const double FluidSolver::updateFlux() {
	double dt, min_dt = 1e10;
	for (int F_idx = 0; F_idx < nFaces; F_idx++) {
		int faceIdx = F2face(F_idx);
		switch(mesh->getFace(faceIdx)->getFluidType()) {
			case Face::FluidType::Interior : dt = updateFluxInterior(faceIdx); break;
			case Face::FluidType::Opening  : dt = updateFluxOpening (faceIdx); break;
			case Face::FluidType::Wall     : dt = updateFluxWall    (faceIdx); break;  
		}  
		min_dt = dt < min_dt ? dt : min_dt;
	}
	return min_dt;
}

void FluidSolver::updateRateOfChange() {
	rate_of_change.setZero();
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		const int cell_idx = U2cell(U_idx);
		const Cell* cell = mesh->getCell(cell_idx);
		assert(cell->getFluidType() == Cell::FluidType::Fluid);
		for (const Face* face : cell->getFaceList()) {
			assert(face->getFluidType() != Face::FluidType::Undefined);
			const int face_idx = face->getIndex();
			const int F_idx = face2F(face_idx);
			const int incidence = cell->getOrientation(face);
			rate_of_change.row(U_idx) -= incidence * face->getArea() * F.row(F_idx);
		}
		rate_of_change.row(U_idx) /= cell->getVolume();
	}
}

const double FluidSolver::updateFluxInterior(const int faceIdx)
{
	const Face* face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::Interior);

	const Eigen::Vector3d faceNormal = face->getNormal();

	// Get left and right cell of this face
	std::vector<Cell*> faceCells = face->getCellList();
	Cell *leftCell, *rightCell;
	assert(faceCells.size() == 2);
	if (faceCells[0]->getOrientation(face) > 0) {
		leftCell  = faceCells[0];
		rightCell = faceCells[1]; 
	}
	else {
		leftCell  = faceCells[1];
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
	const Eigen::VectorXd qL = U.row(cell2U(idxL));
	const Eigen::VectorXd qR = U.row(cell2U(idxR));

	// update flux
	F.row(face2F(faceIdx)) = RusanovFlux(qL, qR, faceNormal);
	const double s = std::max(maxWaveSpeed(qL, faceNormal), maxWaveSpeed(qR, faceNormal));
	
	return dx / s;
}

const double FluidSolver::updateFluxOpening(const int faceIdx)
{
	const Face* face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::Opening);

	const Eigen::Vector3d faceNormal = face->getNormal();
	assert((faceNormal.cross(Eigen::Vector3d::UnitZ())).norm() == 0);  // assert opening face is alway parallel to z-axis

	const Cell* faceCell = face->getCellList()[0];
	Eigen::VectorXd q = U.row(cell2U(faceCell->getIndex()));

	// distance between cell center and face center
	const Eigen::Vector3d fc = face->getCenter();
	const Eigen::Vector3d cc = faceCell->getCenter();
	const double dx = (fc - cc).norm();
	assert(dx > 0);

	F.col(face2F(faceIdx)) = RusanovFlux(q, q, faceNormal);
	const double s = maxWaveSpeed(q, faceNormal);
	
	return dx / s;
}

const double FluidSolver::updateFluxWall(const int faceIdx)
{
	const Face* face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::Wall);

	const Eigen::Vector3d faceNormal = face->getNormal();
	assert(std::abs(faceNormal[2]) < std::numeric_limits<double>::epsilon());  // assert wall face is alway perpendicular to z-axis

	const Cell* faceCell = face->getCellList()[0];
	Eigen::VectorXd qL, qR;
	if (faceCell->getOrientation(face) > 0) {
		qL = U.row(cell2U(faceCell->getIndex()));
		qR = qL;
		qR.segment(1,3) -= 2 * (qR.segment(1,3).dot(faceNormal)) * faceNormal;
	}
	else {
		qR = U.row(cell2U(faceCell->getIndex()));
		qL = qR;
		qL.segment(1,3) -= 2 * (qL.segment(1,3).dot(faceNormal)) * faceNormal;
	}

	// distance between cell center and face center
	const Eigen::Vector3d fc = face->getCenter();
	const Eigen::Vector3d cc = faceCell->getCenter();
	const double dx = (fc - cc).norm();
	assert(dx > 0);

	F.col(face2F(faceIdx)) = RusanovFlux(qL, qR, faceNormal);
	const double s = maxWaveSpeed(qL, faceNormal); // Both sides have the same wave speed

	return dx / s;
}

Eigen::VectorXd FluidSolver::Flux(const Eigen::VectorXd &q, const Eigen::Vector3d &fn) const {
	Eigen::VectorXd flux(5);
	Eigen::VectorXd q_prim = conservative2primitive(q);
	flux[0] = q.segment(1,3).dot(fn);
	flux.segment(1,3) = flux[0]*q_prim.segment(1,3) + q_prim[4] * fn / vareps2;
	flux[4] = (q[4] + q_prim[4] / vareps2) * q_prim.segment(1,3).dot(fn);
	return flux;
}

Eigen::VectorXd FluidSolver::RusanovFlux(const Eigen::VectorXd & qL, 
							             const Eigen::VectorXd & qR, 
										 const Eigen::Vector3d & fn) const {
	double s = std::max(maxWaveSpeed(qL, fn), maxWaveSpeed(qR, fn));
	return 0.5 * (Flux(qL, fn) + Flux(qR, fn)) - 0.5 * s * (qR - qL);
											
}

Eigen::VectorXd FluidSolver::conservative2primitive(const Eigen::VectorXd &q_cons) const {
	Eigen::VectorXd q_prim(5);
	q_prim[0] = q_cons[0];
	q_prim.segment(1,3) = q_cons.segment(1,3) / q_cons[0];
	q_prim[4] = (gamma - 1) * (q_cons[4] - 0.5*q_cons[0]*q_prim.segment(1,3).norm()) * vareps2;
	assert(q_prim[4] > 0);
	return q_prim;
}

Eigen::VectorXd FluidSolver::primitive2conservative(const Eigen::VectorXd &q_prim) const {
	Eigen::VectorXd q_cons(5);
	q_cons[0] = q_prim[0];
	q_cons.segment(1,3) = q_prim.segment(1,3) * q_prim[0];
	q_cons[4] = 0.5*q_prim[0]*q_prim.segment(1,3).norm() + q_prim[4] / (gamma - 1) / vareps2;
	return q_cons;
}

const double FluidSolver::speedOfSound(const Eigen::VectorXd &q_cons) const {
	Eigen::VectorXd q_prim = conservative2primitive(q_cons);
	return std::sqrt(gamma * q_prim[4] / q_cons[0] / vareps2);
}

const double FluidSolver::maxWaveSpeed(const Eigen::VectorXd &q_cons, const Eigen::Vector3d &normal) const {
	const double sos = speedOfSound(q_cons);
	const double velocity = (q_cons.segment(1,3).dot(normal)) / q_cons[0];
	return sos + velocity;
};

void FluidSolver::writeSnapshot(H5Writer & writer) const {
	Eigen::MatrixXd attributeOfCell(mesh->getNumberOfCells(), 5);
	attributeOfCell.setZero();  
	for (int cell_idx = 0; cell_idx < mesh->getNumberOfCells(); cell_idx++) {
		if (mesh->getCell(cell_idx)->getFluidType() == Cell::FluidType::Fluid) {
			attributeOfCell.row(cell_idx) = conservative2primitive(U.row(cell2U(cell_idx)));
		}
	}
	writer.writeDoubleVector(attributeOfCell.col(0), "/density");
	writer.writeDoubleVector(attributeOfCell.col(4), "/pressure");
	writer.writeDoubleMatrix(attributeOfCell.block(0,1,mesh->getNumberOfCells(),3), "/velocity");
}
