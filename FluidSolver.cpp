#include "FluidSolver.h"

FluidSolver::FluidSolver(const DualMesh * mesh) : 
	mesh(mesh), 
 	nFaces(mesh->facet_counts.nF_interior + mesh->facet_counts.nF_opening + mesh->facet_counts.nF_wall + mesh->facet_counts.nF_mixed),
	nCells(mesh->facet_counts.nC_fluid)
{
	U 	= Eigen::MatrixXd::Zero(nCells, 5);
	F 	= Eigen::MatrixXd::Zero(nFaces, 5);
	rhs = Eigen::MatrixXd::Zero(nCells, 5);
	rate_of_change = Eigen::MatrixXd::Zero(nCells, 5);
	S   = Eigen::VectorXd::Zero(nFaces);

	eta = Eigen::MatrixXd::Zero(mesh->getNumberOfCells(), 3);

	// init mapping : mesh index <---> vector component index
	U2cell_map = Eigen::VectorXi(nCells).setConstant(-1);
	cell2U_map = Eigen::VectorXi(mesh->getNumberOfCells()).setConstant(-1);
	F2face_map = Eigen::VectorXi(nFaces).setConstant(-1);
	face2F_map = Eigen::VectorXi(mesh->getNumberOfFaces()).setConstant(-1);

	int count = 0;
	for (const Cell* c : mesh->getCells()) {
		if (c->getFluidType() == Cell::FluidType::Fluid) {
			U2cell_map[count] = c->getIndex();
			cell2U_map[c->getIndex()] = count;
			count++;
		}
	}
	assert(count == nCells);
	count = 0;
	for (const Face* f : mesh->getFaces()) {
		if (f->getFluidType() != Face::FluidType::Undefined) {
			F2face_map[count] = f->getIndex();
			face2F_map[f->getIndex()] = count;
			count++;
		}
	}
	assert(count == nFaces);
}

// Since member variables can not be initialized in constructor delegated to another one, I have to copy and pasta the body.
// TODO: A more elegent way?
FluidSolver::FluidSolver(const DualMesh* mesh, const double gamma, const double mass, const double charge, const std::string name) : 
	gamma(gamma), vareps2(mass), charge(charge), name(name),
  	mesh(mesh), 
  	nFaces(mesh->facet_counts.nF_interior + mesh->facet_counts.nF_opening + mesh->facet_counts.nF_wall + mesh->facet_counts.nF_mixed),
  	nCells(mesh->facet_counts.nC_fluid) 
{
	U 	= Eigen::MatrixXd::Zero(nCells, 5);
	F 	= Eigen::MatrixXd::Zero(nFaces, 5);
	rhs = Eigen::MatrixXd::Zero(nCells, 5);
	rate_of_change = Eigen::MatrixXd::Zero(nCells, 5);
	S   = Eigen::VectorXd::Zero(nFaces);

	eta = Eigen::MatrixXd::Zero(mesh->getNumberOfCells(), 3);

	// init mapping : mesh index <---> vector component index
	U2cell_map = Eigen::VectorXi(nCells).setConstant(-1);
	cell2U_map = Eigen::VectorXi(mesh->getNumberOfCells()).setConstant(-1);
	F2face_map = Eigen::VectorXi(nFaces).setConstant(-1);
	face2F_map = Eigen::VectorXi(mesh->getNumberOfFaces()).setConstant(-1);

	int count = 0;
	for (const Cell* c : mesh->getCells()) {
		if (c->getFluidType() == Cell::FluidType::Fluid) {
			U2cell_map[count] = c->getIndex();
			cell2U_map[c->getIndex()] = count;
			count++;
		}
	}
	assert(count == nCells);
	count = 0;
	for (const Face* f : mesh->getFaces()) {
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

void FluidSolver::timeStepping(const double dt) {
	updateRateOfChange(false);
	U += dt * rate_of_change;
	assert(isValidState());
	if (!isValidState()) {
		std::cout << "******************************" << std::endl;
		std::cout << "*   Fluid State not valid!   *" << std::endl;
		std::cout << "******************************" << std::endl; 
	}
}

void FluidSolver::timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B) {
	rhs.setZero();
	applyLorentzForce(E, B);
	updateRateOfChange(true);
	U += dt * rate_of_change;

	// check_updatedMomentum();
	assert(isValidState());
	if (!isValidState()) {
		std::cout << "******************************" << std::endl;
		std::cout << "*   Fluid State not valid!   *" << std::endl;
		std::cout << "******************************" << std::endl; 
	}
}

void FluidSolver::timeStepping(const double dt, 
							   const Eigen::MatrixXd &E, 
							   const Eigen::MatrixXd &B,
							   const double alpha,
							   const FluidSolver* anotherSpecies) {
	rhs.setZero();
	applyLorentzForce(E, B);
	applyFrictionTerm(anotherSpecies, alpha);
	updateRateOfChange(true);
	U += dt * rate_of_change;

	// check_updatedMomentum();
	assert(isValidState());
	if (!isValidState()) {
		std::cout << "******************************" << std::endl;
		std::cout << "*   Fluid State not valid!   *" << std::endl;
		std::cout << "******************************" << std::endl; 
	}
}

void FluidSolver::applyInitialCondition() {
	Eigen::VectorXd qL(5), qR(5);
	if (name.compare("electron") == 0) {
		qL << 1.0, 0.0, 0.0, 0.0, 1.0;
		qR << 1.0, 0.0, 0.0, 0.0, 1.0;
	}
	else if (name.compare("ion") == 0) {
		qL << 1.0, 0.0, 0.0, 0.0, 1.0;
		qR << 1.0, 0.0, 0.0, 0.0, 1.0;
	}
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

void FluidSolver::applyInitialCondition(const std::string h5_file) {
	H5Reader reader(h5_file);
	Eigen::VectorXd n = reader.readVectorData("/" + name + "-density");
	Eigen::MatrixXd vel = reader.readMatrixData("/" + name + "-velocity");
	Eigen::VectorXd p = reader.readVectorData("/" + name + "-pressure");
	Eigen::MatrixXd temp(mesh->getNumberOfCells(), 5);
	temp.col(0) = n;
	temp.middleCols(1,3) = vel;
	temp.col(4) = p;
	for (int i = 0; i < nCells; i++) {
		U.row(i) = primitive2conservative(temp.row(U2cell(i)));
	} 
}

double FluidSolver::updateFluxExplicit() {
	double dt, min_dt = 1e10;
	for (int F_idx = 0; F_idx < nFaces; F_idx++) {
		const int faceIdx = F2face(F_idx);
		switch(mesh->getFace(faceIdx)->getFluidType()) {
			case Face::FluidType::Interior : dt = updateFluxInterior(faceIdx); break;
			case Face::FluidType::Opening  : dt = updateFluxOpening (faceIdx); break;
			case Face::FluidType::Wall     : dt = updateFluxWall    (faceIdx); break;
			case Face::FluidType::Mixed    : dt = updateFluxMixed   (faceIdx); break;
		}  
		min_dt = dt < min_dt ? dt : min_dt;
	}
	assert(min_dt > 0);
	return min_dt;
}

void FluidSolver::updateMassFluxImplicit() {
	const Eigen::VectorXd n_extended = getExtended_n();
	for (int F_idx = 0; F_idx < nFaces; F_idx++) {
		const int faceIdx = F2face(F_idx);
		const std::vector<Cell*> cells = mesh->getFace(faceIdx)->getCellList();
		const Eigen::Vector3d normal = mesh->getFace(faceIdx)->getNormal();
		switch(mesh->getFace(faceIdx)->getFluidType()) {
			case Face::FluidType::Interior : {
				assert(cells.size() == 2);
				const int cell1_idx = cells[0]->getIndex();
				const int cell2_idx = cells[1]->getIndex();
				Eigen::Vector3d nu1_old = U.row(cell2U(cell1_idx)).segment(1,3);
				Eigen::Vector3d nu2_old = U.row(cell2U(cell2_idx)).segment(1,3);
				Eigen::Vector3d nu1_new = updatedMomentum.row(cell1_idx);
				Eigen::Vector3d nu2_new = updatedMomentum.row(cell2_idx);
				F(F_idx, 0) += 0.5 * (nu1_new + nu2_new - nu1_old - nu2_old).dot(normal); 
				break;
			}
			case Face::FluidType::Opening : {
				assert(cells.size() == 1);
				const int cell_idx = cells[0]->getIndex();
				Eigen::Vector3d nu_old = U.row(cell2U(cell_idx)).segment(1,3);
				Eigen::Vector3d nu_new = updatedMomentum.row(cell_idx);
				F(F_idx, 0) += (nu_new - nu_old).dot(normal);
				break;
			}
			case Face::FluidType::Wall : {
				assert(cells.size() == 1);
				// do nothing, since the advection part is always zero
				break;
			}
			case Face::FluidType::Mixed : {
				assert(cells.size() == 1);
				for (const Face* subf : mesh->getFace(faceIdx)->getSubFaceList()) {
					if (subf->getFluidType() == Face::FluidType::Wall) {
						// do nothing, since the advection part is always zero
					}
					else if (subf->getFluidType() == Face::FluidType::Opening) {
						const int cell_idx = cells[0]->getIndex();
						Eigen::Vector3d nu_old = U.row(cell2U(cell_idx)).segment(1,3);
						Eigen::Vector3d nu_new = updatedMomentum.row(cell_idx);
						F(F_idx, 0) += (nu_new - nu_old).dot(normal) * subf->getArea() / mesh->getFace(faceIdx)->getArea();
					}
					else {
						assert(false && "sub face should not be of type other than <Opening> and <Wall>");
					}
				}
				break;
			}
			default : assert(false); break;
		}
	}	
}

void FluidSolver::updateMomentum(const double dt, const Eigen::MatrixXd &E) {
	updatedMomentum = dt / vareps2 * charge * getExtended_n().asDiagonal() * E + eta;
}


void FluidSolver::updateMomentum(const double dt, const Eigen::MatrixXd &E, const double alpha, const FluidSolver* anotherSpecies) {
	Eigen::VectorXd n_this  = getExtended_n();
	Eigen::VectorXd n_other = anotherSpecies->getExtended_n(); 
	Eigen::VectorXd temp1 = 
			(1 + dt * alpha / anotherSpecies->vareps2 * n_this.array() 
		    + dt * alpha / anotherSpecies->vareps2 * anotherSpecies->charge / charge * n_other.array())  /
			(1 + dt * alpha * (1.0 / anotherSpecies->vareps2 * n_this.array() + 1.0 / vareps2 * n_other.array())) * 
			(dt / vareps2 * charge * n_this.array());
	Eigen::VectorXd temp2 = (dt * alpha / vareps2 * n_this).array() 
				/ (1 + dt * alpha * (n_this / anotherSpecies->vareps2 + n_other / vareps2).array());
	Eigen::VectorXd temp3 = (1 + dt * alpha / anotherSpecies->vareps2 * n_this.array()) 
	            / (1 + dt * alpha * (n_this / anotherSpecies->vareps2 + n_other / vareps2).array());
	updatedMomentum = temp1.asDiagonal() * E 
					+ temp2.asDiagonal() * anotherSpecies->eta 
					+ temp3.asDiagonal() * eta;
}


void FluidSolver::applyLorentzForce(const Eigen::MatrixXd &E, const Eigen::MatrixXd &B) {
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		Eigen::Vector3d m_vec = U.row(U_idx).segment(1,3);  // momentum vector
		Eigen::Vector3d B_vec = B.row(U2cell(U_idx));       // B-field vector
		Eigen::Vector3d E_vec = E.row(U2cell(U_idx));       // E-field vector
		rhs.row(U_idx).segment(1,3) += 1.0 / vareps2 * charge * (U(U_idx, 0) * E_vec + m_vec.cross(B_vec));
		rhs(U_idx, 4) += 1.0 / vareps2 * charge * m_vec.dot(E_vec);
	}
}

void FluidSolver::applyFrictionTerm(const FluidSolver* anotherSpecies, const double alpha) {
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		const Eigen::Vector3d temp = 1.0 / vareps2 * alpha * 
									(U(U_idx, 0) * anotherSpecies->updatedMomentum.row(U2cell(U_idx))
								    - anotherSpecies->U(U_idx, 0) * updatedMomentum.row(U2cell(U_idx)));
		rhs.row(U_idx).segment(1,3) += temp;
		rhs(U_idx, 4) += temp.dot(U.row(U_idx).segment(1,3));
	}
}

void FluidSolver::updateRateOfChange(const bool with_rhs) {
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
			if (!face->isPlane()) {
				for (const Face* subf : face->getSubFaceList()) {
					assert(cell->getOrientation(subf) == incidence);
				}
			}
			rate_of_change.row(U_idx) -= incidence * face->getArea() * F.row(F_idx);
		}
		rate_of_change.row(U_idx) /= cell->getVolume();
	}
	if (with_rhs) {
		rate_of_change += rhs;
	}
}

double FluidSolver::updateFluxInterior(const int faceIdx)
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
	const Eigen::VectorXd qL = U.row(cell2U(idxL));
	const Eigen::VectorXd qR = U.row(cell2U(idxR));

	// update flux
	F.row(face2F(faceIdx)) = RusanovFlux(qL, qR, faceNormal);
	const double s = std::max(maxWaveSpeed(qL, faceNormal), maxWaveSpeed(qR, faceNormal));
	S[face2F(faceIdx)] = s;
	
	assert(dx > 0);
	assert(s > 0);
	
	return dx / s;
}

double FluidSolver::updateFluxOpening(const int faceIdx)
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

	F.row(face2F(faceIdx)) = RusanovFlux(q, q, faceNormal);
	const double s = maxWaveSpeed(q, faceNormal);
	S[face2F(faceIdx)] = s;

	assert(dx > 0);
	assert(s > 0);
	
	return dx / s;
}

double FluidSolver::updateFluxWall(const int faceIdx)
{
	const Face* face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::Wall);

	const Eigen::Vector3d faceNormal = face->getNormal();
	//assert(std::abs(faceNormal[2]) < 100 * std::numeric_limits<double>::epsilon());  // assert wall face is alway perpendicular to z-axis

	const Cell* faceCell;
	if (face->getCellList()[0]->getFluidType() == Cell::FluidType::Fluid) {
		faceCell = face->getCellList()[0];
	}
	else {
		faceCell = face->getCellList()[1];
	}
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

	F.row(face2F(faceIdx)) = RusanovFlux(qL, qR, faceNormal);
	const double s = maxWaveSpeed(qL, faceNormal); // Both sides have the same wave speed
	S[face2F(faceIdx)] = s;
	assert(dx > 0);
	assert(s > 0);

	return dx / s;
}

double FluidSolver::updateFluxMixed(const int faceIdx) {
	const Face* face = mesh->getFace(faceIdx);
	assert(face->getFluidType() == Face::FluidType::Mixed);
	assert(face->getSubFaceList().size() >= 2 && face->getSubFaceList().size() <= 3);
	assert(face->isBoundary()); // In our setup, mixed (non-plane) face is at the boundary.
	assert(face->getCellList().size() == 1);
	int nWalls = 0, nOpen = 0;
	for (const Face* subf : face->getSubFaceList()) {
		if (subf->getFluidType() == Face::FluidType::Opening) {
			nOpen++;
		}
		else if (subf->getFluidType() == Face::FluidType::Wall) {
			nWalls++;
		}
		else {
			assert(false);
		}
	}

	const Cell* faceCell = face->getCellList()[0];
	const Eigen::Vector3d cc = faceCell->getCenter();
	std::vector<double> dt_;  // A list to hold the allowed dt for each sub-face

	Eigen::VectorXd qL, qR;
	F.row(face2F(faceIdx)).setZero();  // reset the flux
	S[face2F(faceIdx)] = 0;            // reset artifical dissipation
	for (const Face* subf : face->getSubFaceList()) {
		const Eigen::VectorXd subFaceNormal = subf->getNormal();
		assert(subFaceNormal.dot(face->getNormal()) > 0);
		const Eigen::Vector3d subfc = subf->getCenter();
		if (subf->getFluidType() == Face::FluidType::Wall) {
			if (faceCell->getOrientation(subf) > 0) {
				qL = U.row(cell2U(faceCell->getIndex()));
				qR = qL;
				qR.segment(1,3) -= 2 * (qR.segment(1,3).dot(subFaceNormal)) * subFaceNormal;
			}
			else {
				qR = U.row(cell2U(faceCell->getIndex()));
				qL = qR;
				qL.segment(1,3) -= 2 * (qL.segment(1,3).dot(subFaceNormal)) * subFaceNormal;
			}
		}
		else if (subf->getFluidType() == Face::FluidType::Opening) {
			qL = U.row(cell2U(faceCell->getIndex()));
			qR = qL;
		}
		else {
			assert(false);
		}
		const double dx = (subfc - cc).norm();
		const double s = maxWaveSpeed(qL, subFaceNormal); // Both sides have the same wave speed
		F.row(face2F(faceIdx)) += RusanovFlux(qL, qR, subFaceNormal) * subf->getArea();
		S[face2F(faceIdx)] += s * subf->getArea();
		dt_.push_back(dx / s);
	}
	F.row(face2F(faceIdx)) /= face->getArea();
	S[face2F(faceIdx)] /= face->getArea();

	return *std::min_element(dt_.begin(), dt_.end()); // return the least allowed dt
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
	q_prim[4] = (gamma - 1) * (q_cons[4] - 0.5*q_cons[0]*q_prim.segment(1,3).squaredNorm()) * vareps2;
	assert(q_prim[4] > 0);
	return q_prim;
}

Eigen::VectorXd FluidSolver::primitive2conservative(const Eigen::VectorXd &q_prim) const {
	Eigen::VectorXd q_cons(5);
	q_cons[0] = q_prim[0];
	q_cons.segment(1,3) = q_prim.segment(1,3) * q_prim[0];
	q_cons[4] = 0.5*q_prim[0]*q_prim.segment(1,3).squaredNorm() + q_prim[4] / (gamma - 1) / vareps2;
	return q_cons;
}

double FluidSolver::speedOfSound(const Eigen::VectorXd &q_cons) const {
	Eigen::VectorXd q_prim = conservative2primitive(q_cons);
	return std::sqrt(gamma * q_prim[4] / q_cons[0] / vareps2);
}

double FluidSolver::maxWaveSpeed(const Eigen::VectorXd &q_cons, const Eigen::Vector3d &normal) const {
	const double sos = speedOfSound(q_cons);
	const double velocity = std::abs((q_cons.segment(1,3).dot(normal))) / q_cons[0];
	return sos + velocity;
};

bool FluidSolver::isValidState() const {
	Eigen::VectorXd q_prim(5);
	for (int i = 0; i < nCells; i++) {
		q_prim = conservative2primitive(U.row(i));
		if (q_prim[0] < 0 || q_prim[4] < 0) return false;
	}
	return true;
}

void FluidSolver::update_eta(const double dt, const Eigen::MatrixXd &B) const {
	Eigen::Vector3d m_vec, B_vec;
	Eigen::MatrixXd temp = U.block(0, 1, nCells, 3);
	temp += dt * rate_of_change.block(0, 1, nCells, 3);	 // At this moment, <rate_of_change> is merely due to fluxes without rhs.
	// TODO: This row-wise cross production needs a more efficient implementation
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		m_vec = U.row(U_idx).segment(1,3);
		B_vec = B.row(U2cell(U_idx));
		temp.row(U_idx) += m_vec.cross(B_vec) * dt * charge / vareps2;
	}
	// We need to extend the row size to number of cells with non-fluid entry being zero vector
	eta.setZero();
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		eta.row(U2cell(U_idx)) = temp.row(U_idx);
	} 
}

Eigen::VectorXd FluidSolver::get_mu(const double dt, 
									const Eigen::MatrixXd &B, 
									const Tensor3& A, 
									const Eigen::SparseMatrix<double>& D) const {
	update_eta(dt, B);
	// First we need extended vector of number density.
	Eigen::VectorXd n_extended = getExtended_n();
	Eigen::VectorXd s_extended = getExtended_s();
	return 0.5 * A.twoContract(eta) - 0.5 * (D * n_extended).cwiseProduct(s_extended);
}

Eigen::VectorXd FluidSolver::get_mu(const double dt, 
		                            const Eigen::MatrixXd &B, 
							   		const Tensor3& A, 
							   		const Eigen::SparseMatrix<double>& D,
							   		const double alpha,
							  		const FluidSolver* anotherSpecies) const {
	assert(anotherSpecies != this);
	update_eta(dt, B);
	anotherSpecies->update_eta(dt, B);  // For easy implementation this function is called twice, which is not necessary.
	Eigen::VectorXd n_this  = getExtended_n();
	Eigen::VectorXd n_other = anotherSpecies->getExtended_n(); 
	Eigen::VectorXd temp1 = (dt * alpha / vareps2 * n_this).array() 
				/ (1 + dt * alpha * (n_this / anotherSpecies->vareps2 + n_other / vareps2).array());
	Eigen::VectorXd temp2 = (1 + dt * alpha / anotherSpecies->vareps2 * n_this.array()) 
	            / (1 + dt * alpha * (n_this / anotherSpecies->vareps2 + n_other / vareps2).array());
	return 0.5 * A.twoContract(temp1.asDiagonal() * anotherSpecies->eta + temp2.asDiagonal() * eta) 
	       - 0.5 * (D * n_this).cwiseProduct(getExtended_s());
}

Eigen::SparseMatrix<double> FluidSolver::get_T(const double dt,
											   const Tensor3& A,
											   const Tensor3& R) const {
	// First we need extended vector of number density.
	Eigen::VectorXd n_extended = getExtended_n();
	n_extended *= dt * charge / vareps2;
	return 0.5 * A.twoContract(R.firstDimWiseProduct(n_extended));
}

Eigen::SparseMatrix<double> FluidSolver::get_T(const double dt,
											   const Tensor3& A,
											   const Tensor3& R,
											   const double alpha,
											   const FluidSolver* anotherSpecies) const {
	assert(anotherSpecies != this);
	const Eigen::VectorXd n_this  = getExtended_n();
	const Eigen::VectorXd n_other = anotherSpecies->getExtended_n(); 
	Eigen::VectorXd temp = 
			(1 + dt * alpha / anotherSpecies->vareps2 * n_this.array() 
		    + dt * alpha / anotherSpecies->vareps2 * anotherSpecies->charge / charge * n_other.array())  /
			(1 + dt * alpha * (1.0 / anotherSpecies->vareps2 * n_this.array() + 1.0 / vareps2 * n_other.array())) * 
			(dt / vareps2 * charge * n_this.array());
	return 0.5 * A.twoContract(R.firstDimWiseProduct(temp));
}

Eigen::VectorXd FluidSolver::getExtended_n() const {
	Eigen::VectorXd n_extended(mesh->getNumberOfCells());
	n_extended.setZero();
	for (int U_idx = 0; U_idx < nCells; U_idx++) {
		n_extended[U2cell(U_idx)] = U(U_idx, 0);
	}
	return n_extended;
}

Eigen::VectorXd FluidSolver::getExtended_s() const {
	Eigen::VectorXd S_extended(mesh->getNumberOfFaces());
	S_extended.setZero();
	for (int F_idx = 0; F_idx < nFaces; F_idx++) {
		S_extended[F2face(F_idx)] = S[F_idx];
	}
	return S_extended;
}

void FluidSolver::writeSnapshot(H5Writer & writer) const {
	Eigen::MatrixXd attributeOfCell(mesh->getNumberOfCells(), 5);
	attributeOfCell.setZero();  
	for (int cell_idx = 0; cell_idx < mesh->getNumberOfCells(); cell_idx++) {
		if (mesh->getCell(cell_idx)->getFluidType() == Cell::FluidType::Fluid) {
			attributeOfCell.row(cell_idx) = conservative2primitive(U.row(cell2U(cell_idx)));
		}
	}
	writer.writeDoubleVector(attributeOfCell.col(0), "/" + name + "-density");
	writer.writeDoubleVector(attributeOfCell.col(4), "/" + name + "-pressure");
	writer.writeDoubleMatrix(attributeOfCell.block(0,1,mesh->getNumberOfCells(),3), "/" + name + "-velocity");
}

Eigen::VectorXd FluidSolver::getNorms() const {
	double u_1norm = 0, u_2norm = 0;
	double n_1norm = 0, n_2norm = 0;
	double p_1norm = 0, p_2norm = 0;
	for (int i = 0; i < nCells; i++) {
		const double vol = mesh->getCell(U2cell(i))->getVolume();
		const Eigen::VectorXd q_prim = conservative2primitive(U.row(i));
		n_1norm += std::abs(q_prim[0]) * vol;
		n_2norm += q_prim[0] * q_prim[0] * vol;
		u_1norm += q_prim.segment(1,3).norm() * vol;
		u_2norm += q_prim.segment(1,3).squaredNorm() * vol;
		p_1norm += std::abs(q_prim[4]) * vol;
		p_2norm += q_prim[4] * q_prim[4] * vol;
	}
	n_2norm = std::sqrt(n_2norm);
	u_2norm = std::sqrt(u_2norm);
	p_2norm = std::sqrt(p_2norm);
	Eigen::VectorXd norms(6);
	norms << n_1norm, n_2norm, u_1norm, u_2norm, p_1norm,  p_2norm;
	return norms;
}

void FluidSolver::check_A_and_D(const Tensor3& A, const Eigen::MatrixXd& D) const {
	Eigen::MatrixXd nu_extended(mesh->getNumberOfCells(), 3);
	nu_extended.setZero();
	for (int i = 0; i < nCells; i++) {
		nu_extended.row(U2cell(i)) = U.row(i).segment(1,3);
	}
	Eigen::VectorXd n_extended = getExtended_n();
	Eigen::VectorXd s_extended = getExtended_s();
	Eigen::VectorXd mass_flux_extended = 0.5 * A.twoContract(nu_extended) - 0.5 * s_extended.cwiseProduct(D * n_extended);
	bool flag = true;
	for (int i = 0; i < nFaces; i++) {
		if (std::abs(mass_flux_extended[F2face(i)] - F(i, 0)) > 1e-10) {
			flag = false;
		}
	}
	if (flag) {
		std::cout << ">>>>>>>>>>>>>>>>>>>>" << name << " A, D checked <<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout << ">>>>>>>>>>>>>>>>>>>>" << name << " A, D checking failed <<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
}

void FluidSolver::check_eta() const {
	update_eta(1.0, Eigen::MatrixXd::Zero(mesh->getNumberOfCells(), 3));
	bool flag = true;
	for (int i = 0; i < nCells; i++) {
		if ((eta.row(U2cell(i)) - U.row(i).segment(1,3) - rate_of_change.row(i).segment(1,3)).norm() > 1e-10) {
			flag = false;
		}
	}
	if (flag) {
		std::cout << ">>>>>>>>>>>>>>>>>>>>" << name << " eta checked <<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout << ">>>>>>>>>>>>>>>>>>>> " << name << " eta checking failed <<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
}

void FluidSolver::check_eta2(const Eigen::MatrixXd &B, const double dt) {
	rhs.setZero();
	applyLorentzForce(Eigen::MatrixXd::Zero(nCells, 3), B);
	updateRateOfChange(true);
	Eigen::MatrixXd temp = U + dt * rate_of_change;
	double error = 0;
	for (int i = 0; i < nCells; i++) {
		error += (eta.row(U2cell(i)) - temp.row(i).segment(1,3)).norm();
	}
	if (error < 1e-10) {
		std::cout << ">>>>>>>>>>>>>>>>>>>>" << name << " eta2 checked <<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout << ">>>>>>>>>>>>>>>>>>>> " << name << " eta2 checking failed <<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	}
}

void FluidSolver::check_updatedMomentum() const {
	bool flag = true;
	for (int i = 0; i < nCells; i++) {
		if ((U.row(i).segment(1,3) - updatedMomentum.row(U2cell(i))).norm() > 1e-12) {
			flag = false;
			//std::cout << "difference = " << (U.row(i).segment(1,3) - updatedMomentum.row(U2cell(i))).norm() << std::endl;
			//std::cout << "cell location = " << mesh->getCell(U2cell(i))->getCenter() << std::endl;
		}
	}
	if (flag) {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum checked <<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum checking failed <<<<<<<<<<<<<<<<<<" << std::endl;
	}
}

void FluidSolver::check_updatedMomentum2(const double dt, const Eigen::MatrixXd &E, const FluidSolver* another) const {
	double error;
	error = ((1 + dt * another->getExtended_n().array()).matrix().asDiagonal() * updatedMomentum 
			- dt * getExtended_n().asDiagonal() * another->updatedMomentum
			- dt * charge * getExtended_n().asDiagonal() * E - eta).norm();
	if (error < 1e-10) {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum2 checked <<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum2 checking failed <<<<<<<<<<<<<<<<<<" << std::endl;
	}
}

void FluidSolver::check_updatedMomentum3(const double dt, const Eigen::MatrixXd &E, const FluidSolver* another) const {
	double error;
	error = (updatedMomentum - dt * charge * getExtended_n().asDiagonal() * E 
		  - dt * (getExtended_n().asDiagonal() * another->updatedMomentum - another->getExtended_n().asDiagonal() * updatedMomentum)
		  - eta).norm();
	if (error < 1e-10) {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum3 checked <<<<<<<<<<<<<<<<<<" << std::endl;
	}
	else {
		std::cout  << ">>>>>>>>>>>>>>>>>>>>" << name << " updatedMomentum3 checking failed <<<<<<<<<<<<<<<<<<" << std::endl;
	} 
}