#include "TwoFluidSolver.h"

TwoFluidSolver::TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh, const Interpolator* interpolator) :
    primal         (primalMesh),
    dual           (dualMesh),
    electron_solver(dualMesh, 5./3., 1e-4, -1.0, "electron"), 
    ion_solver     (dualMesh, 5./3., 1.0,  1.0, "ion"),
    interpolator   (interpolator)
{
    init_A_and_D();
}

TwoFluidSolver::~TwoFluidSolver() {

}

void TwoFluidSolver::applyInitialCondition() {
    electron_solver.applyInitialCondition();
    ion_solver.applyInitialCondition();
}

void TwoFluidSolver::applyInitialCondition(const std::string h5_file) {
    electron_solver.applyInitialCondition(h5_file);
    ion_solver.applyInitialCondition(h5_file);
}

double TwoFluidSolver::updateFluxesExplicit() {
    const double dt_e = electron_solver.updateFluxExplicit();
    const double dt_i = ion_solver.updateFluxExplicit();
    // electron_solver.check_A_and_D(A, D);
    // ion_solver.check_A_and_D(A, D);
    return dt_e < dt_i ? dt_e : dt_i;
}

void TwoFluidSolver::updateMassFluxesImplicit() {
    electron_solver.updateMassFluxImplicit();
    ion_solver.updateMassFluxImplicit();
}

void TwoFluidSolver::updateMassFluxesImplicitLumped(const Eigen::VectorXd e, const Eigen::VectorXd dp, const Eigen::MatrixXd glb2lcl) {
    electron_solver.updateMassFluxImplicitLumped(e, dp, glb2lcl);
    ion_solver.updateMassFluxImplicitLumped(e, dp, glb2lcl);
}

void TwoFluidSolver::updateMomentum(const double dt, const Eigen::MatrixXd& E, const bool with_friction) {
    if (!with_friction) {
        electron_solver.updateMomentum(dt, E);
        ion_solver.updateMomentum(dt, E);
    }
    else {
        electron_solver.updateMomentum(dt, E, alpha, &ion_solver);
        ion_solver.updateMomentum(dt, E, alpha, &electron_solver);
    }
}

void TwoFluidSolver::updateRateOfChange(const bool with_rhs) {
     electron_solver.updateRateOfChange(with_rhs);
     ion_solver.updateRateOfChange(with_rhs);
     // electron_solver.check_eta();
     // ion_solver.check_eta();
}

void TwoFluidSolver::timeStepping(const double dt) {
    electron_solver.timeStepping(dt);
    ion_solver.timeStepping(dt);
}

void TwoFluidSolver::timeStepping(const double dt, const Eigen::MatrixXd& E, const Eigen::MatrixXd& B, const bool with_friction) {
    if (!with_friction) { // no friction
        electron_solver.timeStepping(dt, E, B);
        ion_solver.timeStepping(dt, E, B);
    }
    else { // with friction
        // NOTE: the state of both solvers must be updated simultaneously!!
        // TODO: wrap the timestepping into class FluidSolver
        // electron_solver.timeStepping(dt, E, B, alpha, &ion_solver);
        // ion_solver.timeStepping(dt, E, B, alpha, &electron_solver);
        electron_solver.rhs.setZero();
        electron_solver.applyLorentzForce(E, B);
        electron_solver.applyFrictionTerm(&ion_solver, alpha);
        electron_solver.updateRateOfChange(true);
        ion_solver.rhs.setZero();
        ion_solver.applyLorentzForce(E, B);
        ion_solver.applyFrictionTerm(&electron_solver, alpha);
        ion_solver.updateRateOfChange(true);
        electron_solver.U += dt * electron_solver.rate_of_change;
        ion_solver.U += dt * ion_solver.rate_of_change;

        // check 
        // electron_solver.check_updatedMomentum();
        // ion_solver.check_updatedMomentum();
        if (!electron_solver.isValidState() || !ion_solver.isValidState()) {
            std::cout << "******************************" << std::endl;
            std::cout << "*   Fluid State not valid!   *" << std::endl;
            std::cout << "******************************" << std::endl; 
        }
    }
}

void TwoFluidSolver::writeSnapshot(H5Writer &writer) const {
    electron_solver.writeSnapshot(writer);
    ion_solver.writeSnapshot(writer);
}

Eigen::SparseMatrix<double> TwoFluidSolver::get_M_sigma(const double dt, const bool with_friction, const bool lumpedElectricField) {
    
    Eigen::VectorXd dualFaceArea(dual->getNumberOfFaces());
    for (const Face* face : dual->getFaces()) {
        dualFaceArea[face->getIndex()] = face->getArea();
    }   
    if (!with_friction) {
        assert(lumpedElectricField == false); // Currently, collisionless case does not support lumpedElectricField
        auto T_e = electron_solver.get_T(dt, A, interpolator->get_E_interpolator());
        auto T_i = ion_solver.get_T     (dt, A, interpolator->get_E_interpolator());
        M_sigma = dualFaceArea.asDiagonal() * (electron_solver.charge * T_e + ion_solver.charge * T_i);
    }
    else if (with_friction && !lumpedElectricField) {
        auto T_e = electron_solver.get_T(dt, A, interpolator->get_E_interpolator(), alpha, &ion_solver);
        auto T_i = ion_solver.get_T     (dt, A, interpolator->get_E_interpolator(), alpha, &electron_solver);
        M_sigma = dualFaceArea.asDiagonal() * (electron_solver.charge * T_e + ion_solver.charge * T_i);
    }
    else {
        auto T_e = electron_solver.get_T(dt, _A, alpha, &ion_solver);
        auto T_i = ion_solver.get_T     (dt, _A, alpha, &electron_solver);
        Eigen::SparseMatrix<double> temp(dual->getNumberOfFaces(), dual->getNumberOfFaces());
        std::vector<T> triplets;
        for (int i = 0; i < dual->getNumberOfFaces(); i++) {
            const Face * f = dual->getFace(i);
            if (f->isBoundary()) {
                triplets.emplace_back(i, i, 1.0 / f->getArea());
            }
            else {
                triplets.emplace_back(i, i, 1.0 / primal->getEdge(i)->getLength());
            }
        }
        temp.setFromTriplets(triplets.begin(), triplets.end());
        M_sigma = dualFaceArea.asDiagonal() * (electron_solver.charge * T_e + ion_solver.charge * T_i) * temp;
    }
    // Attention: a small conductivity is added for the sake of stability?
    // M_sigma += Eigen::MatrixXd::Identity(dual->getNumberOfFaces(), dual->getNumberOfFaces()).sparseView() * 1e-4; 

    std::cout << "-- M_sigma assembled" << std::endl;
    
    return M_sigma;
}

Eigen::VectorXd TwoFluidSolver::get_j_aux(const double dt, const Eigen::MatrixXd& B, const bool with_friction) const {
    Eigen::VectorXd dualFaceArea(dual->getNumberOfFaces());
    for (const Face* face : dual->getFaces()) {
        dualFaceArea[face->getIndex()] = face->getArea();
    }

    Eigen::VectorXd mu_e = with_friction ? electron_solver.get_mu(dt, B, A, D, alpha, &ion_solver)
                                         : electron_solver.get_mu(dt, B, A, D);
    Eigen::VectorXd mu_i = with_friction ? ion_solver.get_mu(dt, B, A, D, alpha, &electron_solver)
                                         : ion_solver.get_mu(dt, B, A, D);

    Eigen::VectorXd j_aux = dualFaceArea.asDiagonal() * (electron_solver.charge * mu_e + ion_solver.charge * mu_i);
    std::cout << "- j_aux assembled" << std::endl;
    // j_aux.setZero();
    return j_aux;
}

void TwoFluidSolver::init_A_and_D() {
    A = Tensor3(dual->getNumberOfFaces(), dual->getNumberOfCells());
    _A.resize(dual->getNumberOfFaces(), dual->getNumberOfCells());
    D.resize(dual->getNumberOfFaces(), dual->getNumberOfCells());

    for (const Face* face : dual->getFluidFaces()) {
        const int face_idx = face->getIndex();
        const Eigen::Vector3d normal = face->getNormal();
        const std::vector<Cell*> cells = face->getCellList();
        switch (face->getFluidType()) {
            case Face::FluidType::Interior :
            {
                assert(cells.size() == 2);
                const Cell* leftCell  = cells[0]->getOrientation(face) > 0 ? cells[0] : cells[1];
                const Cell* rightCell = leftCell == cells[0] ? cells[1] : cells[0]; 
                A.insert(face_idx, cells[0]->getIndex(), normal);
                A.insert(face_idx, cells[1]->getIndex(), normal);
                _A.coeffRef(face_idx, cells[0]->getIndex()) = 1.0;
                _A.coeffRef(face_idx, cells[1]->getIndex()) = 1.0;
                D.coeffRef(face_idx, leftCell->getIndex())  = -1.0;
                D.coeffRef(face_idx, rightCell->getIndex()) =  1.0;
                break;
            }
            case Face::FluidType::Opening :
                assert(cells.size() == 1);
                A.insert(face_idx, cells[0]->getIndex(), 2 * normal);
                _A.coeffRef(face_idx, cells[0]->getIndex()) = 2.0;
                break;
            case Face::FluidType::Wall :
                // For the wall boundary, the mass flux is alway zero. (Note the momentum and energy fluxes are not zero)
                // assert(cells.size() == 2); If the fluid fills the whole domain, the wall face has only one adjacent cell.
                break; 
            case Face::FluidType::Mixed :
                assert(cells.size() == 1);
                for (const Face* subf : face->getSubFaceList()) {
                    if (subf->getFluidType() == Face::FluidType::Opening) {
                        A.insert(face_idx, cells[0]->getIndex(), 2*subf->getNormal() * subf->getArea() / face->getArea());
                    }
                }
                break;
            default:
                assert(false);
                break;
        }
    }
}

void TwoFluidSolver::checkChargeConservation(const double dt) {
    for (const Face* face : dual->getFaces()) {
        if (face->getFluidType() == Face::FluidType::Opening) {
            const double area = face->getArea();
            const int isInflow = face->getCenter().dot(face->getNormal()) < 0 ? 1 : -1;
            const int F_idx = electron_solver.face2F(face->getIndex()); 
            netChargeInflow += electron_solver.charge * electron_solver.F(F_idx, 0) * area * isInflow * dt;
            netChargeInflow += ion_solver.charge * ion_solver.F(F_idx, 0) * area * isInflow * dt;
        }
    }
    double chargeDiff = 0;
    for (const Cell* cell : dual->getFluidCells()) {
        const int U_idx = electron_solver.cell2U(cell->getIndex());
        const double vol = cell->getVolume();
        chargeDiff += electron_solver.charge * electron_solver.U(U_idx, 0) * vol;
        chargeDiff += ion_solver.charge * ion_solver.U(U_idx, 0) * vol;
    }
    std::cout << " ---------------- charge error = " << chargeDiff - netChargeInflow << std::endl;
}

std::pair<double, double> TwoFluidSolver::computeCurrent() const {
    double anodeCurrent = 0, cathodeCurrent = 0;
    double anodeArea = 0, cathodeArea = 0;
    for (const Face* face : dual->getFluidFaces()) {
        if (face->getFluidType() == Face::FluidType::Opening) {
            const double area = face->getArea();
            const int isInflow = face->getCenter().dot(face->getNormal()) < 0 ? 1 : -1;
            const int F_idx = electron_solver.face2F(face->getIndex());
            if (face->getCenter()[2] < 0) { // Cathode
                cathodeCurrent += electron_solver.charge * electron_solver.F(F_idx, 0) * area * isInflow;
                cathodeCurrent += ion_solver.charge * ion_solver.F(F_idx, 0) * area * isInflow;
                cathodeArea += area;
            }
            else { // Anode
                anodeCurrent += electron_solver.charge * electron_solver.F(F_idx, 0) * area * isInflow;
                anodeCurrent += ion_solver.charge * ion_solver.F(F_idx, 0) * area * isInflow;
                anodeArea += area;
            }
        }
    }
    return std::pair<double, double>{anodeCurrent / anodeArea, cathodeCurrent / cathodeArea};
}


