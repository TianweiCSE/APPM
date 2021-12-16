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

const double TwoFluidSolver::updateFluxesExplicit() {
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

Eigen::SparseMatrix<double> TwoFluidSolver::get_M_sigma(const double dt, const bool with_friction) const {
    Eigen::VectorXd dualFaceArea(dual->getNumberOfFaces());
    for (const Face* face : dual->getFaces()) {
        dualFaceArea[face->getIndex()] = face->getArea();
    }
    auto T_e = with_friction ? electron_solver.get_T(dt, A, interpolator->get_E_interpolator(), alpha, &ion_solver)
                             : electron_solver.get_T(dt, A, interpolator->get_E_interpolator());
    auto T_i = with_friction ? ion_solver.get_T     (dt, A, interpolator->get_E_interpolator(), alpha, &electron_solver)
                             : ion_solver.get_T     (dt, A, interpolator->get_E_interpolator());

    Eigen::SparseMatrix<double> M_sigma = 
        dualFaceArea.asDiagonal() * (electron_solver.charge * T_e + ion_solver.charge * T_i);

    // Attention: a small conductivity is added for the sake of stability?
    // M_sigma += Eigen::MatrixXd::Identity(dual->getNumberOfFaces(), dual->getNumberOfFaces()).sparseView() * 1e-4; 

    std::cout << "- M_sigma assembled" << std::endl;
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
    return j_aux;
}

void TwoFluidSolver::init_A_and_D() {
    A = Tensor3(dual->getNumberOfFaces(), dual->getNumberOfCells());
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
                D.coeffRef(face_idx, leftCell->getIndex())  = -1.0;
                D.coeffRef(face_idx, rightCell->getIndex()) =  1.0;
                break;
            }
            case Face::FluidType::Opening :
                assert(cells.size() == 1);
                A.insert(face_idx, cells[0]->getIndex(), 2 * normal);
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


