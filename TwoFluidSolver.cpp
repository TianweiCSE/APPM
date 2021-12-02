#include "TwoFluidSolver.h"

TwoFluidSolver::TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh, const Interpolator* interpolator) :
    primal         (primalMesh),
    dual           (dualMesh),
    electron_solver(dualMesh, 5./3., 1e-4,  -1.0, "electron"), 
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

const double TwoFluidSolver::updateFluxesExplicit() {
    const double dt_e = electron_solver.updateFluxExplicit();
    const double dt_i = ion_solver.updateFluxExplicit();
    return dt_e < dt_i ? dt_e : dt_i;
}

void TwoFluidSolver::updateMassFluxesImplicit(const double dt, const Eigen::MatrixXd E) {
    electron_solver.updateMassFluxImplicit(dt, E);
    ion_solver.updateMassFluxImplicit(dt, E);
}

void TwoFluidSolver::updateRateOfChange(const bool with_rhs) {
     electron_solver.updateRateOfChange(with_rhs);
     ion_solver.updateRateOfChange(with_rhs);
}

void TwoFluidSolver::timeStepping(const double dt) {
    electron_solver.timeStepping(dt);
    ion_solver.timeStepping(dt);
}

void TwoFluidSolver::timeStepping(const double dt, const Eigen::MatrixXd E, const Eigen::MatrixXd B) {
    
    double max_E = 0, max_B = 0;
    for (int i = 0; i < dual->getNumberOfCells(); i++) {
        if (dual->getCell(i)->getFluidType() == Cell::FluidType::Fluid) {
            max_E = E.row(i).norm() > max_E ? E.row(i).norm() : max_E;
            max_B = B.row(i).norm() > max_B ? B.row(i).norm() : max_B;
        }
    }
    std::cout << " ---------- max E-field :" << E.rowwise().norm().maxCoeff() << std::endl;
    std::cout << " ---------- max E-field in fluid:" << max_E << std::endl;
    std::cout << " ---------- max B-field :" << B.rowwise().norm().maxCoeff() << std::endl;
    std::cout << " ---------- max B-field in fluid:" << max_B << std::endl;
    electron_solver.timeStepping(dt, E, B);
    ion_solver.timeStepping(dt, E, B);
}

void TwoFluidSolver::writeSnapshot(H5Writer &writer) const {
    electron_solver.writeSnapshot(writer);
    ion_solver.writeSnapshot(writer);
}


Eigen::SparseMatrix<double> TwoFluidSolver::get_M_sigma(const double dt) const {
    Eigen::VectorXd dualFaceArea(dual->getNumberOfFaces());
    for (const Face* face : dual->getFaces()) {
        dualFaceArea[face->getIndex()] = face->getArea();
    }
    auto T_e = electron_solver.get_T(dt, A, interpolator->get_E_interpolator());
    auto T_i = ion_solver.get_T     (dt, A, interpolator->get_E_interpolator());

    Eigen::SparseMatrix<double> M_sigma = 
        dualFaceArea.asDiagonal() * (electron_solver.charge * T_e + ion_solver.charge * T_i);
    std::cout << "- M_sigma assembled" << std::endl;
    return M_sigma;
}

Eigen::VectorXd TwoFluidSolver::get_j_aux(const double dt, const Eigen::MatrixXd&& B) const {
    Eigen::VectorXd dualFaceArea(dual->getNumberOfFaces());
    for (const Face* face : dual->getFaces()) {
        dualFaceArea[face->getIndex()] = face->getArea();
    }

    Eigen::VectorXd mu_e = electron_solver.get_mu(dt, B, A, D);
    Eigen::VectorXd mu_i = ion_solver.get_mu(dt, B, A, D);

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
                assert(cells.size() == 2);
                break; 
            default:
                assert(false);
                break;
        }
    }
}


