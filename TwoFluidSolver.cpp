#include "TwoFluidSolver.h"

TwoFluidSolver::TwoFluidSolver() {

}

TwoFluidSolver::TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh) :
    electron_solver(dualMesh, 5./3., 1e-4, -1.0, "electron"), 
    ion_solver     (dualMesh, 5./3., 1.0,   1.0, "ion"),
    interpolator   (primalMesh, dualMesh) {

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

void TwoFluidSolver::updateFluxesImplicit(const Eigen::MatrixXd &E) {
    electron_solver.updateFluxImplicit(E);
    ion_solver.updateFluxImplicit(E);
}

 void TwoFluidSolver::updateRateOfChange(const bool with_rhs) {
     electron_solver.updateRateOfChange(with_rhs);
     ion_solver.updateRateOfChange(with_rhs);
 }


void TwoFluidSolver::timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B) {
    electron_solver.timeStepping(dt, E, B);
    ion_solver.timeStepping(dt, E, B);
}

void TwoFluidSolver::writeSnapshot(H5Writer &writer) const {
    electron_solver.writeSnapshot(writer);
    ion_solver.writeSnapshot(writer);
}


Eigen::SparseMatrix<double>&& TwoFluidSolver::get_M_sigma(const double dt, const Eigen::MatrixXd &B) const {
    Eigen::SparseMatrix<double> M_sigma(3,3);

    return std::move(M_sigma);
}

Eigen::VectorXd&& TwoFluidSolver::get_j_aux(const double dt) const {
    Eigen::VectorXd j_aux(3);

    return std::move(j_aux);
}

