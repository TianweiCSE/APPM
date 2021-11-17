#pragma once
#include "FluidSolver.h"
#include "Interpolator.h"


class TwoFluidSolver
{
    public:
        TwoFluidSolver();
        TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh);
        ~TwoFluidSolver();

        void applyInitialCondition();

        const double updateFluxesExplicit();
        void updateFluxesImplicit(const Eigen::MatrixXd &E);
        void updateRateOfChange(const bool with_rhs = false);

        void timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);

        void writeSnapshot(H5Writer &writer) const;

        Eigen::SparseMatrix<double>&& get_M_sigma(const double dt, const Eigen::MatrixXd &B) const;
        Eigen::VectorXd&& get_j_aux(const double dt) const;

    private:
        FluidSolver electron_solver;
        FluidSolver ion_solver;
        Interpolator interpolator;

};
