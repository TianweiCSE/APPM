#pragma once
#include "FluidSolver.h"
#include "Interpolator.h"

typedef Eigen::Triplet<double> T;

class TwoFluidSolver
{
    public:
        TwoFluidSolver();
        TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh, const Interpolator* interpolator);
        ~TwoFluidSolver();

        void applyInitialCondition();

        const double updateFluxesExplicit();
        void updateFluxesImplicit(const Eigen::MatrixXd &E);
        void updateRateOfChange(const bool with_rhs = false);

        void timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B);

        void writeSnapshot(H5Writer &writer) const;

        
        Eigen::SparseMatrix<double>&& get_M_sigma(const double dt) const;
        Eigen::VectorXd&& get_j_aux(const double dt, const Eigen::MatrixXd &B) const;

    private:
        const PrimalMesh* primal;
        const DualMesh* dual;

        FluidSolver electron_solver;
        FluidSolver ion_solver;
        const Interpolator* interpolator;

        Tensor3* A = nullptr;           //< see definition in (4.39)
        Eigen::SparseMatrix<double> D;  //< see definition in (4.39)

        void init_A_and_D();

};
