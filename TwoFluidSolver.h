#pragma once
#include "FluidSolver.h"
#include "Interpolator.h"

typedef Eigen::Triplet<double> T;

class TwoFluidSolver
{
    public:
        TwoFluidSolver();
        /**
         * @brief Construct a new two-Fluid solver object
         * 
         * @param primalMesh primal mesh
         * @param dualMesh dual mesh
         * @param interpolator object responsible for computing E- and B-fields at dual cell centers
         */
        TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh, const Interpolator* interpolator);
        ~TwoFluidSolver();

        void applyInitialCondition();

        /**
         * @brief Update the explicit fluxes at dual fluid faces for each species
         * 
         * @return time step size restricted by CFL condition
         */
        const double updateFluxesExplicit();
        /**
         * @brief Update the semi-implicit fluxes at dual fluid faces for each species 
         * 
         * @param dt time step size
         * @param E New E-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
         */
        void updateMassFluxesImplicit(const double dt, const Eigen::MatrixXd E);
        /**
         * @brief Update the rate of change of conservative variables for each species
         * 
         * @param with_rhs True if source terms are considered 
         */
        void updateRateOfChange(const bool with_rhs);

        void timeStepping(const double dt);
        /**
         * @brief Do timestepping for each species
         * 
         * @param dt time step size
         * @param E New E-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
         * @param B Old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
         */
        void timeStepping(const double dt, const Eigen::MatrixXd E, const Eigen::MatrixXd B);
        /**
         * @brief Write snapshot of solutions for each species
         * 
         * @param writer h5writer
         */
        void writeSnapshot(H5Writer &writer) const;

        
        Eigen::SparseMatrix<double> get_M_sigma(const double dt) const;
        Eigen::VectorXd get_j_aux(const double dt, const Eigen::MatrixXd&& B) const;

    private:
        const PrimalMesh* primal;
        const DualMesh* dual;

        FluidSolver electron_solver;
        FluidSolver ion_solver;
        const Interpolator* interpolator;

        Tensor3 A;           //< see definition in (4.39)
        Eigen::SparseMatrix<double> D;  //< see definition in (4.39)

        void init_A_and_D();

};
