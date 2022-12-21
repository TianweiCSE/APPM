#pragma once
#include "FluidSolver.h"
#include "Interpolator.h"

typedef Eigen::Triplet<double> T;
class AppmSolver;
class MaxwellSolver;

/**
 * @brief This class is responsible for tackling Euler part in two-fluid plasma model.
 * 
 * This class possesses two basic fluid solvers by which Euler equation is solved.
 * The main task for this class:
 *      - compute the Ohm's law (J = E_sigma * e + J_aux).
 *      - compute the tensor A and metrix D defined in (4.39).
 */
class TwoFluidSolver
{
    public:
        TwoFluidSolver() = delete;
        /**
         * @brief Construct a new two-Fluid solver object
         * 
         * @param primalMesh primal mesh
         * @param dualMesh dual mesh
         * @param interpolator object responsible for computing E- and B-fields at dual cell centers
         */
        TwoFluidSolver(const PrimalMesh* primalMesh, const DualMesh* dualMesh, const Interpolator* interpolator);
        ~TwoFluidSolver();

        /**
         * @brief Assign initial conditions for all electromagnetic variables.
         * 
         */
        void applyInitialCondition();
        void applyInitialCondition(const std::string h5_file);

        /**
         * @brief Update the explicit fluxes at dual fluid faces for each species. 
         * 
         * Routine <updateFluxExplicit> of two child fluid solvers are called. 
         * 
         * @return time step size restricted by CFL condition
         */
        double updateFluxesExplicit();
        /**
         * @brief Update the semi-implicit fluxes at dual fluid faces for each species 
         * 
         * Routine <updateMassFluxImplicit> of two child fluid solvers are called.
         */
        void updateMassFluxesImplicit();
        void updateMassFluxesImplicitLumped(const Eigen::VectorXd& e, const Eigen::VectorXd& dp, const Eigen::SparseMatrix<double>& glb2lcl);

        void updateMomentum(const double dt, const Eigen::MatrixXd& E, const bool with_friction);

        /**
         * @brief Update the rate of change of conservative variables for each species
         * 
         * @param with_rhs True if source terms are considered 
         */
        void updateRateOfChange(const bool with_rhs);

        /**
         * @brief Evolve the conservative fluid variables for each species. Source term is ignored.
         * 
         * @param dt 
         */
        void timeStepping(const double dt);
        /**
         * @brief Evolve the conservative fluid variables for each species. Lorentz force is included.
         * 
         * @param dt time step size
         * @param E New E-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
         * @param B Old B-field defined at each cell center. Entries are indexed by cells with entries of solid cells being ZERO.
         * @param with_friction true if friction term is included
         */
        void timeStepping(const double dt, const Eigen::MatrixXd &E, const Eigen::MatrixXd &B, const bool with_friction);

        /**
         * @brief Write snapshot of solutions for each species
         * 
         * @param writer h5writer
         */
        void writeSnapshot(H5Writer &writer) const;

        /**
         * @brief Compute M_sigma defined in (4.41)
         * 
         * @param dt time step size
         * @param with_friction true if friction is included
         * @param lumpedElectricField true if the current flux is computed by lumped electric field
         * @return M_sigma 
         */
        Eigen::SparseMatrix<double> get_M_sigma(const double dt, const bool with_friction, const bool lumpedElectricField);

        /**
         * @brief Compute j_aux defined in (4.41)
         * 
         * @param dt time step size
         * @param B Intepolated B-field at dual cell center
         * @param with_friction true if friction is included
         * @return j_aux
         */
        Eigen::VectorXd get_j_aux(const double dt, const Eigen::MatrixXd& B, const bool with_friction) const;

        void checkChargeConservation(const double dt);
        std::pair<double, double> computeCurrent() const;

    private:
        const PrimalMesh* primal;
        const DualMesh* dual;

        FluidSolver electron_solver;
        FluidSolver ion_solver;
        const Interpolator* interpolator;

        double alpha;

        Tensor3 A;                      //< see definition in (4.39)
        Eigen::SparseMatrix<double> _A; // test
        Eigen::SparseMatrix<double> D;  //< see definition in (4.39)
        Eigen::SparseMatrix<double> M_sigma;

        double netChargeInflow = 0;

        void init_A_and_D();

        friend class AppmSolver;
        friend class MaxwellSolver;

};
