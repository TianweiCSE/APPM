#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include <Eigen/Dense>
#include "Tensor3.h"

typedef std::tuple<int, int, Eigen::Vector3d> tuple3;

class Interpolator {

    public:
        Interpolator();
        Interpolator(const PrimalMesh* primal, const DualMesh* dual);
        ~Interpolator();

        // Get the tensor that interpolates [e, d^partial]^T into E-field at each fluid cell center.
        // Each tensor contains (fluid cell idx, vector component idx, 3d vector)
        const Tensor3& get_E_interpolator() const;
        // Get the tensor that interpolates [e, d^partial]^T into E-field at each fluid cell center.
        // Each tensor contains (fluid cell idx, vector component idx, 3d vector)
        const Tensor3& get_B_interpolator() const; 
        
    private:
        const PrimalMesh* primal = nullptr;
        const DualMesh* dual     = nullptr;
        Tensor3* E_interpolator = nullptr;
        Tensor3* B_interpolator = nullptr;

        // Init the tensor that interpolates [e, d^partial]^T into E-field at each fluid cell center.
        // Each tensor contains (fluid cell idx, vector component idx, 3d vector)
        // Note: the first index corresponds to the fluid cell idx, not fluid variable idx.
        void initInterpolateTensor_E();
        // Init the tensor that interpolates [e, d^partial]^T into E-field at each fluid cell center.
        // Each tensor contains (fluid cell idx, vector component idx, 3d vector)
        // Note: the first index corresponds to the fluid cell idx, not fluid variable idx.
        void initInterpolateTensor_B();

};
