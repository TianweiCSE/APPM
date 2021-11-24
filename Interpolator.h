#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include <Eigen/Dense>
#include "Tensor3.h"

class Interpolator {

    public:
        Interpolator();
        Interpolator(const PrimalMesh* primal, const DualMesh* dual);
        ~Interpolator();

        // Get the tensor that interpolates [e, d^partial]^T into E-field vector at each fluid cell center.
        const Tensor3& get_E_interpolator() const;
        // Get the tensor that interpolates [e, d^partial]^T into B-field vector at each fluid cell center.
        const Tensor3& get_B_interpolator() const; 
        
    private:
        const PrimalMesh* primal = nullptr;
        const DualMesh* dual     = nullptr;
        Tensor3 E_interpolator;
        Tensor3 B_interpolator;

        // Init the tensor that interpolates [e, d^partial]^T into E-field vector at each fluid cell center.
        // Note: the first index corresponds to the fluid cell idx, not fluid variable idx.
        void initInterpolateTensor_E();
        // Init the tensor that interpolates [e, d^partial]^T into B-field vector at each fluid cell center.
        // Note: the first index corresponds to the fluid cell idx, not fluid variable idx.
        void initInterpolateTensor_B();

};
