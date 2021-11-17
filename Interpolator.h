#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include <Eigen/Dense>

class Interpolator {
    public:
        Interpolator();
        Interpolator(const PrimalMesh* primal, const DualMesh* dual);
        ~Interpolator();
        Eigen::SparseMatrix<double>&& getInterpolateMat_E();
        Eigen::SparseMatrix<double>&& getInterpolateMat_B();
    private:
        const PrimalMesh* primal = nullptr;
        const DualMesh* dual = nullptr;

};
