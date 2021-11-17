#include "Interpolator.h"

Interpolator::Interpolator() {

}

Interpolator::Interpolator(const PrimalMesh* primal, const DualMesh* dual)
: primal(primal), dual(dual) {

}

Interpolator::~Interpolator() {

}

Eigen::SparseMatrix<double>&& Interpolator::getInterpolateMat_E() {

    Eigen::SparseMatrix<double> R(3,3);

    return std::move(R);
}

Eigen::SparseMatrix<double>&& Interpolator::getInterpolateMat_B() {
    Eigen::SparseMatrix<double> R(3,3);

    return std::move(R);
}