#include "Tensor3.h"

Tensor3::Tensor3() {

}

Tensor3::Tensor3(const int size_first, const int size_last) 
: size_first(size_first), size_last(size_last) {

    data = new Eigen::SparseMatrix<double>[size_first];
    for (int i = 0; i < size_first; i++) {
        data[i].resize(3, size_last);
    }
}

Tensor3::~Tensor3() {
    delete[] data;
}

Eigen::SparseMatrix<double>& Tensor3::operator[](const int row) const{
    return data[row];
}

Tensor3&& Tensor3::operator=(const Tensor3& t) const {
    Tensor3 new_t(t.size_first, t.size_last);
    for (int i = 0; i < size_first; i++) {
        new_t.data[i] = t.data[i];
    }
    return std::move(new_t);
}

void Tensor3::insert(const int first_idx, const int last_idx, const Eigen::Vector3d v) {
    data[first_idx].coeffRef(0, last_idx) = v[0];
    data[first_idx].coeffRef(1, last_idx) = v[1];
    data[first_idx].coeffRef(2, last_idx) = v[2];
}

Eigen::MatrixX3d&& Tensor3::oneContract(const Eigen::VectorXd& vec) const {
    assert(vec.size() == size_last);
    Eigen::MatrixX3d output(size_first, 3);
    for (int i = 0; i < size_first; i++) {
        output.row(i) = data[i] * vec;
    }
    return std::move(output);
}
    
Eigen::VectorXd&& Tensor3::twoContract(const Eigen::MatrixX3d& mat) const {
    assert(mat.rows() == size_last);
    Eigen::VectorXd output(size_first);
    for (int i = 0; i < size_first; i++) {
        output[i] = data[i].cwiseProduct(mat.transpose()).sum();
    }
    return std::move(output);
}

Eigen::SparseMatrix<double>&& Tensor3::twoContract(const Tensor3& tensor) const {
    assert(tensor.size_first == size_last);
    Eigen::SparseMatrix<double> output(size_first, tensor.size_last);
    std::vector<T> triplets;
    for (int i = 0; i < size_first; i++) {
        for (int j = 0; j < tensor.size_last; j++) {
            Eigen::MatrixX3d temp(size_last, 3);
            for (int k = 0; k < size_last; k++) {
                temp.row(k) = tensor[k].col(j);
            }
            const double value = data[i].cwiseProduct(temp).sum(); 
            if (std::abs(value) > std::numeric_limits<double>::epsilon()) {
                triplets.emplace_back(T{i, j, value});
            }
        }
    }
    output.setFromTriplets(triplets.begin(), triplets.end());
    output.makeCompressed();
    return std::move(output);
}

Tensor3&& Tensor3::firstDimWiseProduct(const Eigen::VectorXd& vec) const {
    Tensor3 output = (*this);
    for (int i = 0; i < size_first; i++) {
        output[i] *= vec[i];
    }
    return std::move(output);
}

