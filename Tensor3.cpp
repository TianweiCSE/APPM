#include "Tensor3.h"

Tensor3::Tensor3() {

}

Tensor3::Tensor3(const int size_first, const int size_last) 
: size_first(size_first), size_last(size_last) {
    data = new Eigen::SparseMatrix<double>[3];
    for (int i = 0; i < 3; i++) {
        data[i].resize(size_first, size_last);
    }
}

Tensor3::Tensor3(const Tensor3& other) {
    size_first = other.size_first;
    size_last  = other.size_last;
    data = new Eigen::SparseMatrix<double>[3];
    for (int i = 0; i < 3; i++) {
        data[i] = other.data[i];
    }
}

Tensor3::Tensor3(Tensor3&& other) noexcept {
    size_first = other.size_first;
    size_last  = other.size_last;
    data = other.data;
    other.data = nullptr;
}

Tensor3::~Tensor3() {
    delete[] data;
}

Tensor3& Tensor3::operator=(const Tensor3& other) {
    if (this == &other) return *this;
    delete[] data;
    this->size_first = other.size_first;
    this->size_last  = other.size_last;
    this->data = new Eigen::SparseMatrix<double>[3];
    for (int i = 0; i < 3; i++) {
        data[i] = other.data[i];
    }
    return *this;
}

Tensor3& Tensor3::operator=(Tensor3&& other) noexcept {
    if (this == &other) return *this;
    delete[] data;
    this->size_first = other.size_first;
    this->size_last  = other.size_last;
    this->data = other.data;
    other.data = nullptr;
    return *this;
}

Tensor3& Tensor3::operator*(const double c) {
    for (int i = 0; i < 3; i++) {
        data[i] *= c;
    }
    return *this;
}

void Tensor3::insert(const int first_idx, const int last_idx, const Eigen::Vector3d v) {
    assert(first_idx < size_first && last_idx < size_last);
    data[0].coeffRef(first_idx, last_idx) = v[0];
    data[1].coeffRef(first_idx, last_idx) = v[1];
    data[2].coeffRef(first_idx, last_idx) = v[2];
}

Eigen::MatrixXd Tensor3::oneContract(const Eigen::VectorXd& vec) const {
    assert(vec.size() == size_last);
    Eigen::MatrixXd output(size_first, 3);
    output.col(0) = data[0] * vec;
    output.col(1) = data[1] * vec;
    output.col(2) = data[2] * vec;
    return output;
}

Eigen::MatrixXd Tensor3::oneContract(Eigen::VectorXd&& vec) const {
    assert(vec.size() == size_last);
    Eigen::MatrixXd output(size_first, 3);
    output.col(0) = data[0] * vec;
    output.col(1) = data[1] * vec;
    output.col(2) = data[2] * vec;
    return output;
}
    
Eigen::VectorXd Tensor3::twoContract(const Eigen::MatrixX3d& mat) const {
    assert(mat.rows() == size_last);
    Eigen::VectorXd output(size_first);
    output.setZero();
    output += data[0] * mat.col(0);
    output += data[1] * mat.col(1);
    output += data[2] * mat.col(2);
    return output;
}

Eigen::VectorXd Tensor3::twoContract(Eigen::MatrixX3d&& mat) const {
    assert(mat.rows() == size_last);
    Eigen::VectorXd output(size_first);
    output.setZero();
    output += data[0] * mat.col(0);
    output += data[1] * mat.col(1);
    output += data[2] * mat.col(2);
    return output;
}

Eigen::SparseMatrix<double> Tensor3::twoContract(const Tensor3& tensor) const {
    assert(tensor.size_first == size_last);
    Eigen::SparseMatrix<double> output(size_first, tensor.size_last);
    output = data[0] * tensor.data[0] 
           + data[1] * tensor.data[1]
           + data[2] * tensor.data[2]; 
    output.makeCompressed();
    return output;
}

Eigen::SparseMatrix<double> Tensor3::twoContract(Tensor3&& tensor) const {
    assert(tensor.size_first == size_last);
    Eigen::SparseMatrix<double> output(size_first, tensor.size_last);
    output = data[0] * tensor.data[0] 
           + data[1] * tensor.data[1]
           + data[2] * tensor.data[2]; 
    output.makeCompressed();
    return output;
}

Tensor3 Tensor3::firstDimWiseProduct(const Eigen::VectorXd& vec) const {
    assert(vec.size() == size_first);
    Tensor3 output = (*this);
    for (int i = 0; i < 3; i++) {
        output.data[i] = vec.asDiagonal() * output.data[i];
    }
    return output;
}

Tensor3 Tensor3::firstDimWiseProduct(Eigen::VectorXd&& vec) const {
    assert(vec.size() == size_first);
    Tensor3 output = (*this);
    for (int i = 0; i < 3; i++) {
        output.data[i] = vec.asDiagonal() * output.data[i];
    }
    return output;
}
 