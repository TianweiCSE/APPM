#pragma once
#include <Eigen/Sparse>
#include <numeric>

typedef Eigen::Triplet<double> T;

class Tensor3 {
    public:
        Tensor3();
        Tensor3(const int size_first, const int size_last);
        ~Tensor3();

        Eigen::SparseMatrix<double>& operator[](const int row) const;
        Tensor3&&                    operator=(const Tensor3& t) const;
        // Tensor3&&                 operator*(const double c);

        void insert(const int first_idx, const int last_idx, const Eigen::Vector3d v);

        Eigen::MatrixX3d&& oneContract(const Eigen::VectorXd&  vec) const;
        Eigen::VectorXd&&  twoContract(const Eigen::MatrixX3d& mat) const;
        Eigen::SparseMatrix<double>&& twoContract(const Tensor3& tensor) const;
        Tensor3&& firstDimWiseProduct(const Eigen::VectorXd& vec) const;
        
    private:
        const int size_first = 0;
        const int size_last  = 0;
        Eigen::SparseMatrix<double>* data = nullptr; 

};