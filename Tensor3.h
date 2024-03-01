#pragma once
#include <Eigen/Sparse>
#include <numeric>
#include <iostream>

typedef Eigen::Triplet<double> T;

/**
 * @brief A data structure to tackle tensors of rank three. The second dimension is restricted to three.
 * 
 */
class Tensor3 {
    public:

        Tensor3();
        Tensor3(const int size_first, const int size_last);
        Tensor3(const Tensor3& other);
        Tensor3(Tensor3&& other) noexcept;
        ~Tensor3();

        //Eigen::SparseMatrix<double>& operator[](const int row) const;
        Tensor3& operator=(const Tensor3& other);
        Tensor3& operator=(Tensor3&& other) noexcept;
        Tensor3& operator*(const double c);

        void insert(const int first_idx, const int last_idx, const Eigen::Vector3d v);
        void setFromTriplets(std::vector<T> comp1_trip, std::vector<T> comp2_trip, std::vector<T> comp3_trip);

        Eigen::MatrixXd oneContract(const Eigen::VectorXd& vec) const;
        Eigen::MatrixXd oneContract(Eigen::VectorXd&& vec) const;
        Eigen::VectorXd twoContract(const Eigen::MatrixX3d& mat) const;
        Eigen::VectorXd twoContract(Eigen::MatrixX3d&& mat) const;
        Eigen::SparseMatrix<double> twoContract(const Tensor3& tensor) const;
        Eigen::SparseMatrix<double> twoContract(Tensor3&& tensor) const;
        Tensor3 firstDimWiseProduct(const Eigen::VectorXd& vec) const;
        Tensor3 firstDimWiseProduct(Eigen::VectorXd&& vec) const;
        
    private:
        int size_first = 0;
        int size_last  = 0;
        Eigen::SparseMatrix<double>* data = nullptr; 

};