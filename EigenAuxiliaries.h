#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>

namespace Eigen {
	template <typename T>
	std::vector<Eigen::Triplet<T>> sparseMatrixToTriplets(const SparseMatrix<T> & M) {
		std::vector<Eigen::Triplet<T>> triplets;
		for (int i = 0; i < M.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it) {
				triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
			}
		}
		return triplets;
	}

	template <typename T>
	void sparseMatrixToFile(const SparseMatrix<T> & M, const std::string & filename) {
		std::cout << "Write matrix to file: " << filename << std::endl;
		std::ofstream file(filename);
		for (int i = 0; i < M.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it) {
				file << it.row() << "," << it.col() << "," << it.value() << std::endl;
			}
		}
	}



	/** 
	* Get a sparse identitiy matrix of non-square shape. (See also speye() in Matlab)
	*/
	//template <typename T>
	//Eigen::SparseMatrix<double> speye(const int rows, const int cols) {
	//	assert(rows > 0);
	//	assert(cols > 0);
	//	const int n = std::min(rows, cols);

	//	typedef Eigen::Triplet<double> TripletType;
	//	std::vector<TripletType> triplets;
	//	triplets.reserve(n);
	//	for (int i = 0; i < n; i++) {
	//		triplets.emplace_back(TripletType(i, i, 1));
	//	}
	//	Eigen::SparseMatrix<double> M(rows, cols);
	//	M.setFromTriplets(triplets.begin(), triplets.end());
	//	return M;
	//}


	/** 
	* Concatenate two sparse matrices vertically.
	*/
	template <typename T>
	Eigen::SparseMatrix<T> vertcat(const Eigen::SparseMatrix<T> & U, const Eigen::SparseMatrix<T> & L) {
		assert(U.cols() == L.cols());
		
		typedef Eigen::Triplet<T> TripletType;
		std::vector<TripletType> triplets;
		triplets.reserve(U.nonZeros() + L.nonZeros());

		// get triplets from upper matrix
		

		for (int k = 0; k < U.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(U, k); it; ++it) {
				triplets.emplace_back(it.row(), it.col(), it.value());
			}
		}
		// get triplets from lower matrix
		for (int k = 0; k < L.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(L, k); it; ++it) {
				triplets.emplace_back(U.rows() + it.row(), it.col(), it.value());
			}
		}
		// set matrix from triplets
		Eigen::SparseMatrix<T> M(U.rows() + L.rows(), U.cols());
		M.setFromTriplets(triplets.begin(), triplets.end());
		return M;
	}

	template <typename T>
	void blockFill(std::vector<Eigen::Triplet<T>> &triplets, const int start_row, const int start_col, Eigen::SparseMatrix<T> &&block) {
		for (int i = 0; i < block.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(block, i); it; ++it) {
				triplets.emplace_back(it.row() + start_row, it.col() + start_col, it.value());
			}
		}
	}

	template <typename T>
	void blockFill(std::vector<Eigen::Triplet<T>> &triplets, const int start_row, const int start_col, const Eigen::SparseMatrix<T> &block) {
		for (int i = 0; i < block.outerSize(); i++) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(block, i); it; ++it) {
				triplets.emplace_back(it.row() + start_row, it.col() + start_col, it.value());
			}
		}
	}

	template <typename T>
	void countRowNNZ(const Eigen::SparseMatrix<T> &M) {
		int avg = M.nonZeros() / M.rows();
		const int rows = M.rows();
		const int cols = M.cols();
		int maxRowNNZ = 0;
		for (int i = 0; i < rows; i++) {
			const Eigen::VectorXd rowVec = M.row(i);
			int count = 0;
			for (int j = 0; j < cols; j++) {
				if (std::abs(M.coeff(i,j)) > 1e-10) {
					count++;
				}
			}
			if (count > maxRowNNZ) {
				maxRowNNZ = count;
			}
		}
		std::cout << "avgNNz = " << avg << ", maxNNz = " << maxRowNNZ << std::endl;
		return;
	}

	// Function to save an Eigen matrix to a binary file
	template<typename Matrix>
	void saveMatrix(const std::string& filename, const Matrix& matrix) {
		std::ofstream file(filename, std::ios::out | std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: Unable to open the file for writing.\n";
			return;
		}
		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
		file.write(reinterpret_cast<const char*>(&rows), sizeof(typename Matrix::Index));
		file.write(reinterpret_cast<const char*>(&cols), sizeof(typename Matrix::Index));
		file.write(reinterpret_cast<const char*>(matrix.data()), rows * cols * sizeof(typename Matrix::Scalar));
		file.close();
	}

	// Function to load an Eigen matrix from a binary file
	template<typename Matrix>
	void loadMatrix(const std::string& filename, Matrix& matrix) {
		std::ifstream file(filename, std::ios::in | std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: Unable to open the file for reading.\n";
			return;
		}
		typename Matrix::Index rows, cols;
		file.read(reinterpret_cast<char*>(&rows), sizeof(typename Matrix::Index));
		file.read(reinterpret_cast<char*>(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		file.read(reinterpret_cast<char*>(matrix.data()), rows * cols * sizeof(typename Matrix::Scalar));
		file.close();
	}


}

