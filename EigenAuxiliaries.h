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
}

