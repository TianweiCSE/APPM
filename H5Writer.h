#pragma once

#include "H5Cpp.h"
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
using namespace H5;


class H5Writer
{
public:
	H5Writer();
	H5Writer(const std::string & filename);
	~H5Writer();
	
	template <typename T, int cols, int rows>
	void writeData(const Eigen::Matrix<T, cols, rows> & matrix, const std::string & dataname);

	// Specialization for type int (-> no, modified)
	template <int cols, int rows> 
	void writeDataInt(const Eigen::Matrix<int,cols,rows> & matrix, const std::string & dataname);

	// Specialization for type double (-> no, modified)
	template <int cols, int rows>
	void writeDataDouble(const Eigen::Matrix<double, cols, rows> & matrix, const std::string & dataname);

	void writeData(const std::vector<int> & vector, const std::string & dataname);

private:
	std::string filename;
	H5File file;
};

template<typename T, int cols, int rows>
inline void H5Writer::writeData(const Eigen::Matrix<T, cols, rows>& matrix, const std::string & dataname)
{
	std::cout << "Generic function for Eigen matrices; do nothing" << std::endl;
}

inline void H5Writer::writeData(const std::vector<int>& vector, const std::string & dataname)
{
	std::cout << "Write int data to " << this->filename << ": " << dataname << std::endl;
	const int rank = 1;
	hsize_t dims[rank];
	dims[0] = vector.size();
	DataSpace dataspace(rank, dims);
	const PredType type = PredType::NATIVE_INT;
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(vector.data(), type);
}

// Specialization for type int
template <int cols, int rows>
inline void H5Writer::writeDataInt(const Eigen::Matrix<int, cols, rows> & matrix, const std::string & dataname)
{
	std::cout << "Write int data to " << this->filename << ": " << dataname << std::endl;
	const int rank = 2;
	hsize_t dims[rank];
	dims[1] = matrix.rows();
	dims[0] = matrix.cols();
	DataSpace dataspace(rank, dims);
	const PredType type = PredType::NATIVE_INT;
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(matrix.transpose().data(), type);
}

// Specialization for type double
template <int cols, int rows>
inline void H5Writer::writeDataDouble(const Eigen::Matrix<double, cols, rows> & matrix, const std::string & dataname)
{
	std::cout << "Write double data to " << this->filename << ": " << dataname << std::endl;
	const int rank = 2;
	hsize_t dims[rank];
	dims[1] = matrix.rows();
	dims[0] = matrix.cols();
	DataSpace dataspace(rank, dims);
	const PredType type = PredType::NATIVE_DOUBLE;
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(matrix.transpose().data(), type);
}

