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
	void writeData(const Eigen::Matrix<T,cols,rows> & matrix, const std::string & dataname);

	void writeDataVector(const std::vector<int> & vec, const std::string & dataname);

private:
	std::string filename;
	H5File file;
};

template <typename T, int cols, int rows>
inline void H5Writer::writeData(const Eigen::Matrix<T, cols, rows> & matrix, const std::string & dataname)
{
	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
	if (matrix.rows() > 1){
		const int rank = 2;
		hsize_t dims[rank];
		dims[0] = matrix.cols();
		dims[1] = matrix.rows();
	}
	else{
		const int rank = 1;
		hsize_t dims[rank];
		dims[0] = matrix.cols();
	}
	DataSpace dataspace(rank, dims);

	if std::numeric_limits<T>::is_integer{
		const PredType type = PredType::NATIVE_INT;
	}
	else{
		const PredType type = PredType::NATIVE_DOUBLE;
	}
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(matrix.transpose().data(), type);
}

void H5Writer::writeDataVector(const std::vector<int> & vec, const std::string & dataname){
	std::cout << "Write int data to " << this->filename << ": " << dataname << std::endl;
	const int rank = 1;
	hsize_t dims[rank];
	dims[0] = vec.size();
	DataSpace dataspace(rank, dims);
	const PredType type = PredType::NATIVE_INT;
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(vec.data(), type);
}


