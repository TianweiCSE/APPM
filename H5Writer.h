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
	
	void writeStdVector(const std::vector<int> & vector, const std::string & dataname);
	void writeIntMatrix(const Eigen::MatrixXi & matrix, const std::string & dataname);
	void writeDoubleMatrix(const Eigen::MatrixXd & matrix, const std::string & dataname);
	void writeIntVector(const Eigen::VectorXi & vector, const std::string & dataname);
	void writeDoubleVector(const Eigen::VectorXd & vector, const std::string & dataname);

private:
	std::string filename;
	H5File file;
};






