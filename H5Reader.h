#pragma once

#include "H5Cpp.h"
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
using namespace H5;

class H5Reader
{
public:
	H5Reader();
	H5Reader(const std::string & filename);
	~H5Reader();

	int readVectorDataSize(const std::string & dataname) const;
	int readMatrixDataSize(const std::string & dataname) const;
	Eigen::VectorXd readVectorData(const std::string & dataname) const;
	Eigen::MatrixXd readMatrixData(const std::string & dataname) const;

private:
	std::string filename;
	H5File file;
};

