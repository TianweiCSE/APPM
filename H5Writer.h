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




