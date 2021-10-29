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

	const int readDataSize(const std::string & dataname) const;

private:
	std::string filename;
	H5File file;
};

