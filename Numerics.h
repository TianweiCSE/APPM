#pragma once

#include <Eigen/Dense>

class Numerics
{
public:
	Numerics();
	~Numerics();

	/** 
	* Raviart-Thomas basis functions for a triangular prism.
	*/
	static Eigen::Vector3d raviartThomasBasis(const int idx, const Eigen::Vector3d & x);
};

