#include "Numerics.h"



Numerics::Numerics()
{
}


Numerics::~Numerics()
{
}

Eigen::Vector3d Numerics::raviartThomasBasis(const int idx, const Eigen::Vector3d & x)
{
	switch (idx) {
	case 0:
		return Eigen::Vector3d(x(0), x(1), 0);
	case 1:
		return Eigen::Vector3d(x(0) - 1, x(1), 0);
	case 2:
		return Eigen::Vector3d(x(0), x(1) - 1, 0);
	case 3:
		return Eigen::Vector3d(0, 0, x(2) - 1);
	case 4:
		return Eigen::Vector3d(0, 0, x(2));
	default:	
		return Eigen::Vector3d::Zero();
	}
	return Eigen::Vector3d::Zero();
}


