#pragma once

#include <cassert>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>

class GeometryItem
{
public:
	GeometryItem();
	GeometryItem(const int index);
	~GeometryItem();

	void setIndex(const int index);
	void setTag(const int tag);

	const int getIndex() const;
	const int getTag() const;

private:
	int index = -1;
	int tag   = -1;
};

