#pragma once
#include "Cell.h"
class TriPrism :
	public Cell
{
public:
	TriPrism();
	TriPrism(const std::vector<Face*> & faces);
	~TriPrism();
};

