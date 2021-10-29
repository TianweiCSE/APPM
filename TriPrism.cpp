#include "TriPrism.h"



TriPrism::TriPrism()
{
}

TriPrism::TriPrism(const std::vector<Face*> & faces)
	: Cell(faces)
{
	// TODO: 
	// Check that the frist (n-2) faces are quadrilaterals, and the last two are triangles
	// Check that the edge sequence in triangles are consistent with quadrilaterals

	// Create Piola map
}


TriPrism::~TriPrism()
{
}
