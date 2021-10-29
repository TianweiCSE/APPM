#include "GeometryItem.h"



GeometryItem::GeometryItem()
{
}

GeometryItem::GeometryItem(const int index)
{
	this->index = index;
}


GeometryItem::~GeometryItem()
{
}

void GeometryItem::setIndex(const int index)
{
	assert(index >= 0);
	this->index = index;
}

const int GeometryItem::getIndex() const
{
	return this->index;
}
