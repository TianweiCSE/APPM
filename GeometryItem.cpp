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

void GeometryItem::setTag(const int tag)
{
	this->tag = tag;
}

int GeometryItem::getIndex() const
{
	return this->index;
}

int GeometryItem::getTag() const
{
	return this->tag;
}