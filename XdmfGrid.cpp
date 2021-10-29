#include "XdmfGrid.h"



XdmfGrid::XdmfGrid()
{
	this->setStartTag("<Grid>");
	this->setEndTag("</Grid>");
}

XdmfGrid::XdmfGrid(const XdmfGrid::Tags & gridTags)
{
	setTags(gridTags);
}


XdmfGrid::~XdmfGrid()
{
}


void XdmfGrid::setTags(const XdmfGrid::Tags & gridTags)
{
	std::stringstream ss;
	ss << "<Grid";
	if (gridTags.name.size() > 0) {
		ss << " Name=\"" << gridTags.name << "\"";
	}
	if (gridTags.gridType != XdmfGrid::GridType::Uniform) {
		ss << " GridType=\"" << gridTags.gridType << "\"";
	
		if (gridTags.gridType == XdmfGrid::GridType::Collection) {
			ss << " CollectionType=\"" << gridTags.collectionType << "\"";
		}
	}
	ss << ">";
	this->setStartTag(ss.str());
	this->setEndTag("</Grid>");
}


std::ostream & operator<<(std::ostream & os, const XdmfGrid::CollectionType & obj) 
{
	switch (obj) {
	case XdmfGrid::CollectionType::Spatial:
		os << "Spatial"; break;
	case XdmfGrid::CollectionType::Temporal:
		os << "Temporal"; break;
	default:
		os << "CollectionType not defined"; assert(false); break;
	}
	return os;
}

std::ostream & operator<<(std::ostream & os, const XdmfGrid::GridType & obj)
{
	switch (obj) {
	case XdmfGrid::GridType::Uniform:
		os << "Uniform"; break;
	case XdmfGrid::GridType::Collection:
		os << "Collection"; break;
	case XdmfGrid::GridType::Tree:
		os << "Tree"; break;
	default:
		os << "GridType not defined"; assert(false); break;
	}
	return os;
}

XdmfGrid::Tags::Tags()
{
	this->name = "";
	this->gridType = GridType::Uniform;
	this->collectionType = CollectionType::Spatial;
}

XdmfGrid::Tags::Tags(const std::string & name, const GridType & gridType, const CollectionType & collectionType)
{
	this->name = name;
	this->gridType = gridType;
	this->collectionType = collectionType;
}
