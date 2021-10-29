#include "XdmfGeometry.h"



XdmfGeometry::XdmfGeometry() 
	: XdmfGeometry(GeometryType())
{
}

XdmfGeometry::XdmfGeometry(const XdmfGeometry::GeometryType & type)
{
	std::stringstream ss;
	ss << "<Geometry";
	ss << " GeometryType=\"" << type << "\"";
	ss << ">";
	this->setStartTag(ss.str());
	this->setEndTag("</Geometry>");
}


XdmfGeometry::~XdmfGeometry()
{
}

std::ostream & operator<<(std::ostream & os, const XdmfGeometry::GeometryType & obj)
{
	switch (obj) {
	case XdmfGeometry::GeometryType::XYZ:
		os << "XYZ"; break;
	default:
		os << "GeometryType not defined"; assert(false); break;
	}
	return os;
}
