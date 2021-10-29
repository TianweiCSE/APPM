#include "XdmfAttribute.h"
#include "XdmfDataItem.h"



XdmfAttribute::XdmfAttribute()
{
	this->setStartTag("<Attribute>");
	this->setEndTag("</Attribute>");
}

XdmfAttribute::XdmfAttribute(const XdmfAttribute::Tags & tags)
{
	std::stringstream ss;
	ss << "<Attribute";
	ss << " Name=\"" << tags.name << "\"";
	ss << " AttributeType=\"" << tags.type << "\"";
	ss << " Center=\"" << tags.center << "\"";
	ss << ">";
	this->setStartTag(ss.str());
	this->setEndTag("</Attribute>");
}


XdmfAttribute::~XdmfAttribute()
{
}

std::ostream & operator<<(std::ostream & os, const XdmfAttribute::Type & obj)
{
	switch (obj) {
	case XdmfAttribute::Type::Scalar:
		os << "Scalar"; break;
	case XdmfAttribute::Type::Vector:
		os << "Vector"; break;
	default:
		os << "XdmfAttribute::Type not implemented"; assert(false); break;
	}
	return os;
}

std::ostream & operator<<(std::ostream & os, const XdmfAttribute::Center & obj)
{
	switch (obj) {
	case XdmfAttribute::Center::Node:
		os << "Node"; break;
	case XdmfAttribute::Center::Cell:
		os << "Cell"; break;
	default:
		os << "XdmfAttribute::Center not implemented"; assert(false); break;
	}
	return os;
}

XdmfAttribute::Tags::Tags()
{
	this->name = "";
	this->type = XdmfAttribute::Type::Scalar;
	this->center = XdmfAttribute::Center::Node;
}

XdmfAttribute::Tags::Tags(const std::string & name, const Type & type, const Center & center)
{
	this->name = name;
	this->type = type;
	this->center = center;
}
