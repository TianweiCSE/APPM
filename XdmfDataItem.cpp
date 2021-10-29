#include "XdmfDataItem.h"



XdmfDataItem::XdmfDataItem()
{
	this->setStartTag("<DataItem>");
	this->setEndTag("</DataItem>");
}

XdmfDataItem::XdmfDataItem(const XdmfDataItem::Tags & tag)
{
	std::stringstream ss;
	ss << "<DataItem";
	if (tag.dimensions.size() > 0) {
		ss << " Dimensions=\"";
		for (int i = 0; i < tag.dimensions.size(); i++) {
			ss << tag.dimensions[i];
			if (i < (tag.dimensions.size()-1)) {
				ss << " ";
			}
		}
		ss << "\"";
	}
	ss << " DataType=\"" << tag.numberType << "\"";
	ss << " Precision=\"" << getPrecision(tag);
	ss << "\"";
	ss << " Format=\"" << tag.format << "\"";
	ss << ">";
	this->setStartTag(ss.str());
	this->setEndTag("</DataItem>");
}

XdmfDataItem::XdmfDataItem(const XdmfDataItem::Tags & tag, const std::string & body)
	: XdmfDataItem(tag)
{
	this->setBody(body);
}


XdmfDataItem::~XdmfDataItem()
{
}

int XdmfDataItem::getPrecision(const XdmfDataItem::Tags & tag)
{
	switch (tag.numberType) {
	case NumberType::Int:
		return 4; 
	case NumberType::Float:
		return 8; break;
	default: assert(false);
	}
	return 0;
}

std::ostream & operator<<(std::ostream & os, const XdmfDataItem::NumberType & obj)
{
	switch (obj) {
	case XdmfDataItem::NumberType::Float:
		os << "Float"; break;
	case XdmfDataItem::NumberType::Int:
		os << "Int"; break;
	default:
		os << "DataItem::NumberType not defined"; assert(false); break;
	}
	return os;
}

std::ostream & operator<<(std::ostream & os, const XdmfDataItem::Format & obj)
{
	switch (obj) {
	case XdmfDataItem::Format::XML:
		os << "XML"; break;
	case XdmfDataItem::Format::HDF:
		os << "HDF"; break;
	default:
		os << "DataItem::Format not defined"; assert(false); break;
	}
	return os;
}

XdmfDataItem::Tags::Tags() 
{
	this->dimensions = { 0 };
	this->numberType = XdmfDataItem::NumberType::Float;
	this->format = XdmfDataItem::Format::XML;
}


XdmfDataItem::Tags::Tags(const std::vector<int> dims, const XdmfDataItem::NumberType & type, const XdmfDataItem::Format & format)
{
	this->dimensions = dims;
	this->numberType = type;
	this->format = format;
}

