#pragma once
#include "XmlElement.h"
#include <cassert>
#include <sstream>
#include <vector>

class XdmfDataItem :
	public XmlElement
{
public:
	enum class NumberType {
		Float, Int
	};
	friend std::ostream & operator<<(std::ostream & os, const XdmfDataItem::NumberType & obj);

	enum class Format {
		XML, HDF
	};

	friend std::ostream & operator<<(std::ostream & os, const XdmfDataItem::Format & obj);

	struct Tags {
		std::vector<int> dimensions = {0};
		NumberType numberType;
		Format format;

		Tags();
		Tags(const std::vector<int> dims, const NumberType & type, const Format & format);
	};

public:
	XdmfDataItem();
	XdmfDataItem(const XdmfDataItem::Tags & tag);
	XdmfDataItem(const XdmfDataItem::Tags & tag, const std::string & body);
	~XdmfDataItem();

private:
	int getPrecision(const XdmfDataItem::Tags & tag);
};

