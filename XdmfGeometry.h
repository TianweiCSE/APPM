#pragma once
#include "XmlElement.h"
#include <cassert>
#include <sstream>

class XdmfGeometry :
	public XmlElement
{
	enum class GeometryType {
		XYZ
	};

	friend std::ostream & operator<<(std::ostream & os, const XdmfGeometry::GeometryType & obj);

public:
	XdmfGeometry();
	XdmfGeometry(const XdmfGeometry::GeometryType & type);
	~XdmfGeometry();
};

