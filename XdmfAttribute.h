#pragma once
#include "XmlElement.h"
#include <cassert>
#include <sstream>

class XdmfAttribute :
	public XmlElement
{
public:

	enum class Type {
		Scalar, Vector
	};
	friend std::ostream & operator<<(std::ostream & os, const XdmfAttribute::Type & obj);

	enum class Center {
		Node, Cell
	};
	friend std::ostream & operator<<(std::ostream & os, const XdmfAttribute::Center & obj);
	
	struct Tags {
		std::string name;
		Type type;
		Center center;

		Tags();
		Tags(const std::string & name, const Type & type, const Center & center);
	};

	XdmfAttribute();
	XdmfAttribute(const XdmfAttribute::Tags & tags);
	~XdmfAttribute();
};

