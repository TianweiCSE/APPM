#pragma once
#include "XmlElement.h"
class XdmfRoot :
	public XmlElement
{
public:
	XdmfRoot();
	~XdmfRoot();

	friend std::ostream & operator<<(std::ostream & os, const XdmfRoot & obj);
};

