#pragma once
#include "XmlElement.h"
#include <sstream>

class XdmfTime :
	public XmlElement
{
public:
	XdmfTime();
	XdmfTime(const double time);
	~XdmfTime();
private:
	double time = 0.0;
};

