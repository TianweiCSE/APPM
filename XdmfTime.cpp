#include "XdmfTime.h"



XdmfTime::XdmfTime()
{
}

XdmfTime::XdmfTime(const double time)
{
	this->time = time;
	std::stringstream ss;
	ss << "<Time Value=\"" << time << "\" />";
	this->setStartTag(ss.str());
}


XdmfTime::~XdmfTime()
{
}
