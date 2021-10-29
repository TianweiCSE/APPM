#include "XdmfDomain.h"



XdmfDomain::XdmfDomain()
{
	this->setStartTag("<Domain>");
	this->setEndTag("</Domain>");
}


XdmfDomain::~XdmfDomain()
{
}
