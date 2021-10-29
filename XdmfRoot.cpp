#include "XdmfRoot.h"



XdmfRoot::XdmfRoot()
{
	setStartTag(std::string("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"));
	setEndTag("</Xdmf>");
}


XdmfRoot::~XdmfRoot()
{
}

std::ostream & operator<<(std::ostream & os, const XdmfRoot & obj) {
	os << "<?xml version = \"1.0\" ?>" << std::endl;
	os << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;

	os << obj.getStartTag() << std::endl;
	if (obj.getBody().size() > 0) {
		os << obj.getBody() << std::endl;
	}
	if (obj.getChildren().size() > 0) {
		for (auto child : obj.getChildren()) {
			os << child << std::endl;
		}
	}
	os << obj.getEndTag();
	return os;
}