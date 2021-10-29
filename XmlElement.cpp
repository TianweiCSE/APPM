#include "XmlElement.h"



XmlElement::XmlElement()
{
}

XmlElement::~XmlElement()
{
}

void XmlElement::addChild(const XmlElement & child)
{
	children.push_back(child);
}

void XmlElement::setStartTag(const std::string & tag)
{
	this->startTag = tag;
}

void XmlElement::setEndTag(const std::string & tag)
{
	this->endTag = tag;
}

void XmlElement::setBody(const std::string & text)
{
	this->body = text;
}


const std::string & XmlElement::getStartTag() const
{
	return startTag;
}

const std::string & XmlElement::getEndTag() const
{
	return endTag;
}

const std::string & XmlElement::getBody() const
{
	return body;
}

const std::vector<XmlElement> & XmlElement::getChildren() const
{
	return children;
}


std::ostream & operator<<(std::ostream & os, const XmlElement & obj)
{
	os << obj.getStartTag();
	if (obj.getBody().size() > 0) {
		os << std::endl;
		os << obj.getBody();
	}
	for (auto child : obj.children) {
		os << std::endl;
		os << child;
	}
	if (obj.getEndTag().size() > 0) {
		os << std::endl;
		os << obj.endTag;
	}
	return os;
}
