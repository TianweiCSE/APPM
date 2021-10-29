#pragma once

#include <string>
#include <vector>
#include <ostream>

class XmlElement
{
public:
	XmlElement();
	~XmlElement();

	void setStartTag(const std::string & tag);
	void setEndTag(const std::string & tag);
	void setBody(const std::string & text);

	const std::string & getStartTag() const;
	const std::string & getEndTag() const;
	const std::string & getBody() const;

	void addChild(const XmlElement & child);
	const std::vector<XmlElement> & getChildren() const;

	friend std::ostream & operator<<(std::ostream & os, const XmlElement & obj);

private:
	std::string startTag;
	std::string endTag;
	std::string body;
	std::vector<XmlElement> children;
};

