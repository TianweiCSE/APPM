#pragma once
#include "XmlElement.h"
#include <cassert>
#include <sstream>

class XdmfTopology :
	public XmlElement
{
public:
	enum class TopologyType {
		Mixed, Polyvertex, Polyline, Polygon
	};
	friend std::ostream & operator<<(std::ostream & os, const XdmfTopology::TopologyType & obj);

	struct Tags {
		TopologyType type;
		int numberOfElements = 0;
		int nodesPerElement = 0;

		Tags();
		Tags(const TopologyType & type, const int nElems, const int nodesPerElem = 0);

	};


public:
	XdmfTopology();
	XdmfTopology(const XdmfTopology::Tags & tags);
	~XdmfTopology();
};

