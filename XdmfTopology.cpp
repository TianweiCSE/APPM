#include "XdmfTopology.h"



XdmfTopology::XdmfTopology()
{
	this->setStartTag("<Topology>");
	this->setEndTag("</Topology>");
}

XdmfTopology::XdmfTopology(const XdmfTopology::Tags & tags)
{
	std::stringstream ss;
	ss << "<Topology";
	ss << " TopologyType=\"" << tags.type << "\"";
	if (tags.numberOfElements > 0) {
		ss << " NumberOfElements=\"" << tags.numberOfElements << "\"";
	}
	if (tags.type == XdmfTopology::TopologyType::Polyvertex 
		|| tags.type == XdmfTopology::TopologyType::Polyline 
		|| tags.type == XdmfTopology::TopologyType::Polygon) {
		assert(tags.nodesPerElement > 0);
	}
	if (tags.nodesPerElement > 0) {
		ss << " NodesPerElement=\"" << tags.nodesPerElement << "\"";
	}
	ss << ">";
	this->setStartTag(ss.str());
	this->setEndTag("</Topology>");
}


XdmfTopology::~XdmfTopology()
{
}

std::ostream & operator<<(std::ostream & os, const XdmfTopology::TopologyType & obj)
{
	switch (obj) {
	case XdmfTopology::TopologyType::Mixed:
		os << "Mixed"; break;
	case XdmfTopology::TopologyType::Polyvertex:
		os << "Polyvertex"; break;
	case XdmfTopology::TopologyType::Polyline:
		os << "Polyline"; break;
	case XdmfTopology::TopologyType::Polygon:
		os << "Polygon"; break;
	default:
		os << "TopologyType not defined"; assert(false); break;
	}
	return os;
}


XdmfTopology::Tags::Tags()
{
	this->numberOfElements = 0;
	this->nodesPerElement = 0;
}

XdmfTopology::Tags::Tags(const XdmfTopology::TopologyType & type, const int nElems, const int nodesPerElem)
{
	this->type = type;
	this->numberOfElements = nElems;
	this->nodesPerElement = nodesPerElem;
}