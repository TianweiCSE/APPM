#pragma once
#include "XmlElement.h"
#include <sstream>
#include <cassert>

class XdmfGrid :
	public XmlElement
{

public:
	enum class CollectionType {
		Spatial, Temporal
	};

	enum class GridType {
		Uniform, Tree, Collection
	};

	struct Tags {
		std::string name;
		XdmfGrid::GridType gridType;
		XdmfGrid::CollectionType collectionType;

		Tags();
		Tags(const std::string & name, const GridType & gridType = GridType::Uniform, const CollectionType & collectionType = CollectionType::Spatial);
	};

	friend std::ostream & operator<<(std::ostream & os, const XdmfGrid::CollectionType & obj);
	friend std::ostream & operator<<(std::ostream & os, const XdmfGrid::GridType & obj);

public:
	XdmfGrid();
	XdmfGrid(const XdmfGrid::Tags & gridTags);
	~XdmfGrid();

	void setTags(const XdmfGrid::Tags & gridTags);

};

