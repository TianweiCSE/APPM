#pragma once
#include "GeometryItem.h"

#include "Face.h"
#include <limits>
class Face;

class Cell :
	public GeometryItem
{
public:

	enum class FluidType {
		DEFAULT, FLUID, SOLID
	};

	Cell();
	Cell(const std::vector<Face*> & faces);
	~Cell();

	const double getVolume() const;
	const std::vector<Face*> getFaceList() const;
	bool hasFace(const Face * face) const;

	/// return True if the face normal coincide with the outer normal of the cell
	const int getOrientation(const Face * face) const;

	const Eigen::Vector3d & getCenter() const;
	const Eigen::Matrix3Xd getVertexCoordinates() const;

	void setFluidType(const FluidType & type);
	const FluidType getFluidType() const;

private:
	double volume = 0;
	Eigen::Vector3d center;
	std::vector<Face*> faceList;

	FluidType fluidType = FluidType::DEFAULT;
};

