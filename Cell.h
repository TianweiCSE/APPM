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
		Undefined, Fluid, Solid
	};

	Cell();
	Cell(const std::vector<Face*> & faces);
	~Cell();

	const std::vector<Face*> getFaceList() const;
	const std::vector<Edge*> getEdgeList() const;

	bool hasFace(Face * face) const;

	// Return true if the face normal has positive inner product with the outer normal of the cell
	const int getOrientation(const Face * face) const;

	const double          computeVolume();
	const Eigen::Vector3d computeCenter();
	const double          getVolume() const;
	const Eigen::Vector3d getCenter() const;

	const Eigen::MatrixXd getVertexCoordinates() const;

	void setFluidType(const FluidType & type);
	const FluidType getFluidType() const;

private:
	double volume = 0;
	Eigen::Vector3d center;
	std::vector<Face*> faceList;
	mutable std::vector<Edge*> edgeList;

	FluidType fluidType = FluidType::Undefined;
};

