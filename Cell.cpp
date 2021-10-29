#include "Cell.h"



Cell::Cell()
{
}

Cell::Cell(const std::vector<Face*>& faces)
{
	//std::cout << "Create cell from faces: {";
	//for (auto face : faces) {
	//	std::cout << face->getIndex() << ",";
	//}
	//std::cout << "}" << std::endl;
	//std::cout << "Face vertices: " << std::endl;
	//for (auto face : faces) {
	//	std::cout << "face idx " << face->getIndex() << ": ";
	//	const std::vector<Vertex*> f2v = face->getVertexList();
	//	for (auto v : f2v) {
	//		std::cout << v->getIndex() << ",";
	//	}
	//	std::cout << std::endl;
	//}

	this->faceList = faces;
	for (auto face : faceList) {
		face->setAdjacient(this);
	}

	// determine cell center
	std::vector<Face*> zFaces;
	for (auto face : faceList) {
		const Eigen::Vector3d fn = face->getNormal();
		if (fn.cross(Eigen::Vector3d(0, 0, 1)).norm() < 100 * std::numeric_limits<double>::epsilon()) {
			zFaces.push_back(face);
		}
	}
	const int n_zFaces = zFaces.size();
	// ytw: this block is for degugging.
	if (n_zFaces < 2) {
		for (auto face : faceList) {
			const Eigen::Vector3d fn = face->getNormal();
			const double fn_x_zUnit_norm = fn.cross(Eigen::Vector3d(0, 0, 1)).norm();
			std::cout << "Face " << face->getIndex() << ": " << "(fn x zUnit).norm() = " << fn_x_zUnit_norm << std::endl;
			if (fn.cross(Eigen::Vector3d(0, 0, 1)).norm() < 100 * std::numeric_limits<double>::epsilon()) {
				zFaces.push_back(face);
			}
		}
	}
	assert(zFaces.size() == 2);
	//std::cout << "z-Faces: ";
	//for (auto f : zFaces) {
	//	std::cout << f->getIndex() << ",";
	//}
	//std::cout << std::endl;
	this->center = Eigen::Vector3d::Zero();
	center = 1. / 2. * (zFaces[0]->getCenter() + zFaces[1]->getCenter());

	// Set cell volume
	volume = 0;
	for (auto face : faceList) {
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		double h = std::abs(fn.dot(center - face->getCenter()));
		double dV = 1. / 3. * fA * h;
		assert(dV > 0);
		volume += dV;
	}
	assert(volume > 0);

	//Eigen::Matrix3d M;
	//M.col(0) = zFaces[0]->getCenter();
	//M.col(1) = zFaces[1]->getCenter();
	//M.col(2) = this->center;
	//std::cout << "face/cell center: " << std::endl << M << std::endl;
}


Cell::~Cell()
{
}

const double Cell::getVolume() const
{
	return volume;
}

const std::vector<Face*> Cell::getFaceList() const
{
	return faceList;
}

bool Cell::hasFace(const Face * face) const
{
	assert(face != nullptr);
	for (auto f : faceList) {
		if (f == face) {
			return true;
		}
	}
	return false;
}

const int Cell::getOrientation(const Face * face) const
{
	assert(face != nullptr);
	bool isMember = false;
	for (auto f : faceList) {
		if (f == face) {
			isMember |= true;
			break;
		}
	}
	assert(isMember);
	const Eigen::Vector3d a = face->getCenter() - this->center;
	const Eigen::Vector3d fn = face->getNormal();
	const int orientation = (a.dot(fn) > 0) ? 1 : -1;
	return orientation;
}

const Eigen::Vector3d & Cell::getCenter() const
{
	return center;
}

const Eigen::Matrix3Xd Cell::getVertexCoordinates() const
{
	std::vector<Vertex*> vertexList;
	for (auto face : faceList) {
		for (auto vertex : face->getVertexList()) {
			vertexList.push_back(vertex);
		}
	}
	std::sort(vertexList.begin(), vertexList.end());    // ytw: How does it work? vertexList contains Vertex*!
	std::vector<Vertex*>::iterator it;
	it = std::unique(vertexList.begin(), vertexList.end());
	int nVertices = std::distance(vertexList.begin(), it);
	assert(nVertices >= 4);

	Eigen::Matrix3Xd coords(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		const Vertex * vertex = vertexList[i];
		coords.col(i) = vertex->getPosition();
	}
	return coords;
}

void Cell::setFluidType(const FluidType & type)
{
	this->fluidType = type;
}

const Cell::FluidType Cell::getFluidType() const
{
	return this->fluidType;
}
