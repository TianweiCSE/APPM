#include "Face.h"

/** ytw:
 * - Dual face center is located where the corresponding primal vertex lies?
 * - Face normal may not form a right-handed system with the order of vertices in vertexList.
 */

Face::Face()
{
}

Face::Face(const std::vector<Vertex*> & faceVertices)
{
	const int nVertices = faceVertices.size();
	assert(nVertices >= 3);
	this->vertexList = faceVertices;
	init();
}

Face::Face(const std::vector<Edge*> & faceEdges)
{
	for (auto edge : faceEdges) {
		assert(edge != nullptr);
	}
	this->edgeList = faceEdges;
	init();
}


Face::~Face()
{
}

std::vector<Edge*> Face::getEdgeList() const
{
	return this->edgeList;
}

std::vector<Vertex*> Face::getVertexList() const
{
	assert(vertexList.size() > 0);
	return vertexList;
	//std::vector<Vertex*> faceVertices;

	//Edge * e = edgeList[0];
	//Vertex * v = e->getVertexB();

	//const int nEdges = edgeList.size();
	//for (int i = 0; i < nEdges; i++) {
	//	faceVertices.push_back(v);
	//	Vertex * otherVertex = e->getOppositeVertex(v);
	//	assert(otherVertex != nullptr);
	//	Edge * otherEdge = nullptr;
	//	for (auto edge : edgeList) {
	//		if (edge != e && edge->hasVertex(otherVertex)) {
	//			otherEdge = edge;
	//		}
	//	}
	//	assert(otherEdge != nullptr);
	//	e = otherEdge;
	//	v = otherVertex;
	//}

	//assert(faceVertices.size() == edgeList.size());
	//return faceVertices;
}

std::vector<Cell*> Face::getCellList() const
{
	return cellList;
}


bool Face::hasFaceEdges(const std::vector<Edge*> faceEdges) const
{
	std::vector<int>::iterator it;
	std::vector<int> edgeIdx;
	for (auto edge : faceEdges) {
		edgeIdx.push_back(edge->getIndex());
	}
	std::sort(edgeIdx.begin(), edgeIdx.end());
	// remove duplicate edge indices 
	it = std::unique(edgeIdx.begin(), edgeIdx.end());
	edgeIdx.resize( std::distance(edgeIdx.begin(), it) );

	if (faceEdges.size() != edgeList.size()) {
		return false;
	}

	std::vector<int> thisFaceEdgeIdx;
	for (auto edge : edgeList) {
		thisFaceEdgeIdx.push_back(edge->getIndex());
	}
	std::sort(thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end());
	
	std::vector<int> result(edgeIdx.size());
	it = std::set_intersection(edgeIdx.begin(), edgeIdx.end(), thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end(), result.begin());
	const int distance = std::distance(result.begin(), it);
	if (distance < edgeList.size()) {
		return false;
	}

	return true;
}


const int Face::getOrientation(const Edge * edge)
{
	assert(edge != nullptr);
	// check if edge is in edgeList of this face
	bool isMember = false;
	for (auto e : edgeList) {
		isMember |= (edge == e);
	}
	assert(isMember);

	// check if (a x b) is parallel or anti-parallel with face normal
	const Eigen::Vector3d posA = edge->getVertexA()->getPosition();
	const Eigen::Vector3d posB = edge->getVertexB()->getPosition();
	const Eigen::Vector3d a = (posA - center).normalized();
	const Eigen::Vector3d b = (posB - center).normalized();
	const Eigen::Vector3d n = (a.cross(b)).normalized();
	const double dotProduct = faceNormal.dot(n);
	assert(dotProduct != 0);	
	const int orientation = (dotProduct > 0) ? 1 : -1;
	return orientation;
}

bool Face::isBoundary() const
{
	if (cellList.size() == 0) { // in 2D
		return false; 
	}
	else { // in 3D
		return cellList.size() == 1;
	}

	//if (cellList.size() > 0) {
	//	return cellList.size() == 1;
	//}
	//else {
	//	for (auto edge : edgeList) {
	//		if (edge->isBoundary()) {
	//			return true;
	//		}
	//	}
	//	return false;
	//}
}

bool Face::hasBoundaryEdges() const
{
	for (auto edge : edgeList) {
		if (edge->isBoundary()) {
			return true;
		}
	}
	return false;
}

const Eigen::Vector3d Face::getCenter() const
{
	return center;
}

const double Face::getArea() const
{
	return area;
}

// ytw:
// True if other is adjacent to this, or other = this
bool Face::hasCommonCell(const Face * other) const
{
	std::vector<int> thisCellIdx;
	for (auto cell : cellList) {
		thisCellIdx.push_back(cell->getIndex());
	}
	std::vector<int> otherCellIdx;
	for (auto cell : other->cellList) {
		otherCellIdx.push_back(cell->getIndex());
	}
	std::vector<int> result(thisCellIdx.size());
	std::vector<int>::iterator it;
	it = std::set_intersection(thisCellIdx.begin(), thisCellIdx.end(), otherCellIdx.begin(), otherCellIdx.end(), result.begin());

	return std::distance(result.begin(), it) > 0;
}


Cell * Face::getCommonCell(const Face * other) const
{
	assert(hasCommonCell(other));
	std::vector<int> thisCellIdx;
	for (auto cell : cellList) {
		thisCellIdx.push_back(cell->getIndex());
	}
	std::vector<int> otherCellIdx;
	for (auto cell : other->cellList) {
		otherCellIdx.push_back(cell->getIndex());
	}
	std::vector<int> result(thisCellIdx.size());
	std::vector<int>::iterator it;
	it = std::set_intersection(thisCellIdx.begin(), thisCellIdx.end(), otherCellIdx.begin(), otherCellIdx.end(), result.begin());
	const int dist = std::distance(result.begin(), it);
	assert(dist > 0);
	assert(dist == 1);    // ytw: rule out the case when other = this
	const int cellIdx = result[0];
	Cell * cell = nullptr;
	for (auto c : cellList) {
		if (c->getIndex() == cellIdx) {
			return c;
		}
	}
	assert(cell != nullptr);  // ytw: seems to be redundent
	return cell;
}

void Face::setAdjacient(Cell * cell)
{
	assert(cell != nullptr);
	this->cellList.push_back(cell);
}

const Eigen::Vector3d Face::getNormal() const
{
	return faceNormal;
}

void Face::setNormal(const Eigen::Vector3d & fn)
{
	assert(fn.norm() > 0);
	this->faceNormal = fn;
}

void Face::init()
{
	assert((edgeList.size() > 0) ^ (vertexList.size() > 0));  // boolean XOR operator: a ^ b

	if (edgeList.size() == 0) {
		const int nVertices = vertexList.size();
		assert(nVertices >= 3);
		for (int i = 0; i < nVertices; i++) {
			Vertex * A = vertexList[i];
			Vertex * B = vertexList[(i + 1) % nVertices];
			assert(A != nullptr);
			assert(B != nullptr);
			Edge * edge = A->getAdjacientEdge(B);
			assert(edge != nullptr);
			edgeList.push_back(edge);
		}
	}

	if (vertexList.size() == 0) {
		// Determine face vertices:
		//     e0       e1      e2 
		// A ------ B ----- C ------ ...
		vertexList = std::vector<Vertex*>();
		// Choose initial vector appropriately
		//Vertex * V = edgeList[0]->getVertexA();
		Vertex * V = edgeList.front()->getCoincidentVertex(edgeList.back());
		for (auto edge : edgeList) {
			vertexList.push_back(V);
			V = edge->getOppositeVertex(V);
		}
		assert(vertexList.size() == edgeList.size());
	}

	// If this face is a triangle: 
	// - check if its vertices form a right-handed system,    ytw: ????
	// - and their normal vector is parallel to z-axis        ytw: But it is not rightly checked.
	if (vertexList.size() == 3) {
		Eigen::Vector3d center;
		center.setZero();
		for (int i = 0; i < 3; i++) {
			center += vertexList[i]->getPosition();
		}
		center *= 1./3.;

		Eigen::Matrix3d vertexPos;
		for (int i = 0; i < 3; i++) {
			vertexPos.col(i) = vertexList[i]->getPosition();
		}
		//std::cout << "Face vertex pos: " << std::endl << vertexPos.transpose() << std::endl;
		//std::cout << "Face center: " << center.transpose() << std::endl;

		/// test if triangles face normals are all parallel to z-axis.
		for (int i = 0; i < 3; i++) {
			const Eigen::Vector3d v0 = vertexList[i]->getPosition();
			const Eigen::Vector3d v1 = vertexList[(i + 1) % 3]->getPosition();

			const Eigen::Vector3d a = v0 - center;
			const Eigen::Vector3d b = v1 - v0;
			const Eigen::Vector3d n = a.normalized().cross(b.normalized());
			//std::cout << "n = " << n.transpose() << std::endl;
			const double n_dot_z = n.dot(Eigen::Vector3d::UnitZ());
			assert(n_dot_z > 0);
			//assert(n.segment(0, 2).norm() < 16 * std::numeric_limits<double>::epsilon());
		}
	}

	assert(edgeList.size() >= 3);
	assert(vertexList.size() >= 3);
	assert(edgeList.size() == vertexList.size());
	assert(isListOfVerticesUnique());
	assert(isListOfEdgesUnique());

	// Register this face being adjacient to edges
	for (auto edge : edgeList) {
		edge->setAdjacient(this);
	}

	if (edgeList.size() == 3) {
		center = getCircumCenter();
	}
	else {
		// arithmetic mean of face vertices     ytw: Assuming the face is dual, does this center coincide with primal vertex? 
		center = Eigen::Vector3d::Zero();
		std::vector<Vertex*> vertexList = getVertexList();
		for (auto v : vertexList) {
			center += v->getPosition();
		}
		center /= vertexList.size();   
	}

	// Determine face normal from face center and vector of first edge
	// ytw: The face normal may not form a right-handed system with the order of vertices in vertexList.
	const Eigen::Vector3d posA = edgeList[0]->getVertexA()->getPosition();
	const Eigen::Vector3d posB = edgeList[0]->getVertexB()->getPosition();
	const Eigen::Vector3d a = (posA - center).normalized();
	const Eigen::Vector3d b = (posB - center).normalized();
	faceNormal = (a.cross(b)).normalized();
	assert(faceNormal.norm() > 0);

	this->area = computeArea();
}

/** 
 * Check if vertices are unique. 
 */
bool Face::isListOfVerticesUnique() const
{
	assert(vertexList.size() >= 3);
	std::vector<int> vertexIdx;
	std::vector<int>::iterator it;
	for (auto v : vertexList) {
		vertexIdx.push_back(v->getIndex());
	}
	std::sort(vertexIdx.begin(), vertexIdx.end());
	it = std::unique(vertexIdx.begin(), vertexIdx.end());
	return it == vertexIdx.end();
}

/** 
 * Check if list of edges is unique.
 */
bool Face::isListOfEdgesUnique() const
{
	assert(edgeList.size() >= 3);
	std::vector<int> edgeIdx;
	std::vector<int>::iterator it;
	for (auto e : edgeList) {
		edgeIdx.push_back(e->getIndex());
	}
	std::sort(edgeIdx.begin(), edgeIdx.end());
	it = std::unique(edgeIdx.begin(), edgeIdx.end());
	return it == edgeIdx.end();
}

const Eigen::Vector3d Face::getCircumCenter() const
{
	assert(edgeList.size() == 3);
	const std::vector<Vertex*> vertexList = getVertexList();

	// Definition of circumcenter: http://mathworld.wolfram.com/Circumcircle.html
	double z = 0;
	Eigen::MatrixXd D(3, 4);
	for (int i = 0; i < 3; i++) {
		const Eigen::Vector3d pos = vertexList[i]->getPosition();
		D(i, 0) = pow(pos(0), 2) + pow(pos(1), 2);
		D(i, 1) = pos(0);
		D(i, 2) = pos(1);
		D(i, 3) = 1.0;
		z += pos(2);
	}
	z /= 3;

	Eigen::Matrix3d Bx, By;
	Bx.col(0) = D.col(0);
	Bx.col(1) = D.col(2);
	Bx.col(2) = D.col(3);

	By.col(0) = D.col(0);
	By.col(1) = D.col(1);
	By.col(2) = D.col(3);
	
	const double a = D.rightCols(3).determinant();
	assert(abs(a) > 0);
	const double bx = -1 * Bx.determinant();
	const double by =      By.determinant();
	return Eigen::Vector3d(-bx / (2*a), -by / (2*a), z);
}

const double Face::computeArea() const
{
	double area = 0;
	const Eigen::Vector3d fc = getCenter();
	for (auto edge : edgeList) {
		const Eigen::Vector3d posA = edge->getVertexA()->getPosition();
		const Eigen::Vector3d posB = edge->getVertexB()->getPosition();
		const Eigen::Vector3d a = posA - fc;
		const Eigen::Vector3d b = posB - posA;
		const double temp = 0.5 * a.cross(b).norm();
		assert(temp > 0);
		area += temp;
	}
	return area;
}

void Face::setFluidType(const Face::FluidType & fluidType)
{
	this->fluidType = fluidType;
}

const Face::FluidType Face::getFluidType() const
{
	return this->fluidType;
}

std::ostream & operator<<(std::ostream & os, const Face & obj)
{
	os << "Face " << obj.getIndex() << ": {";
	for (auto edge : obj.edgeList) {
		os << edge->getIndex() << ",";
	}
	os << "}";
	return os;
}
