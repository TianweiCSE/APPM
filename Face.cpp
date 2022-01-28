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
	vertexList = faceVertices;
	vertexListExtended = faceVertices;
	init();
}

Face::Face(const std::vector<Edge*> & faceEdges)
{
	for (auto edge : faceEdges) {
		assert(edge != nullptr);
	}
	edgeList = faceEdges;
	init();
}


Face::~Face()
{
}

const std::vector<Edge*> Face::getEdgeList() const
{
	return this->edgeList;
}

const std::vector<Vertex*> Face::getVertexList() const
{
	assert(vertexList.size() > 0);
	return vertexList;
}

const std::vector<Vertex*> Face::getVertexListExtended() const
{
	assert(vertexListExtended.size() > 0);
	return vertexListExtended;
}

const std::vector<Cell*> Face::getCellList() const
{
	return cellList;
}

const std::vector<Face*> Face::getSubFaceList() const {
	return subFaceList;
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


int Face::getOrientation(const Edge * edge) const
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
	const Eigen::Vector3d posB = edge->getVertexMid()!= nullptr ? edge->getVertexMid()->getPosition()
															    : edge->getVertexB()->getPosition();
	const Eigen::Vector3d a = (posA - center).normalized();
	const Eigen::Vector3d b = (posB - center).normalized();
	const Eigen::Vector3d n = (a.cross(b)).normalized();
	const double dotProduct = faceNormal.dot(n);
	// assert(dotProduct != 0);	
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

double Face::getArea() const
{
	return area;
}

double Face::getProjectedArea() const {
	if (isPlane()) {
		return area;
	}
	else {
		double area_ = 0;
		const Eigen::Vector3d fc = getCenter();
		const int nVertices = vertexList.size();
		for (int i = 0; i < vertexList.size(); i++) {
			const Eigen::Vector3d posA = vertexList[i]->getPosition();
			const Eigen::Vector3d posB = vertexList[(i+1) % nVertices]->getPosition();
			const Eigen::Vector3d a = posA - fc;
			const Eigen::Vector3d b = posB - posA;
			const double temp = 0.5 * a.cross(b).norm();
			assert(temp > 0);
			area_ += temp;
		}
		return area_;
	}
}

// ytw:
// True if other is adjacent to this, or other = this
bool Face::hasCommonCell(Face * other) const
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


Cell * Face::getCommonCell(Face * other) const
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

void Face::reverseNormal()
{
	faceNormal *= -1;
}

const Eigen::Vector3d Face::computeNormal() {

	if (isPlane()) {
		const Eigen::Vector3d posA = vertexListExtended[0]->getPosition();
		const Eigen::Vector3d posB = vertexListExtended[1]->getPosition();
		const Eigen::Vector3d a = (posA - getCenter()).normalized();
		const Eigen::Vector3d b = (posB - getCenter()).normalized();
		faceNormal = (a.cross(b)).normalized();
	}
	else { // non plane
		const Eigen::Vector3d posA = vertexList[0]->getPosition();
		const Eigen::Vector3d posB = vertexList[1]->getPosition();
		const Eigen::Vector3d a = (posA - getCenter()).normalized();
		const Eigen::Vector3d b = (posB - getCenter()).normalized();
		faceNormal = (a.cross(b)).normalized();
	}
	assert(faceNormal.norm() > 0);
	return faceNormal;
}

const Eigen::Vector3d Face::computeCenter() {
	if (isPlane()) {
		if (vertexListExtended.size() == 3) {  // If triangle, take circumcenter.
			center = getCircumCenter();
		}
		else { // If polygon (and could possibly contain non-straight edges), take the average of ALL vertices.
			center.setZero();
			for (auto v : vertexListExtended) {
				center += v->getPosition();
			}
			center /= vertexListExtended.size(); 
		}
	}
	else { // If not plane, take the average of topological vertices. (Namely the center lies on its projection)
		center.setZero();
		for (auto v : vertexList) {
			center += v->getPosition();
		}
		center /= vertexList.size();
	}
	return center;
}

double Face::computeArea()
{	
	area = 0;
	if (isPlane()) {
		const Eigen::Vector3d fc = getCenter();
		const int nVertices = vertexListExtended.size();
		for (int i = 0; i < vertexListExtended.size(); i++) {
			const Eigen::Vector3d posA = vertexListExtended[i]->getPosition();
			const Eigen::Vector3d posB = vertexListExtended[(i+1) % nVertices]->getPosition();
			const Eigen::Vector3d a = posA - fc;
			const Eigen::Vector3d b = posB - posA;
			const double temp = 0.5 * a.cross(b).norm();
			assert(temp > 0);
			area += temp;
		}
	}
	else {
		assert(subFaceList.size() != 0);
		for (const Face* f : subFaceList) {
			area += f->getArea();
		}
	}
	return area;
}

void Face::init()
{
	assert((edgeList.size() > 0) ^ (vertexList.size() > 0));  // boolean XOR operator: a ^ b

	// If the face is constructed from ordered vertexList
	if (edgeList.size() == 0) {
		// std::cout << "building from vertexList" << std::endl;
		const int nVertices = vertexList.size();
		assert(nVertices >= 3);
		for (int i = 0; i < nVertices; i++) {
			Vertex * A = vertexList[i];
			Vertex * B = vertexList[(i + 1) % nVertices];
			assert(A != nullptr);
			assert(B != nullptr);
			Edge * edge = A->getAdjacientEdge(B);
			if (edge == nullptr) { // If the edge is not built in advance, e.g. when it is a sub-face
				edge = new Edge(A, B);  
			}
			// assert(edge != nullptr); 					 
			edgeList.push_back(edge);
		}
	}

	// If the face is constructed from ordered edgeList
	if (vertexList.size() == 0) {
		// Determine face vertices:
		//     e0       e1      e2 
		// A ------ B ----- C ------ ...
		// Choose initial vector appropriately
		// std::cout << "building from edgeList" << std::endl;
		assert(isContinuousLoop(edgeList));
		Vertex * V = edgeList.front()->getCoincidentVertex(edgeList.back());
		assert(V != nullptr);  // check that the input edgeList is well ordered.
		for (auto edge : edgeList) {
			vertexList.push_back(V);
			vertexListExtended.push_back(V);
			if (edge->getVertexMid() != nullptr) {
				vertexListExtended.push_back(edge->getVertexMid());
			}
			V = edge->getOppositeVertex(V);
		}
		assert(vertexList.size() == edgeList.size());  // The number of topological vertices should match that of edges.
		assert(vertexListExtended.size() >= edgeList.size());
	}

	// Determine face normal from face center and vector of first edge
	// ytw: The face normal may not form a right-handed system with the order of vertices in vertexList.
	computeCenter();
	computeNormal();
	addSubFaces();
	computeArea();
	// If this face is a triangle: 
	// - check if its vertices form a right-handed system,    ytw: ????
	// - and their normal vector is parallel to z-axis        ytw: But it is not rightly checked.
	/*
	if (vertexListExtended.size() == 3) {
		Eigen::Vector3d center;
		center.setZero();
		for (int i = 0; i < 3; i++) {
			center += vertexListExtended[i]->getPosition();
		}
		center *= 1./3.;

		Eigen::Matrix3d vertexPos;
		for (int i = 0; i < 3; i++) {
			vertexPos.col(i) = vertexListExtended[i]->getPosition();
		}
		//std::cout << "Face vertex pos: " << std::endl << vertexPos.transpose() << std::endl;
		//std::cout << "Face center: " << center.transpose() << std::endl;

		/// test if triangles face normal is parallel to z-axis.
		const double n_dot_z = getNormal().dot(Eigen::Vector3d::UnitZ());
		if (n_dot_z < std::numeric_limits<double>::epsilon()) {   // If not right handed, switch the first and second vertices
			Vertex* temp = vertexList[0];
			vertexList[0] = vertexList[1];
			vertexList[1] = temp;
			vertexListExtended = vertexList;
		}
		updateNormal();
	}*/
	
	assert(area > 0);
	assert(edgeList.size() >= 3);
	assert(vertexList.size() >= 3);
	assert(vertexListExtended.size() >= 3);
	assert(edgeList.size() == vertexList.size());
	assert(isListOfVerticesUnique());
	assert(isListOfEdgesUnique());

	// Register this face being adjacient to edges
	for (auto edge : edgeList) {
		edge->setAdjacient(this);
	}	
}

bool Face::isPlane() const {
	int nNonStraightEdges = 0;
	for (const Edge* e : edgeList) {
		if (e->getVertexMid() != nullptr) nNonStraightEdges++;
	}
	return nNonStraightEdges <= 1;
}

void Face::addSubFaces() {
	if (!isPlane()) {
		// std::cout << "adding subFaces" << std::endl;
		assert(edgeList.size() == 4); // In our case, non-plane faces always have four (straight or non-straight) edges.
		assert(vertexList.size() == 4);
		if (vertexList.size() == vertexListExtended.size() - 2) {
			Vertex *f1v1, *f1v2, *f1v3, *f1v4, *f2v1, *f2v2, *f2v3, *f2v4;
			if (edgeList[0]->getVertexMid() == nullptr) {
				f1v1 = edgeList[0]->getCoincidentVertex(edgeList[1]);
				f1v2 = edgeList[1]->getVertexMid();
				f1v3 = edgeList[3]->getVertexMid();
				f1v4 = edgeList[3]->getCoincidentVertex(edgeList[0]);
				f2v1 = f1v2;
				f2v2 = edgeList[1]->getOppositeVertex(f1v1);
				f2v3 = edgeList[2]->getOppositeVertex(f2v2);
				f2v4 = f1v3;
			}
			else {
				f1v1 = edgeList[0]->getCoincidentVertex(edgeList[1]);
				f1v2 = edgeList[1]->getOppositeVertex(f1v1);
				f1v3 = edgeList[2]->getVertexMid();
				f1v4 = edgeList[0]->getVertexMid();
				f2v1 = f1v4;
				f2v2 = f1v3;
				f2v3 = edgeList[2]->getOppositeVertex(f1v2);
				f2v4 = edgeList[3]->getOppositeVertex(f2v3);
			}
			assert(f1v1 != nullptr && f1v2 != nullptr && f1v3 != nullptr && f1v4 != nullptr &&
				   f2v1 != nullptr && f2v2 != nullptr && f2v3 != nullptr && f2v4 != nullptr);
			// Note that the vertices at the same location need rebuilding, otherwise the original vertices 
			// will have adjacent edges that are not belonging to the mesh. 
			Face* f1 = new Face({f1v1->copy(), f1v2->copy(), f1v3->copy(), f1v4->copy()});
			Face* f2 = new Face({f2v1->copy(), f2v2->copy(), f2v3->copy(), f2v4->copy()});
			assert(std::abs(f1->getNormal()[2]) < 1e-12);
			assert(std::abs(f2->getNormal()[2]) < 1e-12);
			if (f1->getNormal().dot(this->getNormal()) < 0) f1->reverseNormal();
			if (f2->getNormal().dot(this->getNormal()) < 0) f2->reverseNormal();
			subFaceList.push_back(f1);
			subFaceList.push_back(f2);
		}
		else if (vertexList.size() == vertexListExtended.size() - 3) {
			int tmp_idx = 0;
			while (edgeList[tmp_idx]->getVertexMid() != nullptr) tmp_idx++;
			assert(edgeList[tmp_idx]->getVertexMid() == nullptr);
			Vertex* v1 = edgeList[tmp_idx]->getCoincidentVertex(edgeList[(tmp_idx + 1) % 4]);
			Vertex* v2 = edgeList[(tmp_idx + 1) % 4]->getVertexMid();
			Vertex* v3 = edgeList[(tmp_idx + 1) % 4]->getOppositeVertex(v1);
			Vertex* v4 = edgeList[(tmp_idx + 2) % 4]->getVertexMid();
			Vertex* v5 = edgeList[(tmp_idx + 2) % 4]->getOppositeVertex(v3);
			Vertex* v6 = edgeList[(tmp_idx + 3) % 4]->getVertexMid();
			Vertex* v7 = edgeList[(tmp_idx + 3) % 4]->getOppositeVertex(v5);
			Vertex* v8 = new Vertex(v4->getPosition() + v6->getPosition() - v5->getPosition());
			
			assert(v1 != nullptr && v2 != nullptr && v3 != nullptr && v4 != nullptr && v5 != nullptr && v6 != nullptr && v7 != nullptr && v8 != nullptr);
			// Note that the vertices at the same location need rebuilding, otherwise the original vertices 
			// will have adjacent edges that are not belonging to the mesh. 
			Face* f1 = new Face({v1->copy(), v2->copy(), v8->copy(), v6->copy(), v7->copy()});
			Face* f2 = new Face({v2->copy(), v3->copy(), v4->copy(), v8->copy()});
			Face* f3 = new Face({v8->copy(), v4->copy(), v5->copy(), v6->copy()});
			assert(f1->getNormal().segment(0,2).norm() < 1e-12); // f1 is supposed to be parallel to z-axis
			assert(std::abs(f2->getNormal()[2]) < 1e-12);
			assert(std::abs(f3->getNormal()[2]) < 1e-12);
			if (f1->getNormal().dot(this->getNormal()) < 0) f1->reverseNormal();
			if (f2->getNormal().dot(this->getNormal()) < 0) f2->reverseNormal();
			if (f3->getNormal().dot(this->getNormal()) < 0) f3->reverseNormal();
			subFaceList.push_back(f1);
			subFaceList.push_back(f2);
			subFaceList.push_back(f3);
		}
		else assert(false);
	}
}

/** 
 * Check if vertices are unique. 
 */
bool Face::isListOfVerticesUnique() const
{	
	/*
	assert(vertexList.size() >= 3);
	std::vector<int> vertexIdx;
	std::vector<int>::iterator it;
	for (auto v : vertexList) {
		vertexIdx.push_back(v->getIndex());
	}
	std::sort(vertexIdx.begin(), vertexIdx.end());
	it = std::unique(vertexIdx.begin(), vertexIdx.end());
	return it == vertexIdx.end();
	*/
	// The reason why checking indices is not working is that for sub-faces the indices are all set to -1.
	// Thus, we adopt a more straightforward way, namely checking the uniqueness of pointers
	for (int i = 0; i < vertexList.size(); i++) {
		for (int j = i + 1; j < vertexList.size(); j++) {
			if (vertexList[i] == vertexList[j]) {
				return false;
			}
		}
	}
	return true;
}

/** 
 * Check if list of edges is unique.
 */
bool Face::isListOfEdgesUnique() const
{	
	/*
	assert(edgeList.size() >= 3);
	std::vector<int> edgeIdx;
	std::vector<int>::iterator it;
	for (auto e : edgeList) {
		edgeIdx.push_back(e->getIndex());
	}
	std::sort(edgeIdx.begin(), edgeIdx.end());
	it = std::unique(edgeIdx.begin(), edgeIdx.end());
	return it == edgeIdx.end();*/
	for (int i = 0; i < edgeList.size(); i++) {
		for (int j = i + 1; j < edgeList.size(); j++) {
			if (edgeList[i] == edgeList[j]) {
				return false;
			}
		}
	}
	return true;

}

const Eigen::Vector3d Face::getCircumCenter() const
{	
	//std::cout << " --------------computing circum center ----------------" << std::endl;
	assert(vertexListExtended.size() == 3);  // This must be a triangle geometrically.
	const std::vector<Vertex*> vertexList = getVertexListExtended();

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
	//std::cout << "Face idx: " << getIndex() << std::endl;
	//std::cout << vertexListExtended[0]->getIndex() << " " << vertexListExtended[0]->getPosition() << std::endl;
	//std::cout << vertexListExtended[1]->getIndex() << " " << vertexListExtended[1]->getPosition() << std::endl;
	//std::cout << vertexListExtended[2]->getIndex() << " " << vertexListExtended[2]->getPosition() << std::endl;

	Eigen::Matrix3d Bx, By;
	Bx.col(0) = D.col(0);
	Bx.col(1) = D.col(2);
	Bx.col(2) = D.col(3);

	By.col(0) = D.col(0);
	By.col(1) = D.col(1);
	By.col(2) = D.col(3);
	
	const double a = D.rightCols(3).determinant();
	assert(std::abs(a) > 0);
	const double bx = -1 * Bx.determinant();
	const double by =      By.determinant();
	return Eigen::Vector3d(-bx / (2*a), -by / (2*a), z);
}

void Face::setFluidType(const Face::FluidType & fluidType)
{
	this->fluidType = fluidType;
}

const Face::FluidType Face::getFluidType() const
{
	return this->fluidType;
}

bool Face::isContinuousLoop(std::vector<Edge*> edges) const {
	bool flag = true;
	for (int i = 0; i < edges.size(); i++) {
		flag = edges[i]->hasCoincidentVertex(edges[(i+1) % edges.size()]);
	}
	return flag;
}

std::ostream & operator<<(std::ostream & os, const Face & obj)
{
	os << "Face " << obj.getIndex() << ": {";
	for (auto edge : obj.getEdgeList()) {
		os << edge->getIndex() << ",";
	}
	os << "}";
	return os;
}
