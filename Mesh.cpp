#include "Mesh.h"

Mesh::Mesh()
{
	//std::cout << "Call to Mesh()" << std::endl;
}

Mesh::Mesh(const std::string & meshPrefix)
{
	//std::cout << "Call to Mesh(string)" << std::endl;
	this->meshPrefix = meshPrefix;
}

Mesh::~Mesh()
{
	//std::cout << "Call to ~Mesh()" << std::endl;
	if (vertexList.size() > 0) {
		for (auto v : vertexList) {
			delete v;
		}
		vertexList.resize(0);
	}
	if (edgeList.size() > 0) {
		for (auto e : edgeList) {
			delete e;
		}
		edgeList.resize(0);
	}
	if (faceList.size() > 0) {
		for (auto f : faceList) {
			delete f;
		}
		faceList.resize(0);
	}
	//if (cellList.size() > 0) {
	//	for (auto c : cellList) {
	//		delete c;
	//	}
	//	cellList.resize(0);
	//}
}

void Mesh::writeToFile()
{
	std::cout << std::endl;
	std::cout << "Mesh:     " << this->meshPrefix << std::endl;
	std::cout << "Vertices: " << vertexList.size() << std::endl;
	std::cout << "Edges:    " << edgeList.size() << std::endl;
	std::cout << "Faces:    " << faceList.size() << std::endl;
	std::cout << "Cells:    " << cellList.size() << std::endl;

	//std::ofstream file;

	const int nVertices = getNumberOfVertices();
	const int nEdges = getNumberOfEdges();

	// Write vertices to file
	Eigen::Matrix3Xd positions(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		positions.col(i) = vertexList[i]->getPosition();
	}
	//file = std::ofstream(this->meshPrefix + "-vertices.dat");
	//file << positions.transpose() << std::endl;

	const std::string h5filename = this->meshPrefix + "-mesh.h5";
	H5Writer h5writer(h5filename);


	/// ---------------------- write info about Vertex -------------------
	h5writer.writeDoubleMatrix(positions, "/vertexPos");

	Eigen::VectorXi vertexIdx(nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexIdx[i] = getVertex(i)->getIndex();   
	}
	h5writer.writeIntVector(vertexIdx, "/vertexIdx");   /// would end up with [0,1,2,...]

	Eigen::VectorXi isBoundaryVertex(nVertices);
	for (int i = 0; i < nVertices; i++) {
		isBoundaryVertex(i) = getVertex(i)->isBoundary();
	}
	h5writer.writeIntVector(isBoundaryVertex, "/isBoundaryVertex");

	Eigen::VectorXi vertexType = getVertexTypes();
	h5writer.writeIntVector(vertexType, "/vertexType");

	// ------------------Create incidence maps and write them to file---------------------
	createIncidenceMaps();
	Eigen::sparseMatrixToFile(edge2vertexMap, this->meshPrefix + "-e2v.dat");
	Eigen::sparseMatrixToFile(face2edgeMap, this->meshPrefix + "-f2e.dat");
	Eigen::sparseMatrixToFile(cell2faceMap, this->meshPrefix + "-c2f.dat");

	/// -------------------- write info about Edge ----------------------
	Eigen::VectorXd edgeLength(nEdges);
	for (int i = 0; i < nEdges; i++) {
		edgeLength(i) = getEdge(i)->getLength();
	}
	h5writer.writeDoubleVector(edgeLength, "/edgeLength");

	Eigen::Matrix2Xi edge2Vertex(2, nEdges);
	for (int i = 0; i < nEdges; i++) {
		edge2Vertex(0, i) = getEdge(i)->getVertexA()->getIndex();
		edge2Vertex(1, i) = getEdge(i)->getVertexB()->getIndex();
	}
	h5writer.writeIntMatrix(edge2Vertex, "/edge2vertex");

	Eigen::VectorXi edgeIdx(nEdges);
	for (int i = 0; i < nEdges; i++) {
		edgeIdx(i) = getEdge(i)->getIndex();
	}
	h5writer.writeIntVector(edgeIdx, "/edgeIdx"); /// would end up with [0,1,2,...]

	Eigen::VectorXi edgeType(nEdges);
	for (int i = 0; i < nEdges; i++) {
		const Edge * edge = getEdge(i);
		bool isBoundaryA = edge->getVertexA()->isBoundary();
		bool isBoundaryB = edge->getVertexB()->isBoundary();
		const int nBoundaryVertices = isBoundaryA + isBoundaryB;
		edgeType(i) = nBoundaryVertices;
	}
	h5writer.writeIntVector(edgeType, "/edgeType");

	/// -------------------- write info about Face --------------------
	const int nFaces = getNumberOfFaces();
	Eigen::MatrixXi f2v(3, nFaces);
	for (int j = 0; j < nFaces; j++) {
		const Face * face = faceList[j];
		std::vector<Vertex*> faceVertices = face->getVertexList();
		const int nFaceVertices = faceVertices.size();
		if (nFaceVertices > f2v.rows()) {
			const int nRowsOld = f2v.rows();
			f2v.conservativeResize(nFaceVertices, f2v.cols());
			f2v.bottomRows(nFaceVertices - nRowsOld).array() = -1;
		}
		for (int k = 0; k < nFaceVertices; k++) {
			f2v(k, j) = faceVertices[k]->getIndex();
		}
	}
	const int f2v_maxCoeff = f2v.array().abs().maxCoeff();

	std::ofstream file(this->meshPrefix + "-f2v.dat");
	file << f2v.transpose() << std::endl;

	Eigen::Matrix3Xd fc(3, nFaces);
	Eigen::Matrix3Xd fn(3, nFaces);
	Eigen::VectorXi faceIdx(nFaces);
	Eigen::VectorXi faceBoundary(nFaces);
	Eigen::VectorXd faceArea(nFaces);
	for (int i = 0; i < nFaces; i++) {
		fc.col(i) = getFace(i)->getCenter();
		fn.col(i) = getFace(i)->getNormal();
		faceIdx(i) = getFace(i)->getIndex();
		faceBoundary(i) = getFace(i)->isBoundary();
		faceArea(i) = getFace(i)->getArea();
	}
	h5writer.writeDoubleMatrix(fc, "/faceCenter");
	h5writer.writeIntVector(faceIdx, "/faceIndex");  /// would end up with [0,1,2,...]
	h5writer.writeIntVector(faceBoundary, "/isFaceBoundary");
	h5writer.writeDoubleVector(fn, "/faceNormal");
	assert((faceArea.array() > 0).all());
	h5writer.writeDoubleVector(faceArea, "/faceArea");

	const std::vector<int> face2vertexIdx = getXdmfTopology_face2vertexIndices();
	h5writer.writeStdVector(face2vertexIdx, "/face2vertex");


	/// ---------------------- write info about Cell ------------------------
	const int nCells = getNumberOfCells();
	Eigen::Matrix3Xd cellCenters(3, nCells);
	for (int i = 0; i < nCells; i++) {
		cellCenters.col(i) = getCell(i)->getCenter();
	}
	h5writer.writeDoubleMatrix(cellCenters, "/cellCenter");

	Eigen::VectorXi cellIdx(nCells);
	for (int i = 0; i < nCells; i++) {
		cellIdx(i) = getCell(i)->getIndex();
	}
	h5writer.writeIntVector(cellIdx, "/cellIndex");  /// would end up with [0,1,2,...]

	const std::vector<int> c2vIdx = getXdmfTopology_cell2vertexIndices();
	h5writer.writeStdVector(c2vIdx, "/cell2vertex");

	Eigen::VectorXd cellVolume(nCells);
	for (int i = 0; i < nCells; i++) {
		cellVolume(i) = getCell(i)->getVolume();
	}
	h5writer.writeDoubleVector(cellVolume, "/cellVolume");

}

void Mesh::writeXdmf()
{
	XdmfRoot root;
	XdmfDomain domain;

	XdmfGrid treeGrid(XdmfGrid::Tags("Grid of Grids", XdmfGrid::GridType::Tree));

	XdmfGrid vertexGrid = getXdmfVertexGrid();
	treeGrid.addChild(vertexGrid);

	XdmfGrid edgeGrid = getXdmfEdgeGrid();
	treeGrid.addChild(edgeGrid);

	XdmfGrid surfaceGrid = getXdmfSurfaceGrid();
	treeGrid.addChild(surfaceGrid);

	domain.addChild(treeGrid);
	root.addChild(domain);

	std::string filename = this->meshPrefix + "-mesh.xdmf";
	std::ofstream file(filename);
	file << root << std::endl;

	writeXdmfVolumeMesh();
}

Vertex * Mesh::addVertex(const Eigen::Vector3d & position)
{
	const int index = vertexList.size();
	Vertex * vertex = new Vertex(position, index);
	vertexList.push_back(vertex);
	return vertex;
}

Edge * Mesh::addEdge(Vertex * A, Vertex * B)
{
	assert(A != nullptr);
	assert(B != nullptr);
	if (!isConnected(A, B)) {
		Edge * edge = new Edge(A, B);
		edge->setIndex(edgeList.size());
		edgeList.push_back(edge);
		return edge;
	}
	else {
		return getEdge(A, B);
	}
}

Face * Mesh::addFace(const std::vector<Edge*> & faceEdges)
{
	assert(faceEdges.size() >= 3);
	if (!isConnected(faceEdges)) {
		Face * face = new Face(faceEdges);
		face->setIndex(faceList.size());
		faceList.push_back(face);
		return face;
	}
	else {
		Face * face = getFace(faceEdges);
		assert(face != nullptr);
		return face;
	}
	return nullptr;
}

Face * Mesh::addFace(const std::vector<Vertex*> & faceVertices)
{
	assert(faceVertices.size() >= 3);
	std::vector<Edge*> faceEdges;
	const int nVertices = faceVertices.size();
	for (int i = 0; i < nVertices; i++) {
		Vertex * A = faceVertices[i];
		Vertex * B = faceVertices[(i + 1) % nVertices];
		Edge * edge = getEdge(A, B);
		assert(edge != nullptr);
		faceEdges.push_back(edge);
	}
	if (!isConnected(faceEdges)) {
		Face * face = new Face(faceEdges);
		face->setIndex(faceList.size());
		faceList.push_back(face);
		return face;
	}
	else {
		Face * face = getFace(faceEdges);
		assert(face != nullptr);
		return face;
	}
	return nullptr;
}

TriPrism * Mesh::addTriPrism(const std::vector<Face*>& sideFaces, Face * bottomFace, Face * topFace)
{
	for (auto f : sideFaces) {
		assert(f != nullptr);
	}
	assert(bottomFace != nullptr);
	assert(topFace != nullptr);

	std::vector<Face*> cellFaces;
	for (auto f : sideFaces) {
		cellFaces.push_back(f);
	}
	cellFaces.push_back(bottomFace);
	cellFaces.push_back(topFace);

	if (!isConnected(cellFaces)) {
		const int idx = cellList.size();
		TriPrism * prism = new TriPrism(cellFaces);
		prism->setIndex(idx);
		cellList.push_back(prism);
		return prism;
	}
	else {
		Cell * cell = getCell(cellFaces);
		assert(cell != nullptr);
		return static_cast<TriPrism*>(cell);
	}
	return nullptr;
}

Cell * Mesh::addCell(const std::vector<Face*>& cellFaces)
{
	assert(cellFaces.size() >= 4);
	if (!isConnected(cellFaces)) {
		const int idx = cellList.size();
		Cell * cell = new Cell(cellFaces);
		cell->setIndex(idx);
		cellList.push_back(cell);
		return cell;
	}
	else {
		Cell * cell = getCell(cellFaces);
		assert(cell != nullptr);
		return cell;
	}
	return nullptr;
}

Vertex * Mesh::getVertex(const int index) const
{
	assert(index >= 0);
	assert(index < vertexList.size());
	Vertex * vertex = vertexList[index];
	assert(vertex != nullptr);
	assert(vertex->getIndex() == index);
	return vertex;
}

Edge * Mesh::getEdge(const int index) const
{
	assert(index >= 0);
	assert(index < edgeList.size());
	return edgeList[index];
}

Face * Mesh::getFace(const int index) const
{
	assert(index >= 0);
	assert(index < faceList.size());
	return faceList[index];
}

Face * Mesh::getFace(const std::vector<Edge*>& faceEdges)
{
	assert(this->isConnected(faceEdges));
	assert(faceEdges[0]->hasConnectedFace(faceEdges));
	return faceEdges[0]->getConnectedFace(faceEdges);
}

Cell * Mesh::getCell(const int index) const
{
	assert(index >= 0);
	assert(index < cellList.size());
	return cellList[index];
}

Cell * Mesh::getCell(const std::vector<Face*>& cellFaces)
{
	assert(cellFaces.size() >= 2);
	const Face * face0 = cellFaces[0];
	const Face * face1 = cellFaces[1];
	assert(face0 != face1);
	assert(face0->hasCommonCell(face1));
	Cell * cell = face0->getCommonCell(face1);
	for (auto face : cellFaces) {
		assert(cell->hasFace(face));
	}
	return cell;
}

const std::string Mesh::getPrefix() const
{
	return meshPrefix;
}

void Mesh::createIncidenceMaps()
{
	create_vertexCoordinates();   /// useless
	create_edge2vertex_map();
	create_face2edge_map();
	create_cell2face_map();
}

const int Mesh::getNumberOfVertices() const
{
	return vertexList.size();
}

const int Mesh::getNumberOfEdges() const
{
	return edgeList.size();
}

const int Mesh::getNumberOfFaces() const
{
	return faceList.size();
}

const int Mesh::getNumberOfCells() const
{
	return cellList.size();
}

const std::vector<int> Mesh::getXdmfTopology_cell2vertexIndices() const
{
	std::vector<int> data;
	const int nCells = cellList.size();
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = cellList[i];
		data.push_back(16); // polyhedron type
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const int nCellFaces = cellFaces.size();
		data.push_back(nCellFaces); // number of faces
		for (int j = 0; j < nCellFaces; j++) {
			const Face * face = cellFaces[j];
			const std::vector<Vertex*> faceVertices = face->getVertexList();
			const int nFaceVertices = faceVertices.size();
			data.push_back(nFaceVertices); // number of vertices
			for (auto vertex : faceVertices) {
				data.push_back(vertex->getIndex());
			}
		}
	}
	return data;
}

const std::vector<int> Mesh::getXdmfTopology_face2vertexIndices() const
{
	std::vector<int> f2v;
	const int nFaces = getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		const Face * face = faceList[i];
		std::vector<Vertex*> faceVertices = face->getVertexList();
		const int nFaceVertices = faceVertices.size();

		int faceType = 0;
		switch (nFaceVertices) {
		case 3: faceType = 4; break;
		case 4: faceType = 5; break;
		default: faceType = 3; break;
		}

		f2v.push_back(faceType);
		if (faceType == 3) { // polygon face
			f2v.push_back(nFaceVertices);
		}
		for (int j = 0; j < faceVertices.size(); j++) {
			f2v.push_back(faceVertices[j]->getIndex());
		}
	}
	return f2v;
}

const std::vector<Cell*> Mesh::getCells() const
{
	return cellList;
}

const std::vector<Face*> Mesh::getFaces() const
{
	return faceList;
}

void Mesh::check() const
{
	//for (auto cell : cellList) {
	//	std::cout << "Cell idx " << cell->getIndex() << std::endl;
	//	const std::vector<Face*> cellFaces = cell->getFaceList();
	//	for (auto face : cellFaces) {
	//		std::cout << "Face " << face->getIndex() << "; ";
	//		std::cout << "orientation = " << cell->getOrientation(face);
	//		std::cout << std::endl;
	//	}
	//}

	if (cellList.size() > 0) {
		// face is either a boundary face, or it has two adjacient cells; orientation +/-1.
		for (auto face : faceList) {
			std::vector<Cell*> faceCells = face->getCellList();
			if (face->isBoundary()) {
				assert(faceCells.size() == 1);
			}
			else {
				assert(faceCells.size() == 2);
				const Cell * cell0 = faceCells[0];
				const Cell * cell1 = faceCells[1];
				std::vector<int> orientations = { cell0->getOrientation(face), cell1->getOrientation(face) };
				int sumOrientations = 0;
				for (auto value : orientations) {
					sumOrientations += value;
				}
				if (sumOrientations != 0) {
					std::cout << "Face idx " << face->getIndex() << " has wrong set of orientation: "  << orientations[0] << ", " << orientations[1] << std::endl;
					std::cout << "face normal:  " << face->getNormal().transpose() << std::endl;
					std::cout << "face center:  " << face->getCenter().transpose() << std::endl;
					std::cout << "cell0 center: " << cell0->getCenter().transpose() << std::endl;
					std::cout << "cell1 center: " << cell1->getCenter().transpose() << std::endl;
					std::cout << "(c0 - fc):    " << (cell0->getCenter() - face->getCenter()).transpose() << std::endl;
					std::cout << "(c1 - fc):    " << (cell1->getCenter() - face->getCenter()).transpose() << std::endl;
					std::cout << "(c0 - fc).fn:    " << (cell0->getCenter() - face->getCenter()).dot(face->getNormal()) << std::endl;
					std::cout << "(c1 - fc).fn:    " << (cell1->getCenter() - face->getCenter()).dot(face->getNormal()) << std::endl;

				}

				assert(orientations[0] + orientations[1] == 0);
			}
		}
	}

}

const Eigen::SparseMatrix<int>& Mesh::get_f2eMap() const
{
	return face2edgeMap;
}

const Eigen::SparseMatrix<int>& Mesh::get_e2vMap() const
{
	return edge2vertexMap;
}

std::vector<Edge*> Mesh::makeContinuousLoop(std::vector<Edge*> edges)
{
	std::vector<Edge*> loop;
	std::vector<Edge*>::iterator it;
	Edge * edge = edges[0];
	Vertex * V = edge->getVertexB();
	const int nEdges = edges.size();

	for (int i = 0; i < nEdges; i++) {
		// add edge to new list
		loop.push_back(edge);
		// remove edge from old list
		for (it = edges.begin(); it != edges.end(); it++) {
			if (*it == edge) { break; }
		}
		edges.erase(it);
		// update vertex that is searched 
		V = edge->getOppositeVertex(V);
		// find edge that has this vertex
		for (it = edges.begin(); it != edges.end(); it++) {
			edge = *it;
			if (edge->hasVertex(V)) {
				break;
			}
		}
	}
	assert(edges.size() == 0);
	assert(loop.size() == nEdges);
	return loop;
}

bool Mesh::isConnected(const Vertex * A, const Vertex * B) const
{
	return (A->isAdjacientTo(B) && B->isAdjacientTo(A));
}

bool Mesh::isConnected(const std::vector<Edge*> & faceEdges) const
{
	bool isAlreadyDefinedFace = true;
	for (auto edge : faceEdges) {
		isAlreadyDefinedFace &= edge->hasConnectedFace(faceEdges);
		if (!isAlreadyDefinedFace) {
			break;
		}
	}
	return isAlreadyDefinedFace;
}

bool Mesh::isConnected(const std::vector<Face*> & cellFaces) const
{
	bool isAlreadyDefined = true;
	for (auto face : cellFaces) {
		isAlreadyDefined &= face->hasCommonCell(face);
	}
	return isAlreadyDefined;
}

Edge * Mesh::getEdge(Vertex * A, Vertex * B)
{
	assert(A != nullptr);
	assert(B != nullptr);
	if (isConnected(A, B)) {
		return A->getAdjacientEdge(B);
	}
	else {
		return nullptr;
	}
}

void Mesh::create_vertexCoordinates()
{
	const int nVertices = vertexList.size();
	vertexCoordinates = Eigen::Matrix3Xd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}
}

void Mesh::create_edge2vertex_map()
{
	const int nEdges = edgeList.size();
	const int nVertices = vertexList.size();
	this->edge2vertexMap = Eigen::SparseMatrix<int>(nEdges, nVertices);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nEdges; i++) {
		auto e = edgeList[i];
		int idxE = e->getIndex();
		/// i = idxE
		const Vertex * A = e->getVertexA();
		const Vertex * B = e->getVertexB();
		const int idxA = A->getIndex();
		const int idxB = B->getIndex();
		triplets.push_back(T(idxE, idxA, -1));
		triplets.push_back(T(idxE, idxB, +1));
	}
	edge2vertexMap.setFromTriplets(triplets.begin(), triplets.end());
	edge2vertexMap.makeCompressed();
}

void Mesh::create_face2edge_map()
{
	const int nFaces = faceList.size();
	const int nEdges = edgeList.size();
	this->face2edgeMap = Eigen::SparseMatrix<int>(nFaces, nEdges);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nFaces; i++) {
		auto face = faceList[i];
		const int idxF = face->getIndex();
		/// i = idxF
		const std::vector<Edge*> faceEdges = face->getEdgeList();
		for (auto edge : faceEdges) {
			auto idxE = edge->getIndex();
			const int orientation = face->getOrientation(edge);
			triplets.push_back(T(idxF, idxE, orientation));
		}
	}
	face2edgeMap.setFromTriplets(triplets.begin(), triplets.end());
	face2edgeMap.makeCompressed();
}

void Mesh::create_cell2face_map()
{
	const int nFaces = faceList.size();
	const int nCells = cellList.size();
	this->cell2faceMap = Eigen::SparseMatrix<int>(nCells, nFaces);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nCells; i++) {
		auto cell = cellList[i];
		const int idxC = cell->getIndex();
		/// i = idxC
		const std::vector<Face*> cellFaces = cell->getFaceList();
		for (auto face : cellFaces) {
			auto idxF = face->getIndex();
			const int orientation = cell->getOrientation(face);
			triplets.push_back(T(idxC, idxF, orientation));
		}
	}
	cell2faceMap.setFromTriplets(triplets.begin(), triplets.end());
	cell2faceMap.makeCompressed();

	// Check cell-to-face map
	// Each face has either one adjacient cell (boundary face), or 
	// it has two adjacient cells (non-boundary face = inner face). 
	// Moreover, because of positive/negative orientations, 
	// a colwise sum of the matrix elements yields +/-1 (boundary) or 0 (non-boundary).
	//Eigen::VectorXi colwiseSum = cell2faceMap.transpose() * Eigen::VectorXi::Ones(nCells);
	//Eigen::VectorXi isBoundaryFace(nFaces);
	//for (int i = 0; i < nFaces; i++) {
	//	isBoundaryFace(i) = getFace(i)->isBoundary();
	//}
	//assert(colwiseSum.size() == nFaces);
	//Eigen::VectorXi notBoundaryFace = isBoundaryFace.array() - 1;
	//assert((colwiseSum.cwiseProduct(notBoundaryFace).array() == 0).all());
	//Eigen::VectorXi temp = colwiseSum.cwiseAbs().cwiseProduct(isBoundaryFace);
	//assert((temp.array() >= 0).all());
	//assert((temp.array() <= 1).all());
	//assert(temp.sum() == isBoundaryFace.sum());
}

XdmfGrid Mesh::getXdmfVertexGrid() const
{
	XdmfGrid vertexGrid(XdmfGrid::Tags("Vertex Grid"));

	// Topology and topology dataItem
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexIdx";
		XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Polyvertex, getNumberOfVertices(), 1));
		topology.addChild(XdmfDataItem(
			XdmfDataItem::Tags(
				{ getNumberOfVertices()}, 
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
			ss.str()
		));
		vertexGrid.addChild(topology);
	}

	// Geometry and geometry data item
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexPos";
		XdmfGeometry geometry;
		geometry.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices(), 3 },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF),
				ss.str()
			));
		vertexGrid.addChild(geometry);
	}

	// Attribute: vertex index
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexIdx";
		XdmfAttribute attribute(XdmfAttribute::Tags("Vertex Index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		vertexGrid.addChild(attribute);
	}

	// Attribute: vertex type
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexType";
		XdmfAttribute attribute(XdmfAttribute::Tags("Vertex Type", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		vertexGrid.addChild(attribute);
	}
	return vertexGrid;
}

XdmfGrid Mesh::getXdmfEdgeGrid() const
{
	XdmfGrid edgeGrid(XdmfGrid::Tags("Edge grid"));

	// Topology
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/edge2vertex";
		XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Polyline, getNumberOfEdges(), 2));
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ 2*getNumberOfEdges() },   /// why 2 times
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		edgeGrid.addChild(topology);
	}

	// Geometry
	{
		XdmfGeometry geometry;
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexPos";
		geometry.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices(), 3 },   /// why number of vertices
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF), 
				ss.str()
			));
		edgeGrid.addChild(geometry);
	}

	// Attribute: edge index
	XdmfAttribute edgeIdxAttribute(
		XdmfAttribute::Tags("Edge index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	std::stringstream ss;
	ss << this->meshPrefix << "-mesh.h5:/edgeIdx";
	edgeIdxAttribute.addChild(
		XdmfDataItem(
			XdmfDataItem::Tags(
				{ getNumberOfEdges() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
			ss.str()
		)
	);
	edgeGrid.addChild(edgeIdxAttribute);

	// Attribute: edge type
	{
		XdmfAttribute edgeTypeAttribute(
			XdmfAttribute::Tags("Edge Type", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/edgeType";
		edgeTypeAttribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfEdges() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		edgeGrid.addChild(edgeTypeAttribute);
	}


	return edgeGrid;
}

XdmfGrid Mesh::getXdmfSurfaceGrid() const
{
	H5Reader h5reader(this->meshPrefix + "-mesh.h5");

	XdmfGrid surfaceGrid(XdmfGrid::Tags("Face Grid"));

	// Topology
	XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, getNumberOfFaces()));  /// why not polygon

	// Topology DataItem
	{
		const int nElements = h5reader.readDataSize("/face2vertex");
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/face2vertex";
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ nElements }, 
					XdmfDataItem::NumberType::Int, 
					XdmfDataItem::Format::HDF), 
				ss.str()               /// Why the data item is not the reture value of getXdmfTopology_face2vertexIndices()?
			)
		);
		surfaceGrid.addChild(topology);
	}

	// Geometry and Geometry DataItem
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexPos";
		XdmfGeometry geometry;
		geometry.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices(), 3 }, 
					XdmfDataItem::NumberType::Float, 
					XdmfDataItem::Format::HDF), 
				ss.str()
			)
		);
		
		surfaceGrid.addChild(geometry);
	}

	// Attribute: face index

	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/faceIndex";
		XdmfAttribute attributeFaceIdx(XdmfAttribute::Tags("Face Index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::string faceIdxStringAttribute = ss.str();
		attributeFaceIdx.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfFaces() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				faceIdxStringAttribute));
		surfaceGrid.addChild(attributeFaceIdx);
	}

	// Attribute: face area
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/faceArea";
		XdmfAttribute attributeFaceArea(XdmfAttribute::Tags("Face Area", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::string faceAreaStringAttribute = ss.str();
		attributeFaceArea.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfFaces() },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF),
				faceAreaStringAttribute));
		surfaceGrid.addChild(attributeFaceArea);
	}
	
	return surfaceGrid;
}

XdmfGrid Mesh::getXdmfVolumeGrid() const
{
	H5Reader h5reader(this->meshPrefix + "-mesh.h5");
	XdmfGrid volumeGrid(XdmfGrid::Tags("VolumeGrid"));

	// Topology
	XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, getNumberOfCells()));

	// Topology DataItem
	{
		const int nElements = h5reader.readDataSize("/cell2vertex");
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/cell2vertex";
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ nElements },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
					ss.str()
			));
		volumeGrid.addChild(topology);
	}

	// Geometry and Geometry DataItem
	XdmfGeometry geometry;
	std::stringstream ss;
	ss << this->meshPrefix << "-mesh.h5:/vertexPos";
	geometry.addChild(
		XdmfDataItem(
			XdmfDataItem::Tags(
				{ getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
		));

	volumeGrid.addChild(geometry);

	// Attribute: cell index
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Cell Index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/cellIndex";
		std::string bodyString = ss.str();
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfCells() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF),
				bodyString));
		volumeGrid.addChild(attribute);
	}

	// Attribute: cell volume
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Cell Volume", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
        ss << this->meshPrefix << "-mesh.h5:/cellVolume";
		std::string bodyString = ss.str();
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfCells() },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF),
				bodyString));
		volumeGrid.addChild(attribute);
	}

	return volumeGrid;
}

void Mesh::writeXdmfVolumeMesh() const
{
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid volumeGrid = getXdmfVolumeGrid();
	domain.addChild(volumeGrid);
	root.addChild(domain);
	std::ofstream file(this->meshPrefix + "-volume.xdmf");
	file << root << std::endl;
}

const Eigen::VectorXi Mesh::getVertexTypes() const
{
	const int nV = getNumberOfVertices();
	Eigen::VectorXi types(nV);
	for (int i = 0; i < nV; i++) {
		types(i) = static_cast<int>(getVertex(i)->getType());
	}
	return types;
}

const Eigen::VectorXi Mesh::getEdgeTypes() const
{
	const int nE = getNumberOfEdges();
	Eigen::VectorXi types(nE);
	for (int i = 0; i < nE; i++) {
		types(i) = static_cast<int>(getEdge(i)->getType());
	}
	return types;
}

const Eigen::VectorXi Mesh::getFaceTypes() const
{
	const int nF = getNumberOfFaces();
	Eigen::VectorXi types(nF);
	for (int i = 0; i < nF; i++) {
		types(i) = getFace(i)->isBoundary();
	}
	return types;
}

Mesh::MeshInfo Mesh::getMeshInfo() const
{
	MeshInfo meshInfo;
	assert(getNumberOfVertices() > 0);
	assert(getNumberOfEdges() > 0);
	assert(getNumberOfFaces() > 0);
	assert(getNumberOfCells() > 0);

	// Vertices
	const Eigen::VectorXi vertexTypes = getVertexTypes();
	const int boundaryVertexType = static_cast<int>(Vertex::Type::Boundary);
	const int terminalVertexType = static_cast<int>(Vertex::Type::Terminal);
	meshInfo.nVertices = getNumberOfVertices();
	meshInfo.nVerticesTerminal = (vertexTypes.array() == terminalVertexType).count();
	meshInfo.nVerticesBoundary = (vertexTypes.array() == boundaryVertexType).count() + meshInfo.nVerticesTerminal;
	if (meshPrefix.substr(0, 6) == "primal") {
		assert(meshInfo.nVerticesTerminal > 0);
		assert(meshInfo.nVerticesBoundary > 0);
	}

	// Edges
	const Eigen::VectorXi edgeTypes = getEdgeTypes();
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);
	const int interiorToBoundaryEdgeType = static_cast<int>(Edge::Type::InteriorToBoundary);
	const int boundaryEdgeType = static_cast<int>(Edge::Type::Boundary);
	meshInfo.nEdges = getNumberOfEdges();
	meshInfo.nEdgesInner =
		(edgeTypes.array() == interiorEdgeType).count() +
		(edgeTypes.array() == interiorToBoundaryEdgeType).count();

	// Faces
	const Eigen::VectorXi faceTypes = getFaceTypes();
	const int isBoundaryFace = 1;
	meshInfo.nFaces = getNumberOfFaces();
	meshInfo.nFacesInner = (faceTypes.array() != isBoundaryFace).count();

	// Cells
	meshInfo.nCells = getNumberOfCells();

	return meshInfo;
}


//void Mesh::writeXdmf_volume()
//{
//	std::stringstream ss;
//	std::stringstream body;
//	XmlElement * dataItem;
//
//	std::string filename = this->meshPrefix + "-volume.xdmf";
//	std::cout << "Write XDMF file: " << filename << std::endl;
//
//	XmlElement root("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">", "</Xdmf>");
//	XmlElement * domain = new XmlElement("<Domain>", "</Domain>");
//	root.addChild(domain);
//
//	XmlElement * grid = new XmlElement("<Grid Name=\"Appm cells\">", "</Grid>");
//	domain->addChild(grid);
//
//	ss << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" << getNumberOfCells() << "\">";
//	XmlElement * topology = new XmlElement(ss.str(), "</Topology>");
//	grid->addChild(topology);
//	
//	ss = std::stringstream();
//
//	// Version with HDF
//	const std::string h5filename = (this->meshPrefix) + "-mesh.h5";
//	H5Reader h5reader(h5filename);
//	const std::string dataname = "/cell2vertex";
//	const int nElements = h5reader.readDataSize(dataname); 
//	body = std::stringstream();
//	body << (this->meshPrefix + "-mesh.h5:/cell2vertex");
//	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Format=\"HDF\">";
//
//	dataItem = new XmlElement(ss.str(), "</DataItem>", body.str());
//	topology->addChild(dataItem);
//
//	XmlElement * geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
//	grid->addChild(geometry);
//	ss = std::stringstream();
//	ss << "<DataItem Dimensions=\"" << vertexCoordinates.cols() << " " << 3 << "\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
//	body = std::stringstream();
//	//body << vertexCoordinates.transpose();
//	body << (this->meshPrefix + "-mesh.h5:/vertexPos");
//	dataItem = new XmlElement(ss.str(), "</DataItem>", body.str());
//	geometry->addChild(dataItem);
//
//	// Attribute: cell index
//	ss = std::stringstream();
//	ss << "<Attribute Name=\"Cell index\" AttributeType=\"Scalar\" Center=\"Cell\">";
//	XmlElement * attribute = new XmlElement(ss.str(), "</Attribute>");
//	body = std::stringstream();
//	body << (this->meshPrefix + "-mesh.h5:/cellIndex");
//	ss = std::stringstream();
//	ss << "<DataItem Format=\"HDF\" DataType=\"Int\" Precision=\"4\" Dimensions=\"" << getNumberOfCells() << "\">";
//	dataItem = new XmlElement(ss.str(), "</DataItem>", body.str());
//	attribute->addChild(dataItem);
//	grid->addChild(attribute);
//
//
//	std::ofstream file(filename);
//	file << "<?xml version=\"1.0\" ?>" << std::endl;
//	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
//	file << root << std::endl;
//}
