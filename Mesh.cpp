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
	const int nVertices = getNumberOfVertices();
	const int nEdges    = getNumberOfEdges();
	const int nFaces    = getNumberOfFaces();
	const int nCells    = getNumberOfCells();

	const std::string h5filename = meshPrefix + "-mesh.h5";
	H5Writer h5writer(h5filename);

	// ------------------Create incidence maps and write them to file---------------------
	update_vertexCoordinates();
	createIncidenceMaps();
	// Eigen::sparseMatrixToFile(edge2vertexMap, this->meshPrefix + "-e2v.dat");
	// Eigen::sparseMatrixToFile(face2edgeMap,   this->meshPrefix + "-f2e.dat");
	// Eigen::sparseMatrixToFile(cell2faceMap,   this->meshPrefix + "-c2f.dat");


	/// ---------------------- write info about Vertex -------------------
	h5writer.writeDoubleMatrix(vertexCoordinates,  "/vertexPos" );  // coordinates without auxiliary vertices
	h5writer.writeDoubleMatrix(getVertexCoordinatesExtended(),  "/vertexPosExtended" );  // coordinates including auxiliary vertices
	h5writer.writeIntVector   (getVertexIndices(), "/vertexIdx" );   // would end up with [0,1,2,...]
	h5writer.writeIntVector   (getVertexTypes(),   "/vertexType");

	/// ------------------------ write info about Edge ----------------------
	h5writer.writeStdVector(getXdmfTopology_edge2vertexIndices(), "/edge2vertex"); 
	h5writer.writeIntVector(getEdgeIndices(), "/edgeIdx" ); // would end up with [0,1,2,...]
	h5writer.writeIntVector(getEdgeTypes(),   "/edgeType");
	h5writer.writeDoubleVector(getEdgeLengths(), "/edgeLength");

	/// ----------------------- write info about Face --------------------
	Eigen::MatrixXd fc(nFaces, 3);
	Eigen::MatrixXd fn(nFaces, 3);
	for (int i = 0; i < nFaces; i++) {
		fc.row(i)   = getFace(i)->getCenter();
		fn.row(i)   = getFace(i)->getNormal();
	}
	h5writer.writeStdVector(getXdmfTopology_face2vertexIndices(), "/face2vertex");
	h5writer.writeIntVector(getFaceIndices(), "/faceIdx"   );  /// would end up with [0,1,2,...]
	h5writer.writeIntVector(getFaceTypes(),   "/faceType"  );
	h5writer.writeDoubleMatrix(fc,            "/faceCenter");
	h5writer.writeDoubleMatrix(fn,            "/faceNormal");
	h5writer.writeDoubleVector(getFaceAreas(),"/faceArea"  );
	
	/// ---------------------- write info about Cell ------------------------
	Eigen::MatrixXd cellCenters(nCells, 3);
	for (int i = 0; i < nCells; i++) {
		cellCenters.row(i) = getCell(i)->getCenter();
	}
	h5writer.writeStdVector(getXdmfTopology_cell2vertexIndices(), "/cell2vertex");
	h5writer.writeIntVector(getCellIndices(),    "/cellIdx"   );  /// would end up with [0,1,2,...]
	h5writer.writeIntVector(getCellTypes(),      "/cellType"  );
	h5writer.writeDoubleMatrix(cellCenters,      "/cellCenter");
	h5writer.writeDoubleVector(getCellVolumes(), "/cellVolume");
}

void Mesh::writeXdmf()
{	
	// write the grids for vertex, edge and face
	{
		XdmfRoot root;
		XdmfDomain domain;

		XdmfGrid treeGrid(XdmfGrid::Tags("Grid of Grids", XdmfGrid::GridType::Tree));

		XdmfGrid vertexGrid = getXdmfVertexGrid();
		treeGrid.addChild(vertexGrid);

		XdmfGrid edgeGrid = getXdmfEdgeGrid();
		treeGrid.addChild(edgeGrid);

		XdmfGrid faceGrid = getXdmfFaceGrid();
		treeGrid.addChild(faceGrid);

		domain.addChild(treeGrid);
		root.addChild(domain);

		std::string filename = this->meshPrefix + "-mesh.xdmf";
		std::ofstream file(filename);
		file << root << std::endl;
		file.close(); 
	}

	// write grid for cell
	{	
		XdmfRoot root;
		XdmfDomain domain;
		XdmfGrid cellGrid = getXdmfCellGrid();
		domain.addChild(cellGrid);
		root.addChild(domain);
		std::ofstream file(this->meshPrefix + "-cell-mesh.xdmf");
		file << root << std::endl;
		file.close();  
	}
	
}

Vertex * Mesh::addVertex(const Eigen::Vector3d & position)
{
	Vertex * vertex = new Vertex(position, vertexList.size());
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

Edge * Mesh::getEdge(Vertex * A, Vertex * B) const
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
	Face * face0 = cellFaces[0];
	Face * face1 = cellFaces[1];
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

const std::vector<int> Mesh::getXdmfTopology_edge2vertexIndices() const {
	std::vector<int> data;
	const int nEdges = edgeList.size();
	for (int i = 0; i < nEdges; i++) {
		Edge* edge = edgeList[i];
		data.push_back(2); // type indicator of POLYLINE
		if (edge->getVertexMid() != nullptr) {
			data.push_back(3);  // this edge has three vertices
			data.push_back(edge->getVertexA()->getIndex());
			data.push_back(edge->getVertexMid()->getIndex());
			data.push_back(edge->getVertexB()->getIndex());
		}
		else {
			data.push_back(2); // this edge has two vertices
			data.push_back(edge->getVertexA()->getIndex());
			data.push_back(edge->getVertexB()->getIndex());
		}
	}
	return data;
}

const std::vector<int> Mesh::getXdmfTopology_face2vertexIndices() const
{
	std::vector<int> f2v;
	const int nFaces = getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		Face * face = faceList[i];
		std::vector<Vertex*> faceVerticesExtended = face->getVertexListExtended();
		f2v.push_back(3);
		f2v.push_back(faceVerticesExtended.size());
		for (Vertex* v : faceVerticesExtended) {
			f2v.push_back(v->getIndex());
		}
	}
	return f2v;
}

const std::vector<int> Mesh::getXdmfTopology_cell2vertexIndices() const
{
	std::vector<int> data;
	const int nCells = cellList.size();
	for (int i = 0; i < nCells; i++) {
		Cell * cell = cellList[i];
		data.push_back(16); // polyhedron type
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const int nCellFaces = cellFaces.size();
		data.push_back(nCellFaces); // number of faces
		for (int j = 0; j < nCellFaces; j++) {
			Face * face = cellFaces[j];
			const std::vector<Vertex*> faceVerticesExtended = face->getVertexListExtended();
			data.push_back(faceVerticesExtended.size()); // number of vertices
			for (auto vertex : faceVerticesExtended) {
				data.push_back(vertex->getIndex());
			}
		}
	}
	return data;
}

const std::vector<Cell*> Mesh::getCells() const
{
	return cellList;
}

const std::vector<Face*> Mesh::getFaces() const
{
	return faceList;
}

const std::vector<Edge*> Mesh::getEdges() const
{
	return edgeList;
}

const std::vector<Vertex*> Mesh::getVertices() const
{
	return vertexList;
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

void Mesh::facetCounting() {
	int nV_boundary = 0, nV_undefined = 0, nV_interior = 0, nV_electrode = 0, nV_insulating = 0,
	    nE_boundary = 0, nE_undefined = 0, nE_interior = 0, nE_electrode = 0, nE_insulating = 0,
		nF_boundary = 0, nF_undefined = 0, nF_interior = 0, nF_opening = 0, nF_wall = 0, 
		nC_undefined = 0, nC_fluid = 0, nC_solid = 0;
	for (Vertex* v : vertexList) {
		if (v->isBoundary()) nV_boundary++;
		switch(v->getType()) {
			case Vertex::Type::Undefined  : nV_undefined  ++; break;
			case Vertex::Type::Interior   : nV_interior   ++; break;
			case Vertex::Type::Electrode  : nV_electrode  ++; break;
			case Vertex::Type::Insulating : nV_insulating ++; break;
		}
	}
	for (Edge* e : edgeList) {
		if (e->isBoundary()) nE_boundary++;
		switch(e->getType()) {
			case Edge::Type::Undefined  : nE_undefined  ++; break;
			case Edge::Type::Interior   : nE_interior   ++; break;
			case Edge::Type::Electrode  : nE_electrode  ++; break;
			case Edge::Type::Insulating : nE_insulating ++; break;
		}
	}
	for (Face* f : faceList) {
		if (f->isBoundary()) nF_boundary++;
		switch(f->getFluidType()) {
			case Face::FluidType::Undefined : nF_undefined ++; break;
			case Face::FluidType::Interior  : nF_interior  ++; break;
			case Face::FluidType::Opening   : nF_opening   ++; break;
			case Face::FluidType::Wall      : nF_wall      ++; break;
		}
	}
	for (Cell* c : cellList) {
		switch(c->getFluidType()) {
			case Cell::FluidType::Undefined : nC_undefined ++; break;
			case Cell::FluidType::Fluid     : nC_fluid     ++; break;
			case Cell::FluidType::Solid     : nC_solid     ++; break;
		}
	}
	facet_counts.nV_undefined  = nV_undefined;
	facet_counts.nV_boundary   = nV_boundary;
	facet_counts.nV_interior   = nV_interior;
	facet_counts.nV_electrode  = nV_electrode;
	facet_counts.nV_insulating = nV_insulating;
	facet_counts.nE_undefined  = nE_undefined;
	facet_counts.nE_boundary   = nE_boundary;
	facet_counts.nE_interior   = nE_interior;
	facet_counts.nE_electrode  = nE_electrode;
	facet_counts.nE_insulating = nE_insulating;
	facet_counts.nF_undefined  = nF_undefined;
	facet_counts.nF_boundary   = nF_boundary;
	facet_counts.nF_interior   = nF_interior;
	facet_counts.nF_opening    = nF_opening;
	facet_counts.nF_wall       = nF_wall;
	facet_counts.nC_undefined  = nC_undefined;
	facet_counts.nC_fluid      = nC_fluid;
	facet_counts.nC_solid      = nC_solid;
	assert(facet_counts.nV_undefined + facet_counts.nV_interior + facet_counts.nV_electrode + facet_counts.nV_insulating == getNumberOfVertices());
	assert(facet_counts.nE_undefined + facet_counts.nE_interior + facet_counts.nE_electrode + facet_counts.nE_insulating == getNumberOfEdges());
	assert(facet_counts.nF_undefined + facet_counts.nF_interior + facet_counts.nF_opening + facet_counts.nF_wall == getNumberOfFaces());
	assert(facet_counts.nC_undefined + facet_counts.nC_fluid + facet_counts.nC_solid == getNumberOfCells());
	std::cout << "=======================================================" << std::endl;
	std::cout << meshPrefix + "-mesh:" << std::endl;
	std::cout << "      nV = "     << getNumberOfVertices() 
	          << "   boundary = "  << facet_counts.nV_boundary 
			  << "   undefined = " << facet_counts.nV_undefined 
			  << "   interior = "  << facet_counts.nV_interior
			  << "   electrode = " << facet_counts.nV_electrode
			  << "   insulating = "<< facet_counts.nV_insulating << std::endl;
	std::cout << "      nE = "     << getNumberOfEdges() 
	          << "   boundary = "  << facet_counts.nE_boundary 
			  << "   undefined = " << facet_counts.nE_undefined 
			  << "   interior = "  << facet_counts.nE_interior
			  << "   electrode = " << facet_counts.nE_electrode
			  << "   insulating = "<< facet_counts.nE_insulating << std::endl;
	std::cout << "      nF = "     << getNumberOfFaces() 
	          << "   boundary = "  << facet_counts.nF_boundary 
			  << "   undefined = " << facet_counts.nF_undefined 
			  << "   interior = "  << facet_counts.nF_interior
			  << "   opening = "   << facet_counts.nF_opening
			  << "   wall = "      << facet_counts.nF_wall << std::endl;
	std::cout << "      nC = "     << getNumberOfCells()  
			  << "   undefined = " << facet_counts.nC_undefined 
			  << "   fluid = "     << facet_counts.nC_fluid
			  << "   solid = "     << facet_counts.nC_solid << std::endl;
	std::cout << "=======================================================" << std::endl;
}

const Eigen::SparseMatrix<int>& Mesh::get_f2eMap() const
{
	return face2edgeMap;
}

const Eigen::SparseMatrix<int>& Mesh::get_e2vMap() const
{
	return edge2vertexMap;
}

bool Mesh::canBeContinuousLoop(std::vector<Edge*> edges) const {
	for (Edge* edge : edges) {
		int connected = 0;
		for (Edge* other : edges) {
			if (other != edge) {
				if (other->hasCoincidentVertex(edge)) {
					connected ++;
				}
			}
		}
		if (connected != 2) {
			for (Edge* e : edges) {
				std::cout << e->getVertexA()->getIndex() << "-" <<  e->getVertexB()->getIndex() << std::endl;
			}
			return false;
		}
	}
	return true;
}

std::vector<Edge*> Mesh::makeContinuousLoop(std::vector<Edge*> edges) const
{	
	assert(canBeContinuousLoop(edges));
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

bool Mesh::isConnected(Vertex * A, Vertex * B) const
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

void Mesh::update_vertexCoordinates()
{
	const int nVertices = vertexList.size();
	vertexCoordinates = Eigen::MatrixXd(nVertices, 3);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.row(i) = vertexList[i]->getPosition();
	}
}

void Mesh::createIncidenceMaps()
{
	create_edge2vertex_map();
	create_face2edge_map();
	create_cell2face_map();
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
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfVertices()}, 
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
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
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
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
					XdmfDataItem::Format::HDF
				),
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
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		vertexGrid.addChild(attribute);
	}
	return vertexGrid;
}

XdmfGrid Mesh::getXdmfEdgeGrid() const
{	
	H5Reader h5reader(this->meshPrefix + "-mesh.h5");
	XdmfGrid edgeGrid(XdmfGrid::Tags("Edge grid"));

	// Topology
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/edge2vertex";
		XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, getNumberOfEdges()));
		const int nElements = h5reader.readVectorDataSize("/edge2vertex");
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ nElements }, 
					XdmfDataItem::NumberType::Int, 
					XdmfDataItem::Format::HDF
				), 
				ss.str()
			)
		);
		edgeGrid.addChild(topology);
	}

	// Geometry
	{
		XdmfGeometry geometry;
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexPosExtended";
		geometry.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ (int)getVertexCoordinatesExtended().rows(), 3 },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF
					), 
				ss.str()
			)
		);
		edgeGrid.addChild(geometry);
	}

	// Attribute: edge index
	XdmfAttribute attribute(XdmfAttribute::Tags("Edge index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
	std::stringstream ss;
	ss << this->meshPrefix << "-mesh.h5:/edgeIdx";
	attribute.addChild(
		XdmfDataItem(
			XdmfDataItem::Tags(
				{ getNumberOfEdges() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF
			),
			ss.str()
		)
	);
	edgeGrid.addChild(attribute);

	// Attribute: edge type
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Edge Type", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/edgeType";
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfEdges() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		edgeGrid.addChild(attribute);
	}

	// Attribute: edge length
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Edge Length", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/edgeLength";
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfEdges() },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		edgeGrid.addChild(attribute);
	}


	return edgeGrid;
}

XdmfGrid Mesh::getXdmfFaceGrid() const
{
	H5Reader h5reader(this->meshPrefix + "-mesh.h5");

	XdmfGrid faceGrid(XdmfGrid::Tags("Face Grid"));

	// Topology
	XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, getNumberOfFaces()));  

	// Topology DataItem
	{
		const int nElements = h5reader.readVectorDataSize("/face2vertex");
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/face2vertex";
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ nElements }, 
					XdmfDataItem::NumberType::Int, 
					XdmfDataItem::Format::HDF
				), 
				ss.str()
			)
		);
		faceGrid.addChild(topology);
	}

	// Geometry and Geometry DataItem
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/vertexPosExtended";
		XdmfGeometry geometry;
		geometry.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ (int)getVertexCoordinatesExtended().rows(), 3 }, 
					XdmfDataItem::NumberType::Float, 
					XdmfDataItem::Format::HDF
				), 
				ss.str()
			)
		);
		
		faceGrid.addChild(geometry);
	}

	// Attribute: face index

	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/faceIdx";
		XdmfAttribute attribute(XdmfAttribute::Tags("Face Index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfFaces() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		faceGrid.addChild(attribute);
	}

	// Attribute: face type
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/faceType";
		XdmfAttribute attribute(XdmfAttribute::Tags("Face Type", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfFaces() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		faceGrid.addChild(attribute);
	}

	// Attribute: face area
	{
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/faceArea";
		XdmfAttribute attribute(XdmfAttribute::Tags("Face Area", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfFaces() },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		faceGrid.addChild(attribute);
	}
	
	return faceGrid;
}

XdmfGrid Mesh::getXdmfCellGrid() const
{
	H5Reader h5reader(this->meshPrefix + "-mesh.h5");
	XdmfGrid cellGrid(XdmfGrid::Tags("cellGrid"));

	// Topology
	XdmfTopology topology(XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, getNumberOfCells()));

	// Topology DataItem
	{
		const int nElements = h5reader.readVectorDataSize("/cell2vertex");
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/cell2vertex";
		topology.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ nElements },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		cellGrid.addChild(topology);
	}

	// Geometry and Geometry DataItem
	XdmfGeometry geometry;
	std::stringstream ss;
	ss << this->meshPrefix << "-mesh.h5:/vertexPosExtended";
	geometry.addChild(
		XdmfDataItem(
			XdmfDataItem::Tags(
				{ (int)getVertexCoordinatesExtended().rows(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF
			),
			ss.str()
		)
	);

	cellGrid.addChild(geometry);

	// Attribute: cell index
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Cell Index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
		ss << this->meshPrefix << "-mesh.h5:/cellIdx";
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfCells() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		cellGrid.addChild(attribute);
	}

	// Attribute: cell type
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Cell Type", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
        ss << this->meshPrefix << "-mesh.h5:/cellType";
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfCells() },
					XdmfDataItem::NumberType::Int,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		cellGrid.addChild(attribute);
	}

	// Attribute: cell volume
	{
		XdmfAttribute attribute(XdmfAttribute::Tags("Cell Volume", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell));
		std::stringstream ss;
        ss << this->meshPrefix << "-mesh.h5:/cellVolume";
		attribute.addChild(
			XdmfDataItem(
				XdmfDataItem::Tags(
					{ getNumberOfCells() },
					XdmfDataItem::NumberType::Float,
					XdmfDataItem::Format::HDF
				),
				ss.str()
			)
		);
		cellGrid.addChild(attribute);
	}

	return cellGrid;
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
		types(i) = static_cast<int>(getFace(i)->getFluidType());
	}
	return types;
}

const Eigen::VectorXi Mesh::getCellTypes() const
{
	const int nC = getNumberOfCells();
	Eigen::VectorXi types(nC);
	for (int i = 0; i < nC; i++) {
		types(i) = static_cast<int>(getCell(i)->getFluidType());
	}
	return types;
}

const Eigen::VectorXi Mesh::getVertexIndices() const {
	Eigen::VectorXi indices(getNumberOfVertices());
	for (int i = 0; i < getNumberOfVertices(); i++) {
		indices[i] = vertexList[i]->getIndex();
		assert(i == indices[i]);
	}
	return indices;
}

const Eigen::VectorXi Mesh::getEdgeIndices() const {
	Eigen::VectorXi indices(getNumberOfEdges());
	for (int i = 0; i < getNumberOfEdges(); i++) {
		indices[i] = edgeList[i]->getIndex();
		assert(i == indices[i]);
	}
	return indices;
}

const Eigen::VectorXi Mesh::getFaceIndices() const {
	Eigen::VectorXi indices(getNumberOfFaces());
	for (int i = 0; i < getNumberOfFaces(); i++) {
		indices[i] = faceList[i]->getIndex();
		assert(i == indices[i]);
	}
	return indices;
}

const Eigen::VectorXi Mesh::getCellIndices() const {
	Eigen::VectorXi indices(getNumberOfCells());
	for (int i = 0; i < getNumberOfCells(); i++) {
		indices[i] = cellList[i]->getIndex();
		assert(i == indices[i]);
	}
	return indices;
}

const Eigen::VectorXd Mesh::getEdgeLengths() const {
	Eigen::VectorXd lengths(getNumberOfEdges());
	for (int i = 0; i < getNumberOfEdges(); i++) {
		lengths[i] = edgeList[i]->getLength();
		assert(lengths[i] > 0);
	}
	return lengths;
}

const Eigen::VectorXd Mesh::getFaceAreas() const {
	Eigen::VectorXd areas(getNumberOfFaces());
	for (int i = 0; i < getNumberOfFaces(); i++) {
		areas[i] = faceList[i]->getArea();
		assert(areas[i] > 0);
	}
	return areas;
}

const Eigen::VectorXd Mesh::getCellVolumes() const {
	Eigen::VectorXd volumes(getNumberOfCells());
	for (int i = 0; i < getNumberOfCells(); i++) {
		volumes[i] = cellList[i]->getVolume();
		assert(volumes[i] > 0);
	}
	return volumes;
}

const Eigen::MatrixXd Mesh::getVertexCoordinatesExtended() const {
	Eigen::MatrixXd vertexCoordinatesExtended = vertexCoordinates;
	std::vector<Vertex*> aux_vertices;
	for (Edge* edge : edgeList) {
		if (edge->getVertexMid() != nullptr) {   // If this edge is not straight
			assert(edge->getVertexMid()->getIndex() >= getNumberOfVertices());
			aux_vertices.push_back(edge->getVertexMid());
		}
	}
	Eigen::size_t nVerticesExtended = getNumberOfVertices() + aux_vertices.size();
	vertexCoordinatesExtended.conservativeResize(nVerticesExtended, Eigen::NoChange);  // size = (nV + nV_aux, 3)
	for (Vertex* auxVertex : aux_vertices) {
		assert(auxVertex->getIndex() >= getNumberOfVertices());
		assert(auxVertex->getIndex() < getNumberOfVertices() + aux_vertices.size());
		vertexCoordinatesExtended.row(auxVertex->getIndex()) = auxVertex->getPosition();  
	}
	return vertexCoordinatesExtended;
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
//	const int nElements = h5reader.readVectorDataSize(dataname); 
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
