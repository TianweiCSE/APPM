#include "DualMesh.h"



DualMesh::DualMesh(PrimalMesh* primal)
	: Mesh("dual"), primal(primal)
{
}

DualMesh::DualMesh(const std::string & meshPrefix, PrimalMesh* primal) 
	: Mesh(meshPrefix), primal(primal)
{
}

DualMesh::~DualMesh()
{
}

Edge* DualMesh::addEdge(Edge* e1, Edge* e2){
	assert(e1->hasCoincidentVertex(e2));
	Edge* edge = new Edge(e1, e2);
	edge->setIndex(edgeList.size());
	edgeList.push_back(edge);
	return edge; 
}

void DualMesh::init_dualMesh()
{	
	assert(primal->getNumberOfCells() > 0);
	const int nPrimalVertices = primal->getNumberOfVertices();
	const int nPrimalEdges    = primal->getNumberOfEdges();
	const int nPrimalFaces    = primal->getNumberOfFaces();
	const int nPrimalCells    = primal->getNumberOfCells();
	const std::vector<Vertex*> primalVertices = primal->getVertices();
	const std::vector<Edge*>   primalEdges    = primal->getEdges();
	const std::vector<Face*>   primalFaces    = primal->getFaces();
	const std::vector<Cell*>   primalCells    = primal->getCells();

	/************************************************
	 * Add dual vertices
	 *    - at EVERY primal cell center (!identical index)
	 *    - at boundary primal face center
	 ************************************************/   
	// at primal cell center
	for (auto cell : primalCells) {
		addVertex(cell->getCenter());
	}
	// at primal boundary face center
	Eigen::SparseVector<int> primalFaceToDualVertex(nPrimalFaces);
	for (int i = 0; i < nPrimalFaces; i++) {
		const Face * face = primal->getFace(i);
		assert(face->getIndex() == i);
		if (face->isBoundary()) {
			const Vertex * V = addVertex(face->getCenter());
			primalFaceToDualVertex.coeffRef(i) = V->getIndex();
		}
	}

	// std::cout << "Dual vertices added." << std::endl;

	/******************************************************************
	 * Add dual edges
	 *    - penetrating EVERY primal face (!identical index & orientation)
	 *    - intersecting primal boundary edges
	 ******************************************************************/
	// Add dual edges penetrating primal faces 
	for (auto face : primalFaces) {
		const std::vector<Cell*> adjCells = face->getCellList();

		int idxA = -1;
		int idxB = -1;
		if (face->isBoundary()) {
			assert(adjCells.size() == 1);
			const Cell * cell = adjCells[0];
			const int orientation = cell->getOrientation(face);
			assert(abs(orientation) == 1);
			if (orientation > 0) {
				idxA = cell->getIndex();  // equal to the index of dual vertex 
				idxB = primalFaceToDualVertex.coeff(face->getIndex());
			}
			else {
				idxA = primalFaceToDualVertex.coeff(face->getIndex());
				idxB = cell->getIndex();
			}
		}
		else {
			assert(adjCells.size() == 2);
			int orientation = adjCells[0]->getOrientation(face);
			if (orientation > 0) {
				idxA = adjCells[0]->getIndex();
				idxB = adjCells[1]->getIndex();
			}
			else {
				idxA = adjCells[1]->getIndex();
				idxB = adjCells[0]->getIndex();
			}
		}
		assert(idxA >= 0);
		assert(idxB >= 0);
		Vertex * A = getVertex(idxA);
		Vertex * B = getVertex(idxB);
		addEdge(A, B);
	}

	/** Add auxiliary dual edges at the boundary
	 *  Need to tackle two cases:
	 *    - regular edge which is a straight segement
	 *    - irregular edge which has a joint
	 */
	int auxVertexIdx = getNumberOfVertices();   /* Assign indices to auxiliary vertices starting from NumberOfVertices.
												   Note that they will not be appended to vertexList. 
												   They need indexing only for the sake of outputing XDMF file to plot meshes */
	Eigen::SparseVector<int> primalEdgeToBoundaryDualEdge(nPrimalEdges);
	Edge* newEdge;
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primal->getEdge(i);
		assert(edge->getIndex() == i);
		if (edge->isBoundary()) {
			std::vector<Face*> boundaryFaces;
			for (auto face : edge->getFaceList()) {
				if (face->isBoundary()) {
					boundaryFaces.push_back(face);
				}
			}
			assert(boundaryFaces.size() == 2);
			Vertex* v1 = vertexList[primalFaceToDualVertex.coeffRef(boundaryFaces[0]->getIndex())];
			Vertex* v2 = vertexList[primalFaceToDualVertex.coeffRef(boundaryFaces[1]->getIndex())];
			Eigen::Vector3d fn0 = boundaryFaces[0]->getNormal();
			Eigen::Vector3d fn1 = boundaryFaces[1]->getNormal();
			if (fn0.cross(fn1).norm() > 4*std::numeric_limits<double>::epsilon()) {  // irregular edge
				Vertex* mid_v = new Vertex(edge->getHalfwayPosition(), auxVertexIdx++);
				Edge* e1 = new Edge(v1, mid_v);
				Edge* e2 = new Edge(v2, mid_v);
				newEdge = addEdge(e1, e2);
			}
			else {  // regular edge
				newEdge = addEdge(v1, v2);
			}
			primalEdgeToBoundaryDualEdge.coeffRef(i) = newEdge->getIndex(); // index of newly added edge
		}
	}
	// std::cout << "Dual edges added." << std::endl;

	/********************************************************************************
	 * Add dual faces 
	 * 		- being penetrated by EVERY primal edge (!identical index & orientation)
	 *            - composed of only inner edges
	 *            - composed of inner edges and one auxiliary edge
	 * 		- at boundary
	 *            - exclusively composed of auxiliary edges 
	 ********************************************************************************/
	/// Add dual faces being penetrated by primal edge
	std::vector<Edge*> dualEdges;
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * primalEdge = primal->getEdge(i);
		const std::vector<Face*> primalEdgeFaces = primalEdge->getFaceList();
		for (auto primalFace : primalEdgeFaces) {
				const int pFidx = primalFace->getIndex();  // identical to the index of associated dual edge 
				dualEdges.push_back(getEdge(pFidx));
		}

		if (primalEdge->isBoundary()) {  // For boundary case, auxiliary edges should be considered 
			Edge* edge_aux = getEdge(primalEdgeToBoundaryDualEdge.coeffRef(i));
			dualEdges.push_back(edge_aux);
		} 
		std::vector<Edge*> dualEdgesLoop = makeContinuousLoop(dualEdges);
		Face * dualFace = addFace(dualEdgesLoop);
		dualEdges.clear();

		// Check orientation of dual face normal and primal edge; if necessary, flip dual normal
		Eigen::Vector3d fn = dualFace->getNormal();
		const Eigen::Vector3d primalEdgeDir = primalEdge->getDirection();
		if (fn.dot(primalEdgeDir) < 0) {
			dualFace->reverseNormal();
		}
	}

	// Add dual faces at boundary (one to one with primal boundary vertices)
	Face * newFace;
	Eigen::SparseVector<int> primalVertexToDualBoundaryFace(nPrimalVertices);
	for (Vertex* primalVertex : primalVertices) {
		if (primalVertex->isBoundary()) {
			const std::vector<Edge*> vertexEdges = primalVertex->getEdges();
			for (Edge* vertexEdge : vertexEdges) {
				if (vertexEdge->isBoundary()) {
					dualEdges.push_back(getEdge(primalEdgeToBoundaryDualEdge.coeffRef(vertexEdge->getIndex())));
				}
			}
			std::vector<Edge*> dualEdgesLoop = makeContinuousLoop(dualEdges);
			newFace = addFace(dualEdgesLoop);
			dualEdges.clear();
			primalVertexToDualBoundaryFace.coeffRef(primalVertex->getIndex()) = newFace->getIndex(); // index of newly added face
		}
	}

	// std::cout << "Dual faces added." << std::endl;

	/*************************************************************************
	 * Add dual cells (!identical index with associated primal vertices)
	 *************************************************************************/
	std::vector<Face*> dualFaces;
	for (const Vertex* primalVertex : primalVertices) {
		for (const Edge* primalVertexEdge : primalVertex->getEdges()) {
			dualFaces.push_back(getFace(primalVertexEdge->getIndex()));  // Note each primal edge has the same index with its corresponding dual face
		}
		if (primalVertex->isBoundary()) {  // If the vertex is at boundary, auxiliary face is needed to close the dual cell
			dualFaces.push_back(getFace(primalVertexToDualBoundaryFace.coeffRef(primalVertex->getIndex())));
		}
		addCell(dualFaces);
		assert(primalVertex->getIndex() == cellList.size() - 1);
		dualFaces.clear();
	}

	// std::cout << "Dual cells added." << std::endl;

	update_vertexCoordinates();
	createIncidenceMaps();
	init_cellFluidType();
	init_faceFluidType();
	facetCounting();
}

void DualMesh::init_cellFluidType()
{
	const int nCells = this->getNumberOfCells();
	for (int i = 0; i < nCells; i++) {
		Cell::FluidType fluidType;
		Cell * cell = getCell(i);
		const Eigen::Vector3d cellCenter = cell->getCenter();
		const Eigen::Vector2d cellCenter2d = cellCenter.segment(0, 2);

		if (cellCenter2d.norm() < 1) {    /// The fluid radius = 1 
			fluidType = Cell::FluidType::Fluid;
		}
		else {
			fluidType = Cell::FluidType::Solid;
		}
		cell->setFluidType(fluidType);
	}
}

void DualMesh::init_faceFluidType()
{
	//const double terminalRadius = 0.35;
	const int nFaces = this->getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		Face * face = getFace(i);
		Face::FluidType faceFluidType = Face::FluidType::Undefined;

		std::vector<Cell*> faceCells = face->getCellList();
		const int nAdjacientCells = faceCells.size();
		assert(nAdjacientCells >= 1 && nAdjacientCells <= 2);

		int nSolidCells = 0;
		int nFluidCells = 0;
		for (auto cell : faceCells) {
			if (cell->getFluidType() == Cell::FluidType::Fluid) {
				nFluidCells++;
			}
			if (cell->getFluidType() == Cell::FluidType::Solid) {
				nSolidCells++;
			}
		}
		assert(nSolidCells >= 0 && nSolidCells <= faceCells.size());
		assert(nFluidCells >= 0 && nFluidCells <= faceCells.size());
		assert(nFluidCells + nSolidCells == faceCells.size());

		if (nFluidCells == 1 && nSolidCells == 1) {
			faceFluidType = Face::FluidType::Wall;
		}
		else if (nFluidCells == 1 && nSolidCells == 0) {
			faceFluidType = Face::FluidType::Opening;
		}
		else if (nFluidCells == 2) {
			faceFluidType = Face::FluidType::Interior;
		}
		else {
			faceFluidType = Face::FluidType::Undefined;
		}

		face->setFluidType(faceFluidType);
	}
}

void DualMesh::check() const {

	/*****************************************************
	 *  Check:
	 *    - primal_f2e == dual_f2e^t
	 *    - primal_e2v == (-1)*dual_c2f^t
	 *    - numbers of facets
	 *****************************************************/
	const int nPrimalVertices = primal->getNumberOfVertices();
	const int nPrimalEdges    = primal->getNumberOfEdges();
	const int nPrimalFaces    = primal->getNumberOfFaces();
	const int nPrimalCells    = primal->getNumberOfCells();

	// primal face-to-edge map = transpose of dual face-to-edge map  /// Dual mesh has extra edges and faces on the boundary!
	std::cout << "- Checking primal_f2e == dual_f2e^t ------------ ";
	const Eigen::SparseMatrix<int> & p_f2e = primal->get_f2eMap();
	Eigen::SparseMatrix<int> temp = face2edgeMap.topLeftCorner(nPrimalEdges, nPrimalFaces).transpose(); 
	assert(p_f2e.rows() == temp.rows() && p_f2e.cols() == temp.cols());
	Eigen::SparseMatrix<int> delta = (p_f2e - temp).pruned();// pruned() keeps only non-zero matrix entries
	std::cout << (delta.nonZeros() == 0 ? "[PASSED]" : "[FAILED]") << std::endl;
	//assert(delta.nonZeros() == 0);

	// primal edge-to-vertex map = (-1) * transpose of dual cell-to-face map  /// Dual mesh has extra faces on the boundary!
	std::cout << "- Checking primal_e2v == (-1)*dual_c2f^t ------- ";
	const Eigen::SparseMatrix<int> & p_e2v = primal->get_e2vMap();
	temp = -1 * cell2faceMap.topLeftCorner(nPrimalVertices, nPrimalEdges).transpose();
	assert(p_e2v.rows() == temp.rows() && p_e2v.cols() == temp.cols());
	delta = (p_e2v - temp).pruned(); // pruned() keeps only non-zero matrix entries
	std::cout << (delta.nonZeros() == 0 ? "[PASSED]" : "[FAILED]") << std::endl; 
	//assert(delta.nonZeros() == 0);

	// Check the amounts relations between primal and dual facets
	assert(facet_counts.nC_undefined == 0);
	assert(facet_counts.nV_boundary == primal->facet_counts.nF_boundary);
	assert(facet_counts.nE_boundary == primal->facet_counts.nE_boundary);
	assert(facet_counts.nF_boundary == primal->facet_counts.nV_boundary);
	assert(getNumberOfVertices() - facet_counts.nV_boundary == primal->getNumberOfCells());
	assert(getNumberOfEdges() - facet_counts.nE_boundary == primal->getNumberOfFaces());
	assert(getNumberOfFaces() - facet_counts.nF_boundary == primal->getNumberOfEdges());
	assert(getNumberOfCells() == primal->getNumberOfVertices());
	std::cout << "- Facets number checking --------- [PASSED]" << std::endl; 
}
