#include "PrimalMesh.h"



PrimalMesh::PrimalMesh() :
	PrimalMesh("primal")
{
	//std::cout << "Call to PrimalMesh()" << std::endl;
}

PrimalMesh::PrimalMesh(const std::string & meshPrefix)
	: Mesh(meshPrefix)
{
	//std::cout << "Call to PrimalMesh(string)" << std::endl;
}

PrimalMesh::PrimalMesh(const PrimalMeshParams & p)
	: PrimalMesh()
{
	this->params = p;
}


PrimalMesh::~PrimalMesh()
{
	//std::cout << "Call to ~PrimalMesh()" << std::endl;
}

void PrimalMesh::init_cube()
{	
	electrodeGeo = PrimalMesh::ElectrodeGeometry::Square;
	validateParameters();
	const double zmax = params.getZmax();
	const double z0 = -0.5 * zmax;
	const int num_vertex_per_side = 4 * params.getRefinements() + 1;
	const double increment = 3.0 / (num_vertex_per_side - 1);

	Eigen::Vector3d unit_x; unit_x << 1., 0., 0.;
	Eigen::Vector3d unit_y; unit_y << 0., -1., 0.;
	Eigen::Vector3d topLeft; topLeft << -1.5, 1.5, z0;

	for (int i = 0; i < num_vertex_per_side; i++) {
		for (int j = 0; j < num_vertex_per_side; j++) {
			addVertex(topLeft + i*increment*unit_x + j*increment*unit_y);
		}
	}
	for (int i = 0; i < num_vertex_per_side - 1; i++) {
		for (int j = 0; j < num_vertex_per_side - 1; j++) {
			Vertex* topLeft_v  = vertexList[j*num_vertex_per_side + i];
			Vertex* topRight_v = vertexList[j*num_vertex_per_side + i + 1];
			Vertex* botLeft_v  = vertexList[(j+1)*num_vertex_per_side + i];
			Vertex* botRight_v  = vertexList[(j+1)*num_vertex_per_side + i + 1]; 
			auto edge1 = addEdge(topLeft_v,  topRight_v);
			auto edge2 = addEdge(topLeft_v,  botLeft_v);
			auto edge3 = addEdge(botRight_v, botLeft_v);
			auto edge4 = addEdge(botRight_v, topRight_v);
			std::vector<Edge*> faceEdges = {edge1, edge2, edge3, edge4};
			addFace(faceEdges);
		}
	}

	vertexCoordinates = Eigen::MatrixXd(3, getNumberOfVertices());
	for (int k = 0; k < getNumberOfVertices(); k++) {
		vertexCoordinates.col(k) = getVertex(k)->getPosition();
	}
	
	check_zCoord(z0);

	std::cout << "outerMeshExtrude completed." << std::endl; 
	
	const int axialLayers = params.getAxialLayers();
	if (axialLayers == 0) {
		return;
	}
	extrudeMesh(axialLayers, zmax);

	sortVertices();
	sortEdges();
	sortFaces();
	sortCells();
	facetCounting();
}

void PrimalMesh::init_cylinder()
{	
	electrodeGeo = PrimalMesh::ElectrodeGeometry::Round;
	validateParameters();
	const double zmax = params.getZmax();

	const double z0 = -0.5 * zmax;
	init_hexagon(z0);
	check_zCoord(z0);

	assert(getNumberOfVertices() > 0);
	refineMesh(params.getRefinements());
	check_zCoord(z0);
	
	std::cout << "refineMesh completed." << std::endl;

	outerMeshExtrude(params.getOuterLayers());
	check_zCoord(z0);

	std::cout << "outerMeshExtrude completed." << std::endl; 
	
	const int axialLayers = params.getAxialLayers();
	if (axialLayers == 0) {
		return;
	}
	extrudeMesh(axialLayers, zmax);

	sortVertices();
	sortEdges();
	sortFaces();
	sortCells();
	facetCounting();
}

/*
void PrimalMesh::init()
{
	validateParameters();
	const double zmax = params.getZmax();
	const double z0 = -0.5 * zmax;
	const int num_vertex_per_side = 4 * params.getRefinements() + 1;
	const double increment = 1.0 / (num_vertex_per_side + 1);
	
	Eigen::Vector3d unit_x; unit_x << 1., 0., 0.;
	Eigen::Vector3d unit_y; unit_y << 0., 1., 0.;
	
	// Generate the inner part first
	Eigen::Vector3d topLeft; topLeft << -0.5 + increment, 0.5 - increment, z0;
	for (int i = 0; i < num_vertex_per_side; i++) {
		for (int j = 0; j < num_vertex_per_side; j++) {
			addVertex(topLeft + j*increment*unit_x - i*increment*unit_y);
		}
	}
	for (int i = 0; i < num_vertex_per_side - 1; i++) {
		for (int j = 0; j < num_vertex_per_side - 1; j++) {
			Vertex* topLeft_v  = vertexList[i*num_vertex_per_side + j];
			Vertex* topRight_v = vertexList[i*num_vertex_per_side + j + 1];
			Vertex* botLeft_v  = vertexList[(i+1)*num_vertex_per_side + j];
			Vertex* botRight_v  = vertexList[(i+1)*num_vertex_per_side + j + 1]; 
			auto edge1 = addEdge(topLeft_v,  topRight_v);
			auto edge2 = addEdge(topLeft_v,  botLeft_v);
			auto edge3 = addEdge(botRight_v, botLeft_v);
			auto edge4 = addEdge(botRight_v, topRight_v);
			std::vector<Edge*> faceEdges = {edge1, edge2, edge3, edge4};
			addFace(faceEdges);
		}
	}

	// Then, generate the outer layer which has only one layer of cell
	// This is due to the convenction we follow when generating the dual mesh: No corner cell!
	
	addVertex(getVertex(0)->getPosition() - increment*unit_x + increment*unit_y);
	for (int j = 1; j < num_vertex_per_side; j++) {
		Vertex* v = getVertex(j);
		Vertex* v_u;
		if (j != num_vertex_per_side - 1) {
			v_u = addVertex(v->getPosition() + increment * unit_y);
		}
		else {
			v_u = addVertex(v->getPosition() + increment*unit_x + increment*unit_y);
		}
		Vertex* v_l    = getVertex(j-1);
		Vertex* v_lu   = getVertex(v_u->getIndex() - 1);
		Edge* e1 = addEdge(v, v_u);
		Edge* e2 = addEdge(v_u, v_lu);
		Edge* e3 = addEdge(v_lu, v_l);
		Edge* e4 = addEdge(v_l, v);
		addFace(std::vector<Edge*>({e1,e2,e3,e4}));  
	}
	for (int i = 1; i < num_vertex_per_side; i++) {
		Vertex* v = getVertex(i*num_vertex_per_side + num_vertex_per_side - 1);
		Vertex* v_u;
		if (i != num_vertex_per_side - 1) {
			v_u = addVertex(v->getPosition() + increment * unit_x);
		}
		else {
			v_u = addVertex(v->getPosition() + increment*unit_x - increment*unit_y);
		}
		Vertex* v_l    = getVertex((i-1)*num_vertex_per_side + num_vertex_per_side - 1);
		Vertex* v_lu   = getVertex(v_u->getIndex() - 1);
		Edge* e1 = addEdge(v, v_u);
		Edge* e2 = addEdge(v_u, v_lu);
		Edge* e3 = addEdge(v_lu, v_l);
		Edge* e4 = addEdge(v_l, v);
		addFace(std::vector<Edge*>({e1,e2,e3,e4}));  
	}
	for (int j = 1; j < num_vertex_per_side; j++) {
		Vertex* v = getVertex(num_vertex_per_side*num_vertex_per_side - 1 - j);
		Vertex* v_u;
		if (j != num_vertex_per_side - 1) {
			v_u = addVertex(v->getPosition() - increment * unit_y);
		}
		else {
			v_u = addVertex(v->getPosition() - increment*unit_x - increment*unit_y);
		}
		Vertex* v_l    = getVertex(num_vertex_per_side*num_vertex_per_side - j);
		Vertex* v_lu   = getVertex(v_u->getIndex() - 1);
		Edge* e1 = addEdge(v, v_u);
		Edge* e2 = addEdge(v_u, v_lu);
		Edge* e3 = addEdge(v_lu, v_l);
		Edge* e4 = addEdge(v_l, v);
		addFace(std::vector<Edge*>({e1,e2,e3,e4}));  
	}
	for (int i = 1; i < num_vertex_per_side - 1; i++) {
		Vertex* v = getVertex(num_vertex_per_side*(num_vertex_per_side - 1 - i));
		Vertex* v_u;
		v_u = addVertex(v->getPosition() - increment * unit_x);
		Vertex* v_l    = getVertex(num_vertex_per_side*(num_vertex_per_side - i));
		Vertex* v_lu   = getVertex(v_u->getIndex() - 1);
		Edge* e1 = addEdge(v, v_u);
		Edge* e2 = addEdge(v_u, v_lu);
		Edge* e3 = addEdge(v_lu, v_l);
		Edge* e4 = addEdge(v_l, v);
		addFace(std::vector<Edge*>({e1,e2,e3,e4}));  
	}
	{
		Vertex* v   = getVertex(0);
		Vertex* v_u = getVertex(num_vertex_per_side*num_vertex_per_side);
		Vertex* v_l = getVertex(num_vertex_per_side);
		Vertex* v_lu = getVertex(getNumberOfVertices() - 1);
		Edge* e1 = addEdge(v, v_u);
		Edge* e2 = addEdge(v_u, v_lu);
		Edge* e3 = addEdge(v_lu, v_l);
		Edge* e4 = addEdge(v_l, v);
		addFace(std::vector<Edge*>({e1,e2,e3,e4}));
	}

	vertexCoordinates = Eigen::MatrixXd(3, getNumberOfVertices());
	for (int k = 0; k < getNumberOfVertices(); k++) {
		vertexCoordinates.col(k) = getVertex(k)->getPosition();
	}
	

	check_zCoord(z0);

	std::cout << "outerMeshExtrude completed." << std::endl; 
	
	const int axialLayers = params.getAxialLayers();
	if (axialLayers == 0) {
		return;
	}
	extrudeMesh(axialLayers, zmax);

	sortVertices();
	sortEdges();
	sortFaces();
	sortCells();
	facetCounting();
}*/

void PrimalMesh::check() {
	assert(facet_counts.nV_undefined == 0);
	assert(facet_counts.nV_boundary == facet_counts.nV_electrode + facet_counts.nV_insulating);
	assert(facet_counts.nE_undefined == 0);
	assert(facet_counts.nE_boundary == facet_counts.nE_electrode + facet_counts.nE_insulating);
}

void PrimalMesh::init_hexagon(const double zValue)
{
	std::cout << "Initialize with hexagon" << std::endl;
	Vertex * origin = addVertex(Eigen::Vector3d(0, 0, zValue));
	const int corners = 6;
	for (int k = 0; k < corners; k++) {
		const double phi = 2 * M_PI * k / corners;
		const Eigen::Vector3d pos(cos(phi), sin(phi), zValue);
		Vertex * v = addVertex(pos);
		addEdge(origin, v);
	}
	const int nV = getNumberOfVertices();
	vertexCoordinates = Eigen::MatrixXd(3, nV);
	for (int k = 0; k < nV; k++) {
		vertexCoordinates.col(k) = getVertex(k)->getPosition();
	}

	for (int k = 0; k < corners; k++) {
		Vertex * A = getVertex(k + 1);
		Vertex * B = getVertex(((k + 1) % corners) + 1);
		auto edgeA = addEdge(origin, A);
		auto edgeB = addEdge(A, B);
		auto edgeC = addEdge(origin, B);
		std::vector<Edge*> faceEdges = { edgeA, edgeB, edgeC };
		Face * face = addFace(faceEdges);
		//Face * face = addFace({origin, A, B});
	}
}

void PrimalMesh::init_triangle()
{
	assert(false);
	//Vertex * origin = addVertex(Eigen::Vector3d(0, 0, 0));
	//Vertex * A = addVertex(Eigen::Vector3d(1, 0, 0));
	//Vertex * B = addVertex(Eigen::Vector3d(0.5, 1, 0));
	//Edge * edge0 = addEdge(origin, A);
	//Edge * edge1 = addEdge(A, B);
	//Edge * edge2 = addEdge(origin, B);
	//addFace({ edge0, edge1, edge2 });
}

void PrimalMesh::refineMesh(const int nRefinements)
{
	// Strategy: 
	// construct refined f2v-map from current mesh 
	// and reconstruct from that map

	for (int level = 0; level < nRefinements; level++) {
		std::cout << "Mesh refinement level: " << level << std::endl;
		Eigen::Matrix3Xi f2v;

		// Get face-to-vertex list of refined faces
		if (level != 1) {
			f2v = refine_triangles();
		}
		else {
			f2v = refine_triangles_specialCorners();
		}
		//std::cout << "f2v: " << std::endl;
		//std::cout << f2v.transpose() << std::endl;
		assert((f2v.array() >= 0).all());

		// std::ofstream file;
		// std::stringstream ss;
		/*
		{
			std::stringstream ss;
			ss << "level" << level << "-coords.dat";
			std::ofstream file = std::ofstream(ss.str());
			file << vertexCoordinates.transpose() << std::endl;
		}
		// ss = std::stringstream();
		{
			std::stringstream ss;
			ss << "level" << level << "-f2v.dat";
			std::ofstream file = std::ofstream(ss.str());
			file << f2v.transpose() << std::endl;
		}
		*/
		//if (level == 1) { break; }

		const double tol = 16 * std::numeric_limits<double>::epsilon();
		const double zPosMax = vertexCoordinates.row(2).maxCoeff();
		const double zPosMin = vertexCoordinates.row(2).minCoeff();
		assert(std::abs(zPosMax - zPosMin) < tol);

		// clear mesh elements
		for (auto v : vertexList) {
			delete v;
		}
		vertexList.resize(0);
		for (auto e : edgeList) {
			delete e;
		}
		edgeList.resize(0);
		for (auto f : faceList) {
			delete f;
		}
		faceList.resize(0);

		// reconstruct mesh
		int nVertices = vertexCoordinates.cols();
		assert(vertexCoordinates.rows() == 3);
		for (int i = 0; i < nVertices; i++) {
			Eigen::Vector3d pos = vertexCoordinates.col(i);
			addVertex(pos);
		}
		assert(f2v.rows() == 3);
		int nFaces = f2v.cols();
		for (int i = 0; i < nFaces; i++) {
			Eigen::Vector3i triangleIdx = f2v.col(i);
			std::vector<Edge*> faceEdges;
			for (int j = 0; j < 3; j++) {
				assert((triangleIdx.array() >= 0).all());
				const int idxA = triangleIdx(j);
				const int idxB = triangleIdx((j + 1) % 3);

				Vertex * A = getVertex(idxA);
				Vertex * B = getVertex(idxB);
				Edge * edge = addEdge(A, B);
				faceEdges.push_back(edge);
			}
			//std::cout << "face idx: " << i << std::endl;
			//std::cout << "Edges: ";
			//for (auto edge : faceEdges) {
			//	std::cout << *edge << std::endl;
			//}
			Face * face = addFace(faceEdges);
			//Vertex * A = getVertex(triangleIdx(0));
			//Vertex * B = getVertex(triangleIdx(1));
			//Vertex * C = getVertex(triangleIdx(2));
			//Face * face = addFace({A, B, C});
		}

		// Check
		for (auto face : faceList) {
			assert(face->getEdgeList().size() == 3);
			assert(face->getVertexListExtended().size() == 3);
		}
	}

}

void PrimalMesh::outerMeshExtrude(const int nLayers)
{
	for (int layer = 0; layer < nLayers; layer++) {
		//outerMeshExtrude_triangles();
		outerMeshExtrude_prisms();
	}
}

void PrimalMesh::outerMeshExtrude_triangles()
{
	const int nVertices = vertexList.size();
	assert(nVertices == vertexCoordinates.cols());

	const int nEdges = edgeList.size();

	Eigen::SparseVector<int> boundaryVertices(nVertices);
	Eigen::SparseVector<int> boundaryEdges(nEdges);
	Eigen::SparseVector<int> boundaryEdgeVertexIdx(nEdges);   

	// get boundary vertices and boundary edges
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		if (edge->isBoundary()) {
			boundaryEdges.coeffRef(i) = 1;
			const int idxA = edge->getVertexA()->getIndex();
			const int idxB = edge->getVertexB()->getIndex();
			boundaryVertices.coeffRef(idxA) = 1;
			boundaryVertices.coeffRef(idxB) = 1;
		}
	}

	// create outer vertices for each boundary edge
	const int nnzE = boundaryEdges.nonZeros();
	Eigen::Matrix3Xd boundaryEdgeVertexCoord(3, nnzE);
	const int offset_E = vertexCoordinates.cols();
	int idx = 0;
	for (int i = 0; i < nnzE; i++) {
		const int idxE = boundaryEdges.innerIndexPtr()[i];
		const Edge * edge = edgeList[idxE];
		assert(edge != nullptr);
		Eigen::Vector3d pos = edge->getHalfwayPosition();
		pos.segment(0,2) *= 1.25;
		boundaryEdgeVertexCoord.col(idx++) = pos;
		boundaryEdgeVertexIdx.coeffRef(edge->getIndex()) = i + offset_E;   

		const Vertex * V = addVertex(pos);
		assert(V->getIndex() == i + offset_E);
	}
	assert(idx == nnzE);
	vertexCoordinates.conservativeResize(3, offset_E + nnzE);
	vertexCoordinates.rightCols(nnzE) = boundaryEdgeVertexCoord;

	// create outer vertices for each boundary vertex
	const int nnzV = boundaryVertices.nonZeros();
	Eigen::SparseVector<int> boundaryVertexVertexIdx(nVertices);
	Eigen::Matrix3Xd boundaryVertexVertexCoord(3, nnzV);
	const int offset_V = vertexCoordinates.cols();
	idx = 0;
	for (int i = 0; i < nnzV; i++) {
		const int idxV = boundaryVertices.innerIndexPtr()[i];
		boundaryVertexVertexIdx.coeffRef(idxV) = i + offset_V;
		const Vertex * vertex = getVertex(idxV);
		Eigen::Vector3d pos = vertex->getPosition();
		pos.segment(0,2) *= 1.5;
		boundaryVertexVertexCoord.col(idx++) = pos;

		const Vertex * V = addVertex(pos);
		assert(V->getIndex() == i + offset_V);
	}
	assert(idx == nnzV);
	vertexCoordinates.conservativeResize(3, offset_V + nnzV);
	vertexCoordinates.rightCols(nnzV) = boundaryVertexVertexCoord;

	// create faces at boundary edge
	Eigen::Matrix3Xi f2v_edges(3, 2*nnzE);
	f2v_edges.array() = -1;
	idx = 0;
	for (int i = 0; i < nnzE; i++) {
		const int idxE = boundaryEdges.innerIndexPtr()[i];
		Edge * edge = edgeList[idxE];
		const int idxC = boundaryEdgeVertexIdx.coeff(edge->getIndex());
		Vertex * C1 = getVertex(idxC);
		assert(C1 != nullptr);
		Vertex * A = edge->getVertexA();
		Vertex * B = edge->getVertexB();

		// check for face orientation; ensure that the face resulting face has positive orientation in z-direction
		/// ???
		const Eigen::Vector3d a = A->getPosition();
		const Eigen::Vector3d b = B->getPosition() - A->getPosition();
		const Eigen::Vector3d n = a.normalized().cross(b.normalized());
		const double orientation = n.dot(Eigen::Vector3d::UnitZ());
		if (orientation < 0) {
			Vertex * temp = A;
			A = B;
			B = temp;
		}
		f2v_edges.col(idx++) = Eigen::Vector3i(A->getIndex(), C1->getIndex(), B->getIndex());
	
		const int idxA1 = boundaryVertexVertexIdx.coeff(A->getIndex());
		const int idxB1 = boundaryVertexVertexIdx.coeff(B->getIndex());
		Vertex * A1 = getVertex(idxA1);
		Vertex * B1 = getVertex(idxB1);
		f2v_edges.col(idx++) = Eigen::Vector3i(A1->getIndex(), B1->getIndex(), C1->getIndex());
	}
	assert(f2v_edges.cols() == idx);

	// create faces at boundary vertices
	Eigen::Matrix3Xi f2v_vertices(3, 2*nnzV);
	f2v_vertices.array() = -1;
	idx = 0;
	for (int i = 0; i < nnzV; i++) {
		const int idxV = boundaryVertices.innerIndexPtr()[i];
		Vertex * B = getVertex(idxV);

		// get adjacient boundary edges
		const std::vector<Edge*> vertexEdges = B->getEdges();
		std::vector<Edge*> boundaryEdgeList;
		for (auto edge : vertexEdges) {
			if ((edge->getIndex() < nEdges) && edge->isBoundary()) {
				boundaryEdgeList.push_back(edge);
			}
		}
		assert(boundaryEdgeList.size() == 2);

		Vertex * C1 = getVertex(boundaryEdgeVertexIdx.coeff(boundaryEdgeList[0]->getIndex()));
		Vertex * D1 = getVertex(boundaryEdgeVertexIdx.coeff(boundaryEdgeList[1]->getIndex()));
		Vertex * B1 = getVertex(boundaryVertexVertexIdx.coeff(B->getIndex()));

		// check for face orientation
		const Eigen::Vector3d a = C1->getPosition() - B->getPosition();
		const Eigen::Vector3d b = D1->getPosition() - C1->getPosition();
		const Eigen::Vector3d n = a.normalized().cross(b.normalized());
		const double orientation = n.dot(Eigen::Vector3d::UnitZ());
		if (orientation < 0) {
			Vertex * temp = C1;
			C1 = D1;
			D1 = temp;
		}

		f2v_vertices.col(idx++) = Eigen::Vector3i(B->getIndex(), C1->getIndex(), D1->getIndex());
		f2v_vertices.col(idx++) = Eigen::Vector3i(B1->getIndex(), D1->getIndex(), C1->getIndex());
	}
	assert(f2v_vertices.cols() == idx);

	// add faces
	for (int i = 0; i < f2v_edges.cols(); i++) {
		Vertex * A = getVertex(f2v_edges(0, i));
		Vertex * B = getVertex(f2v_edges(1, i));
		Vertex * C = getVertex(f2v_edges(2, i));
		//addFace({ A, B, C });
		addFace({ addEdge(A, B), addEdge(B, C), addEdge(C, A) });
	}
	for (int i = 0; i < f2v_vertices.cols(); i++) {
		Vertex * A = getVertex(f2v_vertices(0, i));
		Vertex * B = getVertex(f2v_vertices(1, i));
		Vertex * C = getVertex(f2v_vertices(2, i));
		addFace({ addEdge(A, B), addEdge(B, C), addEdge(C, A) });
		//addFace({ A, B, C });
	}
}

void PrimalMesh::outerMeshExtrude_prisms()
{
	const int nVertices = getNumberOfVertices();
	assert(nVertices == vertexCoordinates.cols());
	const int nEdges = getNumberOfEdges();

	Eigen::SparseVector<int> boundaryVertices(nVertices);
	Eigen::SparseVector<int> boundaryEdges(nEdges);

	// collect boundary vertices and edges
	for (int i = 0; i < nEdges; i++) {
		const Edge * edge = getEdge(i);
		if (edge->isBoundary()) {
			boundaryEdges.coeffRef(i) = 1;
			const Vertex * A = edge->getVertexA();
			const int idxA = edge->getVertexA()->getIndex();
			const int idxB = edge->getVertexB()->getIndex();
			boundaryVertices.coeffRef(idxA) = 1;
			boundaryVertices.coeffRef(idxB) = 1;
		}
	}

	// create outer vertices for each boundary vertex
	Eigen::SparseVector<int> newVertexIdx(nVertices);
	const int nbV = boundaryVertices.nonZeros();
	Eigen::Matrix3Xd newPosCoords(3, nbV);
	for (int i = 0; i < nbV; i++) {
		const int idxV = boundaryVertices.innerIndexPtr()[i];
		const Vertex * vertex = getVertex(idxV);
		const Eigen::Vector3d posBoundary = vertex->getPosition();
		const Eigen::Vector2d posBoundary_2d = posBoundary.segment(0, 2);
		Eigen::Vector2d pos_2d = posBoundary_2d + 0.5 / params.getOuterLayers() * posBoundary_2d.normalized();
		Eigen::Vector3d newPos(pos_2d(0), pos_2d(1), posBoundary(2));
		newPosCoords.col(i) = newPos;
		Vertex * V = addVertex(newPos);
		newVertexIdx.coeffRef(idxV) = V->getIndex();
		assert(V->getIndex() == nVertices + i);
	}
	const int offset = vertexCoordinates.cols();
	vertexCoordinates.conservativeResize(3, offset + nbV);
	vertexCoordinates.rightCols(nbV) = newPosCoords;

	// create prism face for each boundary edge
	const int nbE = boundaryEdges.nonZeros();
	Eigen::MatrixXi f2v_edges(4, nbE);
	f2v_edges.array() = -1;
	for (int i = 0; i < nbE; i++) {
		Edge * edge = getEdge(boundaryEdges.innerIndexPtr()[i]);
		Vertex * A = edge->getVertexA();
		Vertex * B = edge->getVertexB();

		int idxC = newVertexIdx.coeff(B->getIndex());
		int idxD = newVertexIdx.coeff(A->getIndex());
		assert(idxC >= nVertices);
		assert(idxD >= nVertices);

		Vertex * C = getVertex(idxC);
		Vertex * D = getVertex(idxD);

		// store face-to-vertex indices: (A,B,C,D)
		Eigen::VectorXi col(4);
		col(0) = A->getIndex();
		col(1) = B->getIndex();
		col(2) = C->getIndex();
		col(3) = D->getIndex();
		f2v_edges.col(i) = col;
	}
	
	// add faces
	for (int i = 0; i < f2v_edges.cols(); i++) {
		Vertex * A = getVertex(f2v_edges(0, i));
		Vertex * B = getVertex(f2v_edges(1, i));
		Vertex * C = getVertex(f2v_edges(2, i));
		Vertex * D = getVertex(f2v_edges(3, i));
		Face * face = addFace({ addEdge(A, B), addEdge(B, C), addEdge(D, C), addEdge(A, D)});
	}
}

void PrimalMesh::extrudeMesh(const int nLayers, const double zmax)
{
	if (nLayers <= 0) {
		return;
	}
	std::cout << "Extrude mesh with " << nLayers << " layers" << std::endl;
	assert(zmax > 0);
	const Eigen::Vector3d z_unit(0, 0, 1);
	const int nVertices_2d = vertexList.size();
	const int nEdges_2d = edgeList.size();
	const int nFaces_2d = faceList.size();

	// Reserve memory
	cellList.reserve(nFaces_2d * nLayers);

	// Create vertices
	for (int layer = 1; layer <= nLayers; layer++) {
		for (int i = 0; i < nVertices_2d; i++) {
			Eigen::Vector3d pos = getVertex(i)->getPosition() + layer * ((zmax/nLayers) * z_unit);
			addVertex(pos);
		}
	}

	for (int layer = 1; layer <= nLayers; layer++) {
		std::cout << "Mesh layer: " << layer << std::endl;
		// Create edges: parallel to z_unit
		for (int i = 0; i < nVertices_2d; i++) {
			const Vertex * vRef = getVertex(i);
			Vertex * A = getVertex(vRef->getIndex() + (layer - 1) * nVertices_2d);
			Vertex * B = getVertex(vRef->getIndex() + (layer - 0) * nVertices_2d);
			addEdge(A, B);
		}
		// Create edges: normal to z_unit
		for (int i = 0; i < nEdges_2d; i++) {
			const Edge * edgeRef = getEdge(i);
			const Vertex * Aref = edgeRef->getVertexA();
			const Vertex * Bref = edgeRef->getVertexB();
			Vertex * A = getVertex(Aref->getIndex() + layer * nVertices_2d);
			Vertex * B = getVertex(Bref->getIndex() + layer * nVertices_2d);
			addEdge(A, B);
		}
		// Create cells from 2d mesh
		for (int i = 0; i < nFaces_2d; i++) {
			const Face * faceRef = getFace(i);
			const std::vector<Edge*> edgeListRef = faceRef->getEdgeList();
			//const std::vector<Vertex*> vListRef = faceRef->getVertexList();

			std::vector<Edge*> bottomEdges;
			std::vector<Edge*> topEdges;
			std::vector<Face*> sideFaces;
			for (int k = 0; k < edgeListRef.size(); k++) {
				const Edge * edgeRef = edgeListRef[k];

				Vertex * bottomA = getVertex(edgeRef->getVertexA()->getIndex() + (layer - 1) * nVertices_2d);
				Vertex * bottomB = getVertex(edgeRef->getVertexB()->getIndex() + (layer - 1) * nVertices_2d);
				Edge * bottomEdge = addEdge(bottomA, bottomB);
				bottomEdges.push_back(bottomEdge);

				Vertex * topA = getVertex(edgeRef->getVertexA()->getIndex() + (layer - 0) * nVertices_2d);
				Vertex * topB = getVertex(edgeRef->getVertexB()->getIndex() + (layer - 0) * nVertices_2d);
				Edge * topEdge = addEdge(topA, topB);
				topEdges.push_back(topEdge);

				Edge * sideEdgeA = addEdge(bottomA, topA);
				Edge * sideEdgeB = addEdge(bottomB, topB);
				Face * sideFace = addFace({bottomEdge, sideEdgeB, topEdge, sideEdgeA});
				sideFaces.push_back(sideFace);
			}
			//std::cout << "bottom edges: " << std::endl;
			//for (auto edge : bottomEdges) {
			//	std::cout << *edge << std::endl;
			//}
			Face * bottomFace = addFace(bottomEdges);
			//std::cout << "bottom face center: " << bottomFace->getCenter().transpose() << std::endl;
			Face * topFace = addFace(topEdges);

			assert(bottomFace != nullptr);
			assert(topFace != nullptr);
			addTriPrism(sideFaces, bottomFace, topFace);
		}
	}
}

Eigen::Matrix3Xi PrimalMesh::refine_triangles()
{
	const int nVertices = vertexList.size();
	assert(nVertices > 0);

	// copy vertices
	vertexCoordinates = Eigen::MatrixXd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}

	// define edge midpoint vertices
	const int nEdges = edgeList.size();
	Eigen::Matrix3Xd edgeMidpoints(3, nEdges);
	for (int i = 0; i < nEdges; i++) {
		const Edge * edge = edgeList[i];
		Eigen::Vector3d pos = edge->getHalfwayPosition();
		Eigen::Vector2d pos2d = pos.segment(0, 2);
		if (edge->isBoundary()) {
			//pos.segment(0,2) /= pos.segment(0,2).norm();
			pos2d /= pos2d.norm();
		}
		pos.segment(0, 2) = pos2d;
		edgeMidpoints.col(i) = pos;
	}
	vertexCoordinates.conservativeResize(3, nVertices + nEdges);
	vertexCoordinates.rightCols(nEdges) = edgeMidpoints;

	// Refine faces into four subfaces, defined by edge midpoints
	const int nFaces = faceList.size();
	Eigen::Matrix3Xi f2v(3, 4*nFaces);
	f2v.array() = -1;

	int idx = 0;
	for (int i = 0; i < nFaces; i++) {
		Face * face = faceList[i];
		std::vector<Edge*> faceEdges = face->getEdgeList();
		assert(faceEdges.size() == 3);
		Edge * edge0 = faceEdges[0];
		Edge * edge1 = faceEdges[1];
		Edge * edge2 = faceEdges[2];

		assert(edge1->hasCoincidentVertex(edge2));
		Vertex * A = edge1->getCoincidentVertex(edge2);
		assert(edge2->hasCoincidentVertex(edge0));
		Vertex * B = edge2->getCoincidentVertex(edge0);
		assert(edge0->hasCoincidentVertex(edge1));
		Vertex * C = edge0->getCoincidentVertex(edge1);

		int v0 = A->getIndex();
		int v1 = B->getIndex();
		int v2 = C->getIndex();
		int e0 = faceEdges[0]->getIndex() + nVertices;
		int e1 = faceEdges[1]->getIndex() + nVertices;
		int e2 = faceEdges[2]->getIndex() + nVertices;

		f2v.col(idx++) = Eigen::Vector3i(v0, e2, e1);
		f2v.col(idx++) = Eigen::Vector3i(v1, e0, e2);
		f2v.col(idx++) = Eigen::Vector3i(v2, e1, e0);
		f2v.col(idx++) = Eigen::Vector3i(e0, e1, e2);
	}
	assert(idx == f2v.cols());
	assert((f2v.array() >= 0).all());
	return f2v;
}

Eigen::Matrix3Xi PrimalMesh::refine_triangles_specialCorners()
{
	// TODO: check that vertices are created at zMaxPos! (create in 2d!)

	const int nVertices = vertexList.size();
	assert(nVertices > 0);
	const int nEdges = edgeList.size();
	const int nFaces = faceList.size();
	Eigen::Matrix3Xi f2v(3, 0);


	// copy vertices
	vertexCoordinates = Eigen::MatrixXd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}

	// Indicator for special edges (i.e., those that are split in half instead of thirds)
	Eigen::VectorXi specialEdges = Eigen::VectorXi::Zero(nEdges);
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		const int nA = edge->getVertexA()->getEdges().size();
		const int nB = edge->getVertexB()->getEdges().size();
		specialEdges(i) = (nA == 3 || nB == 3);
	}

	// Create vertices at one-third and two-third of each edge;
	// except for special edges, where we create vertex at midpoint.
	// We keep track of the mapping from edge to new vertices.
	Eigen::Matrix3Xd edgeInnerCoords(3, 2 * nEdges - specialEdges.array().sum());
	assert(edgeInnerCoords.cols() > nEdges);

	Eigen::Matrix2Xi edge2innerVertices(2, nEdges);
	edge2innerVertices.array() = -1;
	int idx = 0;
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		const Eigen::Vector3d A = edge->getVertexA()->getPosition();
		const Eigen::Vector3d B = edge->getVertexB()->getPosition();
		const Eigen::Vector3d edgeVector = edge->getDirection() * edge->getLength();

		if (specialEdges(i)) {
			Eigen::Vector3d pos = A + 0.5 * edgeVector;
			if (edge->isBoundary()) {
				Eigen::Vector2d pos2d = pos.segment(0, 2);
				pos2d.normalize();
				pos.segment(0, 2) = pos2d;
			}
			edgeInnerCoords.col(idx) = pos;
			edge2innerVertices(0, i) = nVertices + idx++;
		}
		else {
			const Eigen::Vector3d pos0 = A + 1. / 3. * edgeVector;
			edgeInnerCoords.col(idx) = pos0;
			edge2innerVertices(0, i) = nVertices + idx++;
			const Eigen::Vector3d pos1 = A + 2. / 3. * edgeVector;
			edgeInnerCoords.col(idx) = pos1;
			edge2innerVertices(1, i) = nVertices + idx++;
		}
	}
	vertexCoordinates.conservativeResize(3, nVertices + edgeInnerCoords.cols());
	vertexCoordinates.rightCols(edgeInnerCoords.cols()) = edgeInnerCoords;

	// new vertices at face center 
	int faceOffset = vertexCoordinates.cols();
	Eigen::Matrix3Xd fc(3, nFaces);
	for (int i = 0; i < nFaces; i++) {
		Eigen::Vector3d center;
		center.setZero();
		const Face * face = faceList[i];

		// face center by arithmetic mean of vertices
		std::vector<Vertex*> faceVertices = face->getVertexListExtended();
		assert(face->getVertexListExtended().size() == face->getVertexList().size());
		assert(faceVertices.size() == 3);
		for (int j = 0; j < 3; j++) {
			center += faceVertices[j]->getPosition();
		}
		center /= 3;
		fc.col(i) = center;
	}
	vertexCoordinates.conservativeResize(3, faceOffset + nFaces);
	vertexCoordinates.rightCols(nFaces) = fc;

	// For each face ...
	f2v = Eigen::Matrix3Xi(3, 9 * nFaces);
	f2v.array() = -1;
	int fidx = 0;
	for (int i = 0; i < nFaces; i++) {
		const Face * face = faceList[i];
		const std::vector<Edge*> faceEdges = face->getEdgeList();
		if (face->hasBoundaryEdges()) {
			// ... boundary face ...
			// find local index of not-special edge
			int idx_notSpecialEdge = -1;
			for (int idx = 0; idx < 3; idx++) {
				if (specialEdges(faceEdges[idx]->getIndex()) == 0) {
					idx_notSpecialEdge = idx;
					break;
				}
			}

			// Vertex index at face center
			int vc = face->getIndex() + faceOffset;
			// Position of face center
			const Eigen::Vector3d center = vertexCoordinates.col(vc);

			assert(idx_notSpecialEdge >= 0);
			// edge0 is not a special edge, but edge1 and edge2 are a special edge
			Edge * edge0 = faceEdges[idx_notSpecialEdge];
			Edge * edge1 = faceEdges[(idx_notSpecialEdge + 1) % 3];
			Edge * edge2 = faceEdges[(idx_notSpecialEdge + 2) % 3];
			assert(specialEdges(edge0->getIndex()) == 0);
			assert(specialEdges(edge1->getIndex()) == 1);
			assert(specialEdges(edge2->getIndex()) == 1);

			// Ensure that the triangle {edge0,edge1,edge2} is a right-handed system
			const Eigen::Vector3d temp1 = edge1->getHalfwayPosition();
			const Eigen::Vector3d temp2 = edge2->getHalfwayPosition();
			const Eigen::Vector3d a = temp1 - center;
			const Eigen::Vector3d b = temp2 - temp1;
			const Eigen::Vector3d n = a.normalized().cross(b.normalized());
			const double tol = 16 * std::numeric_limits<double>::epsilon();
			assert(n.segment(0, 2).norm() < tol);
			// If the normal vector is anti-parallel to z-axis, flip edge1 and edge2
			if (n.dot(Eigen::Vector3d::UnitZ()) < 0) {
				Edge * tempEdge = edge1;
				edge1 = edge2;
				edge2 = tempEdge;
			}

			int v0 = edge1->getCoincidentVertex(edge2)->getIndex();
			int v1 = edge2->getCoincidentVertex(edge0)->getIndex();
			int v2 = edge0->getCoincidentVertex(edge1)->getIndex();

			// points at 1/3 and 2/3 of edge0 
			const int incidence = edge0->getIncidence(getVertex(v1));
			const int e0 = edge0->getIndex();
			int v01 = (incidence == 1) ? edge2innerVertices(0, e0) : edge2innerVertices(1, e0);
			int v02 = (incidence == 1) ? edge2innerVertices(1, e0) : edge2innerVertices(0, e0);
			assert(v01 >= 0);
			assert(v02 >= 0);

			// Midpoint vertices on special edges
			int v11 = edge2innerVertices(0, edge1->getIndex());
			int v22 = edge2innerVertices(0, edge2->getIndex());


			// Define refined faces
			f2v.col(fidx++) = Eigen::Vector3i(vc, v01, v02);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v02, v11);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v11, v0);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v0, v22);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v22, v01);
			f2v.col(fidx++) = Eigen::Vector3i(v1, v01, v22);
			f2v.col(fidx++) = Eigen::Vector3i(v2, v11, v02);
		}
		else {
			// ... not boundary face ...
			Edge * edge0 = faceEdges[0];
			Edge * edge1 = faceEdges[1];
			Edge * edge2 = faceEdges[2];

			const int v0 = edge1->getCoincidentVertex(edge2)->getIndex();
			const int v1 = edge2->getCoincidentVertex(edge0)->getIndex();
			const int v2 = edge0->getCoincidentVertex(edge1)->getIndex();
			
			const int incidence0 = edge0->getIncidence(getVertex(v1));
			const int incidence1 = edge1->getIncidence(getVertex(v0));
			const int incidence2 = edge2->getIncidence(getVertex(v0));

			const int e0 = edge0->getIndex();
			const int e1 = edge1->getIndex();
			const int e2 = edge2->getIndex();
			int v01 = (incidence0 == 1) ? edge2innerVertices(0, e0) : edge2innerVertices(1, e0);
			int v02 = (incidence0 == 1) ? edge2innerVertices(1, e0) : edge2innerVertices(0, e0);
			int v10 = (incidence1 == 1) ? edge2innerVertices(0, e1) : edge2innerVertices(1, e1);
			int v12 = (incidence1 == 1) ? edge2innerVertices(1, e1) : edge2innerVertices(0, e1);
			int v20 = (incidence2 == 1) ? edge2innerVertices(0, e2) : edge2innerVertices(1, e2);
			int v21 = (incidence2 == 1) ? edge2innerVertices(1, e2) : edge2innerVertices(0, e2);

			const int vc = face->getIndex() + faceOffset;
			// Define refined faces
			f2v.col(fidx++) = Eigen::Vector3i(vc, v10, v20);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v20, v21);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v21, v01);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v01, v02);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v02, v12);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v12, v10);
			f2v.col(fidx++) = Eigen::Vector3i(v0, v20, v10);
			f2v.col(fidx++) = Eigen::Vector3i(v1, v01, v21);
			f2v.col(fidx++) = Eigen::Vector3i(v2, v12, v02);
		}
	}
	f2v.conservativeResize(3, fidx);
	assert((f2v.array() >= 0).all());
	return f2v; 
}

void PrimalMesh::test_quadFace()
{
	Vertex * A = addVertex(Eigen::Vector3d(0, 0, 0));
	Vertex * B = addVertex(Eigen::Vector3d(1, 0, 0));
	Vertex * C = addVertex(Eigen::Vector3d(1, 1, 0));
	Vertex * D = addVertex(Eigen::Vector3d(0, 1, 0));
	Edge * edge0 = addEdge(A, B);
	Edge * edge1 = addEdge(B, C);
	Edge * edge2 = addEdge(C, D);
	Edge * edge3 = addEdge(D, A);
	addFace({ edge0, edge1, edge2, edge3 });
	//addFace({ A, B, C, D });
}


/**
 * Sort vertices into order: inner, electrode, insulate vertices.
 */
void PrimalMesh::sortVertices()
{
	std::vector<Vertex*> interiorVertices;
	std::vector<Vertex*> electrodeVertices;
	std::vector<Vertex*> insulatingVertices;

	const int nV = getNumberOfVertices();
	vertexCoordinates = Eigen::MatrixXd(3, nV);
	for (int i = 0; i < nV; i++) {
		vertexCoordinates.col(i) = getVertex(i)->getPosition();
	}

	const double zmax = vertexCoordinates.row(2).maxCoeff();
	const double zmin = vertexCoordinates.row(2).minCoeff();

	std::cout << "z-axis range: [" << zmin << ", " << zmax << "]" << std::endl;

	for (Vertex* vertex : vertexList) {
		if (vertex->isBoundary()) {
			const Eigen::Vector3d pos = vertex->getPosition();
			const Eigen::Vector2d pos_2d(pos.segment(0, 2));

			bool is_zmax = std::abs(pos(2) - zmax) < 100 * std::numeric_limits<double>::epsilon();
			bool is_zmin = std::abs(pos(2) - zmin) < 100 * std::numeric_limits<double>::epsilon();

			if (electrodeGeo == PrimalMesh::ElectrodeGeometry::Square) {
				if ((pos_2d.cwiseAbs().array() < params.getElectrodeRadius()).all() && (is_zmax || is_zmin) ) {
					electrodeVertices.push_back(vertex);
					vertex->setType(Vertex::Type::Electrode);
				}
				else {
					insulatingVertices.push_back(vertex);
					vertex->setType(Vertex::Type::Insulating);
				}
			}
			else { 
				assert(electrodeGeo == PrimalMesh::ElectrodeGeometry::Round);
				if (pos_2d.norm() < params.getElectrodeRadius() && (is_zmax || is_zmin) ) {
					electrodeVertices.push_back(vertex);
					vertex->setType(Vertex::Type::Electrode);
				}
				else {
					insulatingVertices.push_back(vertex);
					vertex->setType(Vertex::Type::Insulating);
				}
			}
		}
		else {
			interiorVertices.push_back(vertex);
			vertex->setType(Vertex::Type::Interior);
		}
	}

	// Define new list of vertices
	std::vector<Vertex*> & sortedVertexList = interiorVertices; 
	sortedVertexList.insert(sortedVertexList.end(), electrodeVertices.begin(),  electrodeVertices.end());
	sortedVertexList.insert(sortedVertexList.end(), insulatingVertices.begin(), insulatingVertices.end());
	// Reassign the index
	for (int i = 0; i < sortedVertexList.size(); i++) {
		sortedVertexList[i]->setIndex(i);
	}
	assert(vertexList.size() == sortedVertexList.size());
	vertexList = sortedVertexList;
}

/** 
* Sort edges into order: Interior, Electrode, Insulating
*/
void PrimalMesh::sortEdges()
{
	std::vector<Edge*> interiorEdges;
	std::vector<Edge*> electrodeEdges;
	std::vector<Edge*> insulatingEdges;

	for (Edge* edge : edgeList) {
		bool isBoundaryA = edge->getVertexA()->isBoundary();
		bool isBoundaryB = edge->getVertexB()->isBoundary();
		Vertex::Type typeA = edge->getVertexA()->getType();
		Vertex::Type typeB = edge->getVertexB()->getType(); 

		// if (isBoundaryA && isBoundaryB) { // This condition is wrong!
		if (edge->isBoundary()) {  // if boundary edge
			if (typeA == Vertex::Type::Electrode && typeB == Vertex::Type::Electrode) {
				electrodeEdges.push_back(edge);
				edge->setType(Edge::Type::Electrode);
			}
			else {
				insulatingEdges.push_back(edge);
				edge->setType(Edge::Type::Insulating);
			}
		}
		else {
			interiorEdges.push_back(edge);
			edge->setType(Edge::Type::Interior);
		}
	}
	
	// Define new list of edges
	std::vector<Edge*> & sortedEdgeList = interiorEdges;
	sortedEdgeList.insert(sortedEdgeList.end(), electrodeEdges.begin(),  electrodeEdges.end());
	sortedEdgeList.insert(sortedEdgeList.end(), insulatingEdges.begin(), insulatingEdges.end());
	// Reassign the indices
	for (int i = 0; i < sortedEdgeList.size(); i++) {
		sortedEdgeList[i]->setIndex(i);
	}
	assert(sortedEdgeList.size() == edgeList.size());
	edgeList = sortedEdgeList;
}


/**
 * Sort faces into order: Interior, Electrode, Insulating
 */
void PrimalMesh::sortFaces()
{
	std::vector<Face*> interiorFaces;
	std::vector<Face*> electrodeFaces;
	std::vector<Face*> insulatingFaces;
	for (Face* face : faceList) {
		if (face->isBoundary()) {
			if (face->getCenter().segment(0,2).norm() < params.getElectrodeRadius()) {
				electrodeFaces.push_back(face);
			}
			else {
				insulatingFaces.push_back(face);
			}
		}
		else {
			interiorFaces.push_back(face);
		}
	}
	std::vector<Face*> & sortedFaceList = interiorFaces;
	sortedFaceList.insert(sortedFaceList.end(), electrodeFaces.begin(), electrodeFaces.end());
	sortedFaceList.insert(sortedFaceList.end(), insulatingFaces.begin(), insulatingFaces.end());
	for (int i = 0; i < sortedFaceList.size(); i++) {
		sortedFaceList[i]->setIndex(i);
	}
	assert(sortedFaceList.size() == faceList.size());
	faceList = sortedFaceList;
}

void PrimalMesh::sortCells()
{
	std::vector<Cell*> innerCells;
	std::vector<Cell*> outerCells;
	const int nCells = getNumberOfCells();
	for (int i = 0; i < nCells; i++) {
		Cell * cell = getCell(i);
		const Eigen::Vector3d cellCenter = cell->getCenter();
		const Eigen::Vector2d cc_2d = cellCenter.segment(0, 2);
		if (cc_2d.norm() < 1) {     /// The radius of plasma domain is one by default
			innerCells.push_back(cell);
		}
		else {
			outerCells.push_back(cell);
		}
	}
	std::vector<Cell*> sortedCells(nCells);
	int offset = 0;
	for (auto cell : innerCells) {
		sortedCells[offset++] = cell;
	}
	for (auto cell : outerCells) {
		sortedCells[offset++] = cell;
	}
	assert(offset == nCells);
	for (int i = 0; i < sortedCells.size(); i++) {
		sortedCells[i]->setIndex(i);
	}
	cellList = sortedCells;
}


void PrimalMesh::validateParameters()
{
	std::cout << "Primal mesh parameters: " << std::endl;
	std::cout << params << std::endl;
	assert(params.getAxialLayers() >= 0);
	if (params.getAxialLayers() == 0) {
		std::cout << "Warning: axialLayers is zero, a 2d mesh is created." << std::endl;
	}
	assert(params.getRefinements() >= 1);
	if (params.getOuterLayers() > 0) {
		assert(params.getRefinements() > 1);
	}
	assert(params.getZmax() > 0);
}

/** 
* Check that all vertices have z-coordinate equal to z0.
*/
void PrimalMesh::check_zCoord(const double z0)
{
	const double zmaxValue = vertexCoordinates.row(2).maxCoeff();
	const double zminValue = vertexCoordinates.row(2).minCoeff();

	const double tol = 16 * std::numeric_limits<double>::epsilon();
	assert(std::abs(zmaxValue - z0) < tol);
	assert(std::abs(zminValue - z0) < tol);
}

int PrimalMesh::count_electrode_vertices() const {
	int count = 0;
	for (Vertex* v : vertexList) {
		if (v->getType() == Vertex::Type::Electrode) {
			count++;
		}
	}
	return count;
}

int PrimalMesh::count_insulating_vertices() const {
	int count = 0;
	for (Vertex* v : vertexList) {
		if (v->getType() == Vertex::Type::Insulating) {
			count++;
		}
	}
	return count;
}

int PrimalMesh::count_interior_faces() const {
	int count = 0;
	for (Face* f : faceList) {
		if (!f->isBoundary()) {
			count++;
		}
	}
	return count;
}

PrimalMesh::PrimalMeshParams::PrimalMeshParams()
{
}

PrimalMesh::PrimalMeshParams::PrimalMeshParams(const std::string & filename)
{
	readParameters(filename);
}

int PrimalMesh::PrimalMeshParams::getRefinements() const
{
	return nRefinements;
}

int PrimalMesh::PrimalMeshParams::getAxialLayers() const
{
	return nAxialLayers;
}

int PrimalMesh::PrimalMeshParams::getOuterLayers() const
{
	return nOuterLayers;
}

double PrimalMesh::PrimalMeshParams::getElectrodeRadius() const
{
	return electrodeRadius;
}

double PrimalMesh::PrimalMeshParams::getZmax() const
{
	return zmax;
}

void PrimalMesh::PrimalMeshParams::readParameters(const std::string & filename)
{
	if (filename.size() <= 0) {
		std::cout << "Filename of PrimalMeshParameters not valid: " << filename << std::endl;
	}
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File is not open: " << filename << std::endl;
	}

	std::string line;
	const char delim = ':';
	while (std::getline(file, line)) {
		std::cout << line << std::endl;

		// find position of delimiter
		int pos = line.find(delim);
		// if delimiter is not found, skip this line
		if (pos == std::string::npos) {
			std::cout << "delimiter (delim = " << delim << ") not found in line: " << line << "; skip this line" << std::endl;
			continue;
		}
		std::string tag = line.substr(0, pos);
		
		if (tag == "axialLayers") {
			std::stringstream(line.substr(pos + 1)) >> this->nAxialLayers;
		}
		if (tag == "refinements") {
			std::stringstream(line.substr(pos + 1)) >> this->nRefinements;
		}
		if (tag == "outerLayers") {
			std::stringstream(line.substr(pos + 1)) >> this->nOuterLayers;
		}
		if (tag == "zMax") {
			std::stringstream(line.substr(pos + 1)) >> this->zmax;
		}
		if (tag == "electrodeRadius") {
			std::stringstream(line.substr(pos + 1)) >> this->electrodeRadius;
		}
	}
}

std::ostream & operator<<(std::ostream & os, const PrimalMesh::PrimalMeshParams & obj)
{
	os << "electrode radius: " << obj.electrodeRadius << std::endl;
	os << "refinements:      " << obj.nRefinements << std::endl;
	os << "axial layers:     " << obj.nAxialLayers << std::endl;
	os << "outer layers:     " << obj.nOuterLayers << std::endl;
	os << "zMax:             " << obj.zmax << std::endl;
	return os;
}
