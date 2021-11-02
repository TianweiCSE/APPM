#include "MaxwellSolver.h"



MaxwellSolver::MaxwellSolver()
{
}

MaxwellSolver::MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual)
{
	this->primal = primal;
	this->dual = dual;
	std::cout << "Dual mesh has " << dual->getNumberOfVertices() << " vertices" << std::endl;

	init();
}

MaxwellSolver::~MaxwellSolver()
{
}

//void MaxwellSolver::updateMaxwellState(const double dt, const double time)
//{
//	std::cout << "You should call the inherited function, not this one" << std::endl;
//}

void MaxwellSolver::writeStates(H5Writer & writer) const
{
	const int nPrimalFaces = primal->getNumberOfFaces();
	const int nDualEdges = dual->getNumberOfEdges();
	const int nDualFaces = dual->getNumberOfFaces();

	writer.writeDoubleVector(B_h, "/bvec");
	writer.writeDoubleVector(E_h, "/evec");
	writer.writeDoubleVector(H_h, "/hvec");
	writer.writeDoubleVector(J_h, "/jvec");

	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd H(3, nDualEdges);
	Eigen::Matrix3Xd J(3, nDualFaces);

	for (int i = 0; i < nPrimalFaces; i++) {
		const Eigen::Vector3d fn = primal->getFace(i)->getNormal();
		const double fA = primal->getFace(i)->getArea();
		B.col(i) = (B_h(i) / fA) * fn;
	}
	for (int i = 0; i < nDualFaces; i++) {
		const Eigen::Vector3d fn = dual->getFace(i)->getNormal();
		const double fA = dual->getFace(i)->getArea();
		J.col(i) = J_h(i) / fA * fn;
	}
	writer.writeDoubleMatrix(B, "/B");
	writer.writeDoubleMatrix(J, "/J");

	const int nPrimalEdges = primal->getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primal->getEdge(i);
		E.col(i) = E_h(i) / edge->getLength() * edge->getDirection();
	}
	writer.writeDoubleMatrix(E, "/E");

	for (int i = 0; i < nDualEdges; i++) {
		const Edge * edge = dual->getEdge(i);
		H.col(i) = H_h(i) / edge->getLength() * edge->getDirection();
	}
	writer.writeDoubleMatrix(H, "/H");
}

const Eigen::VectorXd & MaxwellSolver::getBstate() const
{
	return B_h;
}

void MaxwellSolver::setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector)
{
	const int nFaces = primal->getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		const Face * face = primal->getFace(i);
		const double area = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		B_h(i) = area * fn.dot(fieldVector);
	}
}

void MaxwellSolver::setAzimuthalMagneticFluxField()
{
	const int nFaces = primal->getNumberOfFaces();
	for (int i = 0; i < nFaces; i++) {
		const Face * face = primal->getFace(i);
		const double fa = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		if (fn.cross(Eigen::Vector3d::UnitZ()).norm() > 0.99) {
			// azimuthal field
			Eigen::Vector3d fc = face->getCenter();
			fc(2) = 0;
			fc.normalize();
			B_h(i) = fn.dot(fc.cross(Eigen::Vector3d::UnitZ()));
		}
		else {
			B_h(i) = 0;
		}
	}
}

void MaxwellSolver::setTorusCurrent(const double x1, const double x2, const double z1, const double z2)
{
	std::cout << "Set torus current" << std::endl;
	const Eigen::Vector3d center = Eigen::Vector3d(0.0, 0.0, 0.5);

	const double tol = 8 * std::numeric_limits<double>::epsilon();

	const int nCells = dual->getNumberOfCells();

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dual->getCell(i);
		const Eigen::Vector3d center = cell->getCenter();
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const Eigen::Matrix3Xd cellVertexCoords = cell->getVertexCoordinates();
		const Eigen::ArrayXd xv = cellVertexCoords.row(0).array(); // x-coordinates of cell vertices
		const Eigen::ArrayXd yv = cellVertexCoords.row(1).array(); // 1-coordinates of cell vertices
		const Eigen::ArrayXd zv = cellVertexCoords.row(2).array(); // z-coordinates of cell vertices
		const bool isXZplane = (yv > 0).any() && (yv < 0).any();

		// Skip cells that are not across xz-plane
		if (!isXZplane || (xv < x1).all() || (xv > x2).all() || (zv < z1).all() || (zv > z2).all()) { continue; }

		// For each face of this cell
		for (auto face : cellFaces) {
			const int faceIdx = face->getIndex();
			const Eigen::Vector3d fn = face->getNormal();

			// Skip faces that have a normal vector out of xz plane
			if (abs(Eigen::Vector3d::UnitY().dot(fn)) > tol) { continue; }

			// Get vertex coordinates of this face
			const std::vector<Vertex*> faceVertices = face->getVertexList();
			Eigen::Matrix3Xd faceVertexPos(3, faceVertices.size());
			for (int i = 0; i < faceVertices.size(); i++) {
				faceVertexPos.col(i) = faceVertices[i]->getPosition();
			}
			const Eigen::ArrayXd x = faceVertexPos.row(0).array();
			const Eigen::ArrayXd z = faceVertexPos.row(2).array();
			const bool isXnormal = fn.cross(Eigen::Vector3d::UnitX()).norm() < tol;
			const bool isZnormal = fn.cross(Eigen::Vector3d::UnitZ()).norm() < tol;

			// Determine if face is pierced by one of the loop edges:
			// edge 1: z = z1, x = [x1,x2]
			bool isPiercedByEdge1 = isXnormal && (z <= (z1 + tol)).any() && (x >= x1).all() && (x <= x2).all();
			if (isPiercedByEdge1) {
				J_h(faceIdx) = fn.dot(Eigen::Vector3d::UnitX());
			}

			// edge 2: x = x1, z = [z1,z2]
			bool isPiercedByEdge2 = isZnormal && (x <= (x1 + tol)).any() && (z >= z1).all() && (z <= z2).all();
			if (isPiercedByEdge2) {
				J_h(faceIdx) = -fn.dot(Eigen::Vector3d::UnitZ());
			}

			// edge 3: z = z2, x = [x1,x2]
			bool isPiercedByEdge3 = isXnormal && (z >= (z2 - tol)).any() && (x >= x1).all() && (x <= x2).all();
			if (isPiercedByEdge3) {
				J_h(faceIdx) = -fn.dot(Eigen::Vector3d::UnitX());
			}

			// edge 4: x = x2, z = [z1,z2]
			bool isPiercedByEdge4 = isZnormal && (x >= (x2 - tol)).any() && (z >= z1).all() && (z <= z2).all();
			if (isPiercedByEdge4) {
				J_h(faceIdx) = fn.dot(Eigen::Vector3d::UnitZ());
			}
		}
	}
}

Eigen::VectorXd MaxwellSolver::electricPotentialTerminals(const double time)
{
	const int n = primalMeshInfo.nVerticesTerminal;
	assert(n % 2 == 0);
	const double t0 = 3;
	const double sigma_t = 1;
	const double phiA = 0;
	const double phiB = 0;
	const double phi1 = phiA * 0.5 * (1 + tanh((time - t0) / sigma_t));
	const double phi2 = phiB;
	Eigen::VectorXd phi_terminals(n);
	phi_terminals.topRows(n / 2).array() = phi1;
	phi_terminals.bottomRows(n / 2).array() = phi2;
	return phi_terminals;
}

Eigen::SparseMatrix<double> MaxwellSolver::get_Mnu()
{
	if (Mnu.size() == 0) {
		const int nEdges = primalMeshInfo.nEdges;
		assert(nEdges > 0);

		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		triplets.reserve(nEdges);
		for (int i = 0; i < nEdges; i++) {
			const double fA = dual->getFace(i)->getArea();
			const double eL = primal->getEdge(i)->getLength();
			const double value = fA / eL;
			triplets.emplace_back(T(i, i, value));
		}
		Mnu = Eigen::SparseMatrix<double>(nEdges, nEdges);
		Mnu.setFromTriplets(triplets.begin(), triplets.end());
		Mnu.makeCompressed();
	}
	assert(Mnu.size() > 0);
	assert(Mnu.nonZeros() > 0);
	return Mnu;
}

Eigen::SparseMatrix<double> MaxwellSolver::setupOperatorMeps()
{
	const int nPe = primal->getNumberOfEdges();
	Eigen::SparseMatrix<double> Meps(nPe, nPe);
	Meps.setIdentity();
	for (int i = 0; i < nPe; i++) {
		const double dualFaceArea = dual->getFace(i)->getArea();
		const double primalEdgeLength = primal->getEdge(i)->getLength();
		Meps.coeffRef(i, i) = dualFaceArea / primalEdgeLength;
	}
	return Meps;
}

Eigen::SparseMatrix<double> MaxwellSolver::get_Meps()
{
	if (Meps.size() == 0) {
		const int nFacesInner = primalMeshInfo.nFacesInner;
		assert(nFacesInner > 0);

		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		triplets.reserve(nFacesInner);
		for (int i = 0; i < nFacesInner; i++) {
			const double fA = primal->getFace(i)->getArea();
			const double eL = dual->getEdge(i)->getLength();
			const double value = fA / eL;
			triplets.emplace_back(T(i, i, value));
		}
		Meps = Eigen::SparseMatrix<double>(nFacesInner, nFacesInner);
		Meps.setFromTriplets(triplets.begin(), triplets.end());
		Meps.makeCompressed();
	}
	assert(Meps.size() > 0);
	assert(Meps.nonZeros() > 0);
	return Meps;
}

Eigen::SparseMatrix<double> MaxwellSolver::get_Msigma()
{
	if (Meps.size() == 0) {
		const int n = primalMeshInfo.nEdges;
		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		triplets.reserve(n);
		for (int i = 0; i < n; i++) {
			const double sigma = 1; // TODO <<<<----------------------------------------
			const double fA = dual->getFace(i)->getArea();
			const double eL = primal->getEdge(i)->getLength();
			const double value = fA / eL * sigma;
			triplets.emplace_back(T(i, i, value));
		}
		Meps = Eigen::SparseMatrix<double>(n, n);
		Meps.setFromTriplets(triplets.begin(), triplets.end());
		Meps.makeCompressed();
	}
	return Meps;
}

Eigen::SparseMatrix<double> MaxwellSolver::get_bndryInclOp()
{
	if (bndryInclOp.size() == 0) {
		const int nVertices = primalMeshInfo.nVertices;
		const int nVerticesBoundary = primalMeshInfo.nVerticesBoundary;
		assert(nVerticesBoundary > 0);
		assert(nVertices > nVerticesBoundary);

		typedef Eigen::Triplet<double> TripletType;
		std::vector<TripletType> triplets;
		triplets.reserve(nVerticesBoundary);

		for (int i = 0; i < nVerticesBoundary; i++) {
			const int offset = nVertices - nVerticesBoundary;
			triplets.emplace_back(TripletType(offset + i, i, 1.0));
		}
		bndryInclOp = Eigen::SparseMatrix<double>(nVertices, nVerticesBoundary);
		bndryInclOp.setFromTriplets(triplets.begin(), triplets.end());
		bndryInclOp.makeCompressed();
	}
	return bndryInclOp;
}

void MaxwellSolver::init()
{
	assert(primal->getNumberOfCells() > 0);
	assert(dual->getNumberOfCells() > 0);
	primalMeshInfo = primal->getMeshInfo();
	dualMeshInfo = dual->getMeshInfo();


	// Define data vectors
	B_h = Eigen::VectorXd::Zero(primal->getNumberOfFaces());
	E_h = Eigen::VectorXd::Zero(primal->getNumberOfEdges());
	H_h = Eigen::VectorXd::Zero(dual->getNumberOfEdges());
	J_h = Eigen::VectorXd::Zero(dual->getNumberOfFaces());
}

Eigen::SparseMatrix<int> MaxwellSolver::setupOperatorQ()
{
	const Eigen::VectorXi vertexTypes = primal->getVertexTypes();
	const Eigen::VectorXi edgeTypes = primal->getEdgeTypes();

	const int boundaryVertexType = static_cast<int>(Vertex::Type::Boundary);
	const int terminalVertexType = static_cast<int>(Vertex::Type::Terminal);
	const int innerVertexType = static_cast<int>(Vertex::Type::Inner);

	const int boundaryEdgeType = static_cast<int>(Edge::Type::Boundary);
	const int interiorToBoundaryEdgeType = static_cast<int>(Edge::Type::InteriorToBoundary);
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);

	const int nPv = primal->getNumberOfVertices();
	const int nPvb = (vertexTypes.array() == boundaryVertexType).count() + (vertexTypes.array() == terminalVertexType).count();
	const int nPvi = (vertexTypes.array() == innerVertexType).count();
	const int nPe = primal->getNumberOfEdges();
	const int nPei = (edgeTypes.array() == interiorEdgeType).count() + (edgeTypes.array() == interiorToBoundaryEdgeType).count();
	const int nPeb = (edgeTypes.array() == boundaryEdgeType).count();

	assert(nPe == nPei + nPeb);

	Eigen::SparseMatrix<int> X(nPv, nPvb);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nPvb; i++) {
		triplets.push_back(T(nPvi + i, i, 1));
	}
	X.setFromTriplets(triplets.begin(), triplets.end());
	X.makeCompressed();

	const Eigen::SparseMatrix<int> G = primal->get_e2vMap();
	Eigen::SparseMatrix<int> GX = G * X;

	assert(nPei < nPe);
	Eigen::SparseMatrix<int> id(nPe, nPei);
	//id.setIdentity(); // only for square matrices
	for (int i = 0; i < nPei; i++) {
		id.coeffRef(i, i) = 1;
	}

	// Number of degrees of freedom
	const int nDof = id.cols() + GX.cols();
	assert(nDof == id.cols() + GX.cols());
	Eigen::SparseMatrix<int> Q = Eigen::SparseMatrix<int>(nPe, nDof);
	Q.leftCols(id.cols()) = id;
	Q.rightCols(GX.cols()) = GX;
	assert(Q.rows() == nPe);
	assert(Q.cols() == nDof);
	return Q;
}

Eigen::SparseMatrix<double> MaxwellSolver::setupOperatorMnu()
{
	const int nPf = primal->getNumberOfFaces();
	Eigen::SparseMatrix<double> Mnu(nPf, nPf);
	Mnu.setIdentity();
	for (int i = 0; i < nPf; i++) {
		const double dualEdgeLength = dual->getEdge(i)->getLength();
		const double primalFaceArea = primal->getFace(i)->getArea();
		Mnu.coeffRef(i, i) = dualEdgeLength / primalFaceArea;
	}
	return Mnu;
}
