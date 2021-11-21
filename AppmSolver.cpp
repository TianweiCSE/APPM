#include "AppmSolver.h"

AppmSolver::AppmSolver() 
	: AppmSolver(PrimalMesh::PrimalMeshParams())
{
}

AppmSolver::AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams)
{
	// readParameters("AppmSolverParams.txt");
	init_meshes(primalMeshParams);  // Initialize primal and dual meshes
	interpolator = new Interpolator(primalMesh, dualMesh);

	// MaxwellSolver::MaxwellParams maxwellParams;
	// maxwellParams.lambdaSquare = this->lambdaSquare;
	// maxwellSolver = new MaxwellSolverImplicitEuler(primalMesh, dualMesh, maxwellParams);
	// fluidSolver = new SingleFluidSolver(&dualMesh);
	// fluidSolver = new TwoFluidSolver(&dualMesh);
	
	twofluidSolver = new TwoFluidSolver(primalMesh, dualMesh, interpolator);
	maxwellSolver  = new MaxwellSolver (primalMesh, dualMesh, interpolator);

	// B_vertex = Eigen::Matrix3Xd::Zero(3, primalMesh->getNumberOfVertices());
	// init_RaviartThomasInterpolation();
}

AppmSolver::~AppmSolver()
{
	if (twofluidSolver != nullptr) {
		delete twofluidSolver;
		twofluidSolver = nullptr;
	}

	if (maxwellSolver != nullptr) {
		delete maxwellSolver;
		maxwellSolver = nullptr;
	}
}

void AppmSolver::run()
{
	double time = 0, dt;
	int iteration = 0;

	// For testing of RT interpolation:
	//setAzimuthalMagneticFluxField();
	//setUniformMagneticFluxField(Eigen::Vector3d(1,1,1));
	//interpolateMagneticFluxToPrimalVertices();

	// initialize current flow
	/*
	this->isMaxwellCurrentSource = false;
	maxwellSolver->isMaxwellCurrentSource = isMaxwellCurrentSource;
	if (isMaxwellCurrentSource) {
		const double x1 = -0.5;
		const double x2 = 0.5;
		const double z1 = 0.24;
		const double z2 = 0.76;
		maxwellSolver->setTorusCurrent(x1, x2, z1, z2);
	}*/
	twofluidSolver->applyInitialCondition();
	writeSnapshot(iteration, time);
	while (time < maxTime && iteration < maxIterations) {
		std::cout << "Iteration " << iteration << ",\t time = " << time << std::endl;
		
		const double dt = twofluidSolver->updateFluxesExplicit();  // Compute time step
		twofluidSolver->updateRateOfChange();                      // Compute temporary quantities for later calculations
		maxwellSolver->solveLinearSystem(dt, 
										 twofluidSolver->get_M_sigma(dt, maxwellSolver->getInterpolated_B()), 
										 twofluidSolver->get_j_aux(dt));
		twofluidSolver->updateFluxesImplicit(maxwellSolver->getInterpolated_E());
		twofluidSolver->timeStepping(dt, maxwellSolver->getInterpolated_E(), maxwellSolver->getInterpolated_B());
		maxwellSolver->timeStepping(dt);
	
		// Maxwell equations
		/*
		if (isMaxwellEnabled) {
			maxwellSolver->updateMaxwellState(dt, time);
			interpolateMagneticFluxToPrimalVertices();
		}*/

		iteration++;
		time += dt;
		writeSnapshot(iteration, time);
	}
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	writeSolutionDualCell();

	// test_raviartThomas();

	// Use Paraview (version 5.6.0) to visualize.
	// Light data is given in XDMF files (version 3). 
	// Heavy data is stored in HDF5 files.

	// Note that polyhedrons have been added to XDMF in version 3; version 2 does not support polyhedrons.
	// Moreover, Paraview has 3 readers: XDMF Reader, Xdmf3ReaderS, Xdmf3ReaderT.
	// Xdmf3Reader3: can     read polyhedrons, but not grid-of-grids (i.e., grids with GridType=Tree).
	// XDMF Reader:  can not read polyhedrons, but     grid-of-grids.
	// Therefore, the mesh data for vertices, edges, and faces, are separated from the volume data.
	// writeXdmf();
	// writeXdmfDualVolume();
}

void AppmSolver::interpolateMagneticFluxToPrimalVertices()
{
	B_vertex.setZero();
	// For each cell ...
	//const int nCells = primalMesh->getNumberOfCells();
	// For each cell that has a Piola map
	const int nCells = rt_piolaMatrix.size();  /// the number of triangle prisms
	Eigen::VectorXi countVertexVisits = Eigen::VectorXi::Zero(primalMesh->getNumberOfVertices());

	const Eigen::VectorXd B_h = maxwellSolver->getBstate();

	for (int cidx = 0; cidx < nCells; cidx++) {
		const Cell * cell = primalMesh->getCell(cidx);
		std::vector<Face*> cellFaces = cell->getFaceList();
		std::vector<Vertex*> bottomVertices = cellFaces[3]->getVertexList();
		std::vector<Vertex*> topVertices = cellFaces[4]->getVertexList();

		// Piola map
		const Eigen::Matrix3d & BK = rt_piolaMatrix[cidx];
		const Eigen::Vector3d & bK = rt_piolaVector[cidx];
		const double detBK = BK.determinant();

		// Get prism vertices
		std::vector<Vertex*> cellVertices(6);
		for (int i = 0; i < 3; i++) {
			cellVertices[i] = bottomVertices[i];
			cellVertices[i + 3] = topVertices[i];
		}

		//Eigen::Matrix3cd vertexCoords(3, 6);
		Eigen::MatrixXd refCoords(3, 6);
		Eigen::VectorXi vertexIdx(6);
		for (int i = 0; i < 6; i++) {
			Vertex * v = cellVertices[i];
			const Eigen::Vector3d pos = v->getPosition();
			//vertexCoords.col(i) = pos;
			vertexIdx(i) = v->getIndex();
			refCoords.col(i) = BK.inverse() * (pos - bK);
		}
		const double tolerance = 16 * std::numeric_limits<double>::epsilon();
		// Check that reference coordinates are in [0,1]
		assert((refCoords.array() >= -tolerance).all());
		assert(((refCoords.array() - 1) <= tolerance).all());

		// Get coefficients for RT interpolation
		Eigen::VectorXd coeff(5);
		for (int idx = 0; idx < 5; idx++) {
			Face * face = cellFaces[idx];
			const int faceIncidence = cell->getOrientation(face);
			const int fidx = face->getIndex();
			const int factor = 1;// face->isBoundary() ? 2 : 1;
			coeff(idx) = faceIncidence * B_h(fidx) * factor;
		}

		// Interpolate B-field with RT-basis functions in actual space
		Eigen::Matrix3Xd B_vertex_local = Eigen::Matrix3Xd::Zero(3, 6);
		for (int idx = 0; idx < 5; idx++) {
			for (int i = 0; i < 6; i++) {
				const int vIdx = vertexIdx(i);
				const Eigen::Vector3d refPos = refCoords.col(i);
				const Eigen::Vector3d rt = Numerics::raviartThomasBasis(idx, refPos);
				const int directionFactor = (idx >= 3) ? 2 : 1; // Because: faces in z-direction (i.e., face idx < 3) are quads with normal vector perpendicular to z-axis, and face idx >= 3 are triangles. Hence, det(B) acts differently.
				B_vertex_local.col(i) += directionFactor * coeff(idx) / detBK * BK * rt;
			}
		}
		for (int i = 0; i < 6; i++) {
			const int vIdx = vertexIdx(i);
			countVertexVisits(vIdx) += 1;
			B_vertex.col(vIdx) += B_vertex_local.col(i);
		}
	}
	for (int i = 0; i < primalMesh->getNumberOfVertices(); i++) {
		if (countVertexVisits(i) > 0) {
			B_vertex.col(i) /= countVertexVisits(i);
		}
	}
}

//void AppmSolver::test_raviartThomas()
//{
//	std::cout << "Test Raviart Thomas Interpolation" << std::endl;
//	// setAzimuthalMagneticFluxField();
//
//	B_vertex.setZero();
//
//
//
//	// Primal cells are prisms.
//	// For each cell, 
//	// - compute Piola map, 
//	// - compute reference coordinates of vertices
//	// - evaluate RT basis functions 
//	// - map RT basis functions from reference to actual coordinates
//	// - multiply by magnetic flux B
//	// - write data to file: B at vertices
//	const int nCells = primalMeshInfo.nCells;
//	for (int cidx = 0; cidx < nCells; cidx++) {
//#ifdef _RT_ONECELL
//		std::cout << "cidx = " << cidx << std::endl;
//#endif
//		const Cell* cell = primalMesh->getCell(cidx);
//		const std::vector<Face*> cellFaces = cell->getFaceList();
//		assert(cellFaces.size() == 5); // prism cells
//		assert(cellFaces[3]->getVertexList().size() == 3); // triangle faces are at end of list
//		assert(cellFaces[4]->getVertexList().size() == 3);
//
//		const Face * bottomFace = cellFaces[3];
//		const Face * topFace = cellFaces[4];
//		const std::vector<Vertex*> bottomVertices = bottomFace->getVertexList();
//		const std::vector<Vertex*> topVertices = topFace->getVertexList();
//
//		// check if vertices of triangle faces have same ordering
//		for (int i = 0; i < 3; i++) {
//			Eigen::Vector3d v = topVertices[i]->getPosition() - bottomVertices[i]->getPosition();
//			v.normalize();
//			assert(v.cross(zUnit).norm() < 10 * std::numeric_limits<double>::epsilon());
//		}
//
//		// check if bottom triangle vertices form a right-handed system
//		const Eigen::Vector3d fc = bottomFace->getCenter();
//		for (int i = 0; i < 3; i++) {
//			const Eigen::Vector3d v0 = bottomVertices[i]->getPosition();
//			const Eigen::Vector3d v1 = bottomVertices[(i + 1) % 3]->getPosition();
//			const Eigen::Vector3d a = v0 - fc;
//			const Eigen::Vector3d b = v1 - v0;
//			const Eigen::Vector3d n = a.normalized().cross(b.normalized());
//			assert(n.dot(Eigen::Vector3d::UnitZ()) > 0);
//		}
//
//		std::vector<Vertex*> cellVertices(6);
//		for (int i = 0; i < 3; i++) {
//			cellVertices[i] = bottomVertices[i];
//			cellVertices[i + 3] = topVertices[i];
//		}
//
//		const Eigen::Vector3d v0 = cellVertices[0]->getPosition();
//		const Eigen::Vector3d v1 = cellVertices[1]->getPosition();
//		const Eigen::Vector3d v2 = cellVertices[2]->getPosition();
//		const Eigen::Vector3d v3 = cellVertices[5]->getPosition();
//
//		// Piola map: xRef -> x := BK*xRef + bk
//		Eigen::Matrix3d BK;
//		BK.col(0) = v0 - v2;
//		BK.col(1) = v1 - v2;
//		BK.col(2) = v3 - v2;
//		const double detBK = BK.determinant();
//		const Eigen::Vector3d bk = v2;
//#ifdef _RT_ONECELL
//		std::cout << "BK: " << BK << std::endl;
//		std::cout << "det(BK): " << detBK << std::endl;
//#endif
//		// vertex coordinates in actual space
//		Eigen::Matrix3Xd vertexCoords(3, 6);
//		for (int i = 0; i < 6; i++) {
//			vertexCoords.col(i) = cellVertices[i]->getPosition();
//		}
//		// points in reference coordinates 
//		const int nSamples = 5;
//		Eigen::Matrix3Xd refCoords3d(3, 6);
//		refCoords3d.setZero();
//#ifdef _RT_ONECELL
//		refCoords3d = getPrismReferenceCoords(nSamples);
//#else
//		// Reference coordinates of cell vertices
//		for (int i = 0; i < 6; i++) {
//			const Eigen::Vector3d pos = vertexCoords.col(i);
//			const Eigen::Vector3d refPos = BK.inverse() * (pos - bk);
//			refCoords3d.col(i) = refPos;
//		}
//#endif
//		const double tolerance = 16 * std::numeric_limits<double>::epsilon();
//		assert((refCoords3d.array() >= -1 * std::numeric_limits<double>::epsilon()).all());
//		if ( ( (refCoords3d.array() - 1) > tolerance ).any() ) {
//			std::cout << "max value: " << (refCoords3d.array() - 1).maxCoeff() << std::endl;
//		}
//		assert((refCoords3d.array() - 1 <= tolerance).all());
//
//		for (int i = 0; i < cellFaces.size(); i++) {
//			const Face * face = cellFaces[i];
//			const Eigen::Vector3d fn = face->getNormal();
//			const Eigen::Vector3d fc = face->getCenter();
//
//#ifdef _RT_ONECELL
//			std::cout << "Face " << i << ": ";
//			std::cout << "fc: " << fc.transpose() << "\t ";
//			std::cout << "fn: " << fn.transpose() << "\t "; 
//			std::cout << "incid: " << cell->getOrientation(face) << "\t ";
//			std::cout << std::endl;
//#endif
//		}
//
//		Eigen::VectorXd coeff = Eigen::VectorXd::Zero(5);
//		Eigen::VectorXd B_h_loc(5);
//		for (int i = 0; i < 5; i++) {
//			const Face * face = cellFaces[i];
//			const int fidx = face->getIndex();
//			const int orientation = cell->getOrientation(face);
//			B_h_loc(i) = B_h(fidx);
//			coeff(i) = orientation * B_h(fidx);
//		}
//#ifdef _RT_ONECELL
//		std::cout << "B_h_loc: " << B_h_loc << std::endl;
//		std::cout << "coeff: " << coeff << std::endl;
//#endif
//		
//		// for each vertex
//		Eigen::Matrix3Xd Btest(3, refCoords3d.cols());
//		Eigen::Matrix3Xd xCoords3d(3, refCoords3d.cols());
//		for (int i = 0; i < refCoords3d.cols(); i++) {
//			const Eigen::Vector3d refPos = refCoords3d.col(i);
//			Eigen::Vector3d Bv;
//			Bv.setZero();
//			// for each basis function
//			for (int idx = 0; idx < 5; idx++) {
//				const Eigen::Vector3d rt = Numerics::raviartThomasBasis(idx, refPos);
//				Bv += coeff(idx) / detBK * BK * rt;
//			}
//#ifdef _RT_ONECELL
//			Btest.col(i) = Bv;
//			xCoords3d.col(i) = BK * refPos + bk;
//#else
//			const int vIdx = cellVertices[i]->getIndex();
//			B_vertex.col(vIdx) += Bv;
//#endif
//		}
//#ifdef _RT_ONECELL
//		H5Writer tempwriter("temp.h5");
//		tempwriter.writeData(refCoords3d, "/refCoords3d");
//		tempwriter.writeData(Btest, "/Btest");
//		tempwriter.writeData(xCoords3d, "/xCoords3d");
//		break;
//#endif
//	}
//
//	H5Writer h5writer("rt.h5");
//	h5writer.writeData(B_vertex, "/Bvertex");
//}

const Eigen::Matrix3Xd AppmSolver::getPrismReferenceCoords(const int nSamples)
{
	Eigen::Matrix3Xd refCoords(3, (nSamples + 2) * (nSamples + 1) / 2);
	int sIdx = 0;
	for (int i = 0; i <= nSamples; i++) {
		if (i == 0) {
			refCoords.col(sIdx++) = Eigen::Vector3d(0, 0, 0);
			continue;
		}
		Eigen::Vector3d a, b;
		a = 1.0*i / nSamples * Eigen::Vector3d(1, 0, 0);
		b = 1.0*i / nSamples * Eigen::Vector3d(0, 1, 0);
		Eigen::Vector3d d = b - a;
		const Eigen::VectorXd samples = Eigen::VectorXd::LinSpaced(i + 1, 0., 1.);
		for (int k = 0; k < samples.size(); k++) {
			const Eigen::Vector3d v = a + samples(k) * (b - a);
			refCoords.col(sIdx++) = v;
		}
	}
	assert(sIdx == refCoords.cols());
	const int n = refCoords.cols();
	Eigen::Matrix3Xd refCoords3d = refCoords.replicate(1, nSamples + 1);
	for (int i = 0; i <= nSamples; i++) {
		const int startCol = i * n;
		refCoords3d.block(2, startCol, 1, n).array() = 1.0 * i / nSamples;
	}
	return refCoords3d;
}

void AppmSolver::init_meshes(const PrimalMesh::PrimalMeshParams & primalParams)
{
	std::cout << "============== Init primal mesh ============" << std::endl;

	primalMesh = new PrimalMesh(primalParams); 
	primalMesh->init();
	primalMesh->check();
	primalMesh->writeToFile();
	primalMesh->writeXdmf();

	std::cout << "=============== Init dual mesh =============" << std::endl;
	dualMesh = new DualMesh(primalMesh);
	dualMesh->init();
	dualMesh->check();
	dualMesh->writeToFile();
	dualMesh->writeXdmf();

	std::cout << "=============== Pimal/Dual mesh generated ===============" << std::endl;
}



void AppmSolver::writeSolutionPrimalEdge()
{	
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimesteps = timeStamps.size();
	for (int i = 0; i < nTimesteps; i++) {
		XdmfTime time(timeStamps[i]);
		time_grid.addChild(time);
		time_grid.addChild(getSnapshotPrimalEdge(i));
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_primal_edge.xdmf");
	file << root;
	file.close();
}
/*
void AppmSolver::writeXdmf() {
	const int nTimesteps = timeStamps.size();
	
	const std::string filename = "appm.xdmf";
	std::string gridPrimalEdges;
	std::string gridPrimalFaces;
	std::string gridDualEdges;
	std::string gridDualFaces;

	std::ofstream file(filename);
	file << "<?xml version = \"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
	file << "<Domain>" << std::endl;
	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
	for (int i = 0; i < nTimesteps; i++) {
		const double time = this->timeStamps[i];

		file << "<Grid Name=\"Grid of Grids\" GridType=\"Tree\">" << std::endl;
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridPrimalEdges(i) << std::endl;
		file << xdmf_GridPrimalFaces(i) << std::endl;
		file << xdmf_GridDualEdges(i) << std::endl;
		file << xdmf_GridDualFaces(i) << std::endl;
		file << "</Grid>" << std::endl;
	}
	file << "</Grid>" << std::endl;
	file << "</Domain>" << std::endl;
	file << "</Xdmf>" << std::endl;
}*/

void AppmSolver::writeSolutionDualCell() {
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeSteps = timeStamps.size();
	for (int i = 0; i < nTimeSteps; i++) {
		XdmfTime time(timeStamps[i]);
		time_grid.addChild(time);
		time_grid.addChild(getSnapshotDualCell(i));
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_dual_cell.xdmf");
	file << root;
	file.close();
}

/*
void AppmSolver::writeXdmfDualVolume()
{
	const int nTimesteps = timeStamps.size();
	const std::string filename = "appm-volume.xdmf";
	std::ofstream file(filename);
	file << "<?xml version = \"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
	file << "<Domain>" << std::endl;
	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
	for (int i = 0; i < nTimesteps; i++) {
		const double time = this->timeStamps[i];
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridDualCells(i) << std::endl;
	}
	file << "</Grid>" << std::endl;
	file << "</Domain>" << std::endl;
	file << "</Xdmf>" << std::endl;
}*/

void AppmSolver::writeSnapshot(const int iteration, const double time)
{
	timeStamps.push_back(time);
	std::cout << "Write output at iteration " << iteration << ", time = " << time << std::endl;
	
	std::stringstream ss_filename;
	ss_filename << "appm-" << iteration << ".h5";
	const std::string filename = ss_filename.str();

	H5Writer h5writer(filename);

	// Fluid states
	twofluidSolver->writeSnapshot(h5writer);

	// Maxwell states
	// maxwellSolver->writeStates(h5writer);

	// Interpolated values of B-field to primal vertices
	// h5writer.writeDoubleMatrix(B_vertex, "/Bvertex");
	/*
	Eigen::VectorXd timeVec(1);
	timeVec(0) = time;
	h5writer.writeDoubleVector(timeVec, "/time");
	Eigen::VectorXi iterVec(1);
	iterVec(0) = iteration;
	h5writer.writeIntVector(iterVec, "/iteration");
	*/
}

XdmfGrid AppmSolver::getSnapshotPrimalEdge(const int iteration)
{	
	XdmfGrid grid = primalMesh->getXdmfEdgeGrid();

	// Attribute: Electric field E
	std::stringstream ss;
	ss << "appm-" << iteration << ".h5:/E";
	XdmfAttribute e_field(
		XdmfAttribute::Tags("Electric field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	);
	e_field.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh->getNumberOfEdges(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			ss.str()
		)
	);
	grid.addChild(e_field);

	return grid;
}

XdmfGrid AppmSolver::getSnapshotPrimalFace(const int iteration)
{
	XdmfGrid grid = primalMesh->getXdmfFaceGrid();
	
	// Attribute: Magnetic flux B
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5:/B";
		XdmfAttribute BfieldAttribute(
			XdmfAttribute::Tags("Magnetic Flux", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		BfieldAttribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ primalMesh->getNumberOfFaces(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(BfieldAttribute);
	}
	
	// Attribute: Magnetic flux B
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5:/Bvertex";
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Magnetic Flux Interpolated", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Node)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ primalMesh->getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(attribute);
	}

	return grid;
}

XdmfGrid AppmSolver::getSnapshotDualEdge(const int iteration)
{
	XdmfGrid grid(XdmfGrid::Tags("Dual Edges"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Polyline, dualMesh->getNumberOfEdges(), 2)
	);
	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/edge2vertex";
		topology.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ 2 * dualMesh->getNumberOfEdges() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(topology);
	}

	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/vertexPos";
		XdmfGeometry geometry;
		geometry.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str())
		);
		grid.addChild(geometry);
	}

	// Attribute: Edge index
	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/edgeIdx";
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Edge index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfEdges() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(attribute);
	}

	// Attribute: Magnetic Field H
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5" << ":/H";
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Magnetic field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfEdges(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(attribute);
	}
	return grid;
}

XdmfGrid AppmSolver::getSnapshotDualFace(const int iteration)
{
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	XdmfGrid grid(XdmfGrid::Tags("Dual Faces"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, dualMesh->getNumberOfFaces())
	);
	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/face2vertex";
		topology.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ nElements },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(topology);
	}

	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/vertexPos";
		XdmfGeometry geometry;
		geometry.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str())
		);
		grid.addChild(geometry);
	}
	
	// Attribute: Face index
	{
		std::stringstream ss;
		ss << dualMesh->getPrefix() << "-mesh.h5:/faceIndex";
		XdmfAttribute faceIndexAttribute(
			XdmfAttribute::Tags("Face index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		faceIndexAttribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfFaces() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(faceIndexAttribute);
	}
	

	//// Attribute: Displacement Field D
	//{
	//	XdmfAttribute attribute(
	//		XdmfAttribute::Tags("Displacement Field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	//	);
	//	attribute.addChild(
	//		XdmfDataItem(XdmfDataItem::Tags(
	//			{ dualMesh->getNumberOfFaces(), 3 },
	//			XdmfDataItem::NumberType::Float,
	//			XdmfDataItem::Format::HDF),
	//			(std::stringstream() << dataFilename << ":/D").str()
	//		));
	//	grid.addChild(attribute);
	//}

	// Attribute: Electric current J
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5" << ":/J";
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Electric Current", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfFaces(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(attribute);
	}

	return grid;
}

XdmfGrid AppmSolver::getSnapshotDualCell(const int iteration)
{
	XdmfGrid grid = dualMesh->getXdmfCellGrid();
	
	// Attribute: B-field at primal vertices = dual cell centers

	// Attribute: Density
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5:/density";
		XdmfAttribute density(
			XdmfAttribute::Tags("Number density", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		density.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfCells() },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		grid.addChild(density);
	}
	
	// Attribute: Pressure
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5:/pressure";
		XdmfAttribute pressure(
			XdmfAttribute::Tags("Pressure", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		pressure.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfCells() },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			)
		);
		grid.addChild(pressure);
	}
	
	// Attribute: velocity
	{
		std::stringstream ss;
		ss << "appm-" << iteration << ".h5:/velocity";
		XdmfAttribute velocity(
			XdmfAttribute::Tags("Velocity", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		velocity.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfCells(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(velocity);
	}
	return grid;
}

void AppmSolver::init_RaviartThomasInterpolation()
{
	const int nCells = primalMesh->getNumberOfCells();
	rt_piolaMatrix.reserve(nCells);
	rt_piolaVector.reserve(nCells);

	Eigen::VectorXi isPiolaMapDefined(nCells);
	isPiolaMapDefined.setZero();

	bool isInfoPrinted = true;

	for (int cidx = 0; cidx < nCells; cidx++) {
		const Cell * cell = primalMesh->getCell(cidx);
		const std::vector<Face*> cellFaces = cell->getFaceList();
		if (cellFaces.size() == 5) {
			isPiolaMapDefined(cidx) = 1;
		}
		else {
			isPiolaMapDefined(cidx) = 0;
			if (isInfoPrinted) {
				std::cout << "Raviart-Thomas interpolation is only implemented for triangular prisms!" << std::endl;
				isInfoPrinted = false; // show info only once
			}
			continue;
		}
		assert(cellFaces.size() == 5); // prism cells

		// get top and bottom face of prism
		const Face * bottomFace = cellFaces[3];
		const Face * topFace = cellFaces[4];
		assert(bottomFace->getVertexList().size() == 3); // check if faces are triangular
		assert(topFace->getVertexList().size() == 3);

		// get vertices of triangle faces
		const std::vector<Vertex*> bottomVertices = bottomFace->getVertexList();
		const std::vector<Vertex*> topVertices = topFace->getVertexList();

		// check if vertices of triangle faces have same ordering   /// This is ensured by the way it is constructed.
		for (int i = 0; i < 3; i++) {
			const Eigen::Vector3d v = topVertices[i]->getPosition() - bottomVertices[i]->getPosition();
			assert(v.normalized().cross(Eigen::Vector3d::UnitZ()).norm() < 16 * std::numeric_limits<double>::epsilon());
		}

		// check if bottom triangle vertices form a right-handed system
		const Eigen::Vector3d fc = bottomFace->getCenter();
		for (int i = 0; i < 3; i++) {
			const Eigen::Vector3d v0 = bottomVertices[i]->getPosition();
			const Eigen::Vector3d v1 = bottomVertices[(i + 1) % 3]->getPosition();
			const Eigen::Vector3d a = v0 - fc;
			const Eigen::Vector3d b = v1 - v0;
			const Eigen::Vector3d n = a.normalized().cross(b.normalized());
			assert(n.dot(Eigen::Vector3d::UnitZ()) > 0);
		}

		// Collect prism vertices in standard topology format
		std::vector<Vertex*> cellVertices(6);
		for (int i = 0; i < 3; i++) {
			const int offset = 3;
			cellVertices[i] = bottomVertices[i];
			cellVertices[i + offset] = topVertices[i];
		}

		// Vertex positions of prism
		const Eigen::Vector3d A = cellVertices[0]->getPosition();
		const Eigen::Vector3d B = cellVertices[1]->getPosition();
		const Eigen::Vector3d C = cellVertices[2]->getPosition();
		const Eigen::Vector3d F = cellVertices[5]->getPosition();

		// Define Piola map: xRef  -->  BK * xRef + bK
		const Eigen::Vector3d bK = C;
		Eigen::Matrix3d BK;
		BK.col(0) = A - C;
		BK.col(1) = B - C;
		BK.col(2) = F - C;
		rt_piolaMatrix.emplace_back(BK);
		rt_piolaVector.emplace_back(bK);
	}

	// Check if all Piola maps are defined at beginning of list
	const int nPiolaMapsDefined = isPiolaMapDefined.count();
	//std::cout << "nPiolaMapsDefined: " << nPiolaMapsDefined << std::endl;
	//std::ofstream file("isPiolaMapDefined.dat");
	//file << isPiolaMapDefined << std::endl;
	assert(isPiolaMapDefined.topRows(nPiolaMapsDefined).all());
	assert(isPiolaMapDefined.bottomRows(nCells - nPiolaMapsDefined).any() == 0);

	// Truncate list of Piola maps
	rt_piolaMatrix.resize(nPiolaMapsDefined);
	rt_piolaVector.resize(nPiolaMapsDefined);
}

void AppmSolver::readParameters(const std::string & filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File not opened: " << filename;
		exit(-1);
	}

	std::string line;
	const char delim = ':';

	while (std::getline(file, line)) {
		int pos = line.find(delim);
		std::string tag = line.substr(0, pos);

		if (tag == "maxIterations") {
			std::istringstream(line.substr(pos + 1)) >> this->maxIterations;
		}
		if (tag == "maxTime") {
			std::istringstream(line.substr(pos + 1)) >> this->maxTime;
		}
		if (tag == "isFluidEnabled") {
			std::istringstream(line.substr(pos + 1)) >> this->isFluidEnabled;
		}
		if (tag == "isMaxwellEnabled") {
			std::istringstream(line.substr(pos + 1)) >> this->isMaxwellEnabled;
		}
		if (tag == "timestepSize") {
			std::istringstream(line.substr(pos + 1)) >> this->timestepSize;
		}
		if (tag == "lambdaSquare") {
			std::istringstream(line.substr(pos + 1)) >> this->lambdaSquare;
		}
	}

	std::cout << std::endl;
	std::cout << "Appm Solver parameters:" << std::endl;
	std::cout << "======================="  << std::endl;
	std::cout << "maxIterations:  " << maxIterations << std::endl;
	std::cout << "maxTime:        " << maxTime << std::endl;
	std::cout << "isFluidEnabled: " << isFluidEnabled << std::endl;
	std::cout << "isMaxwellEnabled: " << isMaxwellEnabled << std::endl;
	std::cout << "timestepSize: " << timestepSize << std::endl;
	std::cout << "lambdaSquare: " << lambdaSquare << std::endl;
	std::cout << "=======================" << std::endl;
}

const std::string AppmSolver::xdmf_GridPrimalEdges(const int iteration) const
{
	std::stringstream ss;
	ss << "<Grid Name=\"Primal Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << primalMesh->getNumberOfEdges() << "\""
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * primalMesh->getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfVertices() << " 3\"" 
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Electric field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfEdges() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/E" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	ss << "</Grid>";
	return ss.str();
}

const std::string AppmSolver::xdmf_GridPrimalFaces(const int iteration) const
{
	H5Reader h5reader;
	h5reader = H5Reader("primal-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Primal Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << primalMesh->getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic flux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/B" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic flux interpolated\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/Bvertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	ss << "</Grid>";
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualEdges(const int iteration) const
{
	std::stringstream ss;
	ss << "<Grid Name=\"Dual Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << dualMesh->getNumberOfEdges() << "\"" 
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * dualMesh->getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfEdges() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/H" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualFaces(const int iteration) const
{
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh->getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Electric Current\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/J" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualCells(const int iteration) const
{
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
	const int nElements = h5reader.readDataSize("/cell2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Cells\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh->getNumberOfCells() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/cell2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Cell index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/cellIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Cell Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dualMeshTypes.h5:/cellFluidTypes" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;


	ss << "<Attribute Name=\"Magnetic Flux Interpolated\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh->getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/Bvertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	// ss << fluidSolver->getXdmfOutput(iteration);
	ss << "</Grid>";
	return ss.str();
}
