#include "AppmSolver.h"

AppmSolver::AppmSolver() 
	: AppmSolver(PrimalMesh::PrimalMeshParams())
{
}

AppmSolver::AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams)
{
	readParameters("AppmSolverParams.txt");
	init_meshes(primalMeshParams);  // Initialize primal and dual meshes

	interpolator = new Interpolator(primalMesh, dualMesh);
	std::cout << "- Interpolator ready" << std::endl;

	twofluidSolver = new TwoFluidSolver(primalMesh, dualMesh, interpolator);
	std::cout << "- TwoFluidSolver ready" << std::endl;

	maxwellSolver  = new MaxwellSolver (primalMesh, dualMesh, interpolator);
	maxwellSolver->parameters.lambdaSquare = lambda*lambda;
	std::cout << "- MaxwellSolver ready" << std::endl;

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

	twofluidSolver->applyInitialCondition();
	maxwellSolver->applyInitialCondition();
	writeSnapshot(iteration, time);
	while (time < maxTime && iteration < maxIterations) {
		std::cout << "Iteration " << iteration << ",\t time = " << time << std::endl;
		
		double dt = twofluidSolver->updateFluxesExplicit();  // Compute time step
		double max_B = maxwellSolver->getInterpolated_B().rowwise().norm().maxCoeff();
		dt /= (1e4 * max_B * dt + 1);
		if (time + dt > maxTime) dt = maxTime - time;

		twofluidSolver->updateRateOfChange(false);                 // Compute temporary quantities for later calculations
		maxwellSolver->solveLinearSystem(time,                     // Solve the linear system and update <e> vector
										 dt, 
										 twofluidSolver->get_M_sigma(dt), 
										 twofluidSolver->get_j_aux(dt, maxwellSolver->getInterpolated_B()));
		twofluidSolver->updateMassFluxesImplicit(dt, maxwellSolver->getInterpolated_E());  // Update the flux
		twofluidSolver->timeStepping(dt, maxwellSolver->getInterpolated_E(), maxwellSolver->getInterpolated_B()); // Evolve the fluid variables
		maxwellSolver->evolveMagneticFlux(dt);  // Evolve <b> vector
		verboseDiagnosis();
		
		iteration++;
		time += dt;
		if (iteration % itersPerWrite == 0)  writeSnapshot(iteration, time);
	}
	if (timeStamps.back().first != iteration)  writeSnapshot(iteration, time);
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	writeSolutionPrimalVertex();// boundary electric potential
	writeSolutionDualCell();    // number density, velocity, pressure of all species; E-field; B-field
	writeSolutionDualFace();	// electric current
	writeSolutionPrimalEdge();	// edge voltage
	writeSolutionPrimalFace();  // magnetic flux

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

void AppmSolver::writeSolutionPrimalVertex() {
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeStamps = timeStamps.size();
	for (int i = 0; i < nTimeStamps; i++) {
		XdmfTime time(timeStamps[i].second);
		XdmfGrid snapshotGrid = getSnapshotPrimalVertex(timeStamps[i].first);
		snapshotGrid.addChild(time);
		time_grid.addChild(snapshotGrid);
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_primal_edge.xdmf");
	file << root;
	file.close();
}


void AppmSolver::writeSolutionPrimalEdge()
{	
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeStamps = timeStamps.size();
	for (int i = 0; i < nTimeStamps; i++) {
		XdmfTime time(timeStamps[i].second);
		XdmfGrid snapshotGrid = getSnapshotPrimalEdge(timeStamps[i].first);
		snapshotGrid.addChild(time);
		time_grid.addChild(snapshotGrid);
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_primal_edge.xdmf");
	file << root;
	file.close();
}

void AppmSolver::writeSolutionPrimalFace()
{	
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeStamps = timeStamps.size();
	for (int i = 0; i < nTimeStamps; i++) {
		XdmfTime time(timeStamps[i].second);
		XdmfGrid snapshotGrid = getSnapshotPrimalFace(timeStamps[i].first);
		snapshotGrid.addChild(time);
		time_grid.addChild(snapshotGrid);
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_primal_face.xdmf");
	file << root;
	file.close();
}

void AppmSolver::writeSolutionDualCell() {
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeStamps = timeStamps.size();
	for (int i = 0; i < nTimeStamps; i++) {
		XdmfTime time(timeStamps[i].second);
		XdmfGrid snapshotGrid = getSnapshotDualCell(timeStamps[i].first);
		snapshotGrid.addChild(time);
		time_grid.addChild(snapshotGrid);
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_dual_cell.xdmf");
	file << root;
	file.close();
}

void AppmSolver::writeSolutionDualFace() {
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid time_grid(XdmfGrid::Tags("Time Grid", XdmfGrid::GridType::Collection, XdmfGrid::CollectionType::Temporal));

	const int nTimeStamps = timeStamps.size();
	for (int i = 0; i < nTimeStamps; i++) {
		XdmfTime time(timeStamps[i].second);
		XdmfGrid snapshotGrid = getSnapshotDualFace(timeStamps[i].first);
		snapshotGrid.addChild(time);
		time_grid.addChild(snapshotGrid);
	}
	domain.addChild(time_grid);
	root.addChild(domain);
	std::ofstream file("solutions_dual_face.xdmf");
	file << root;
	file.close();
}

void AppmSolver::writeSnapshot(const int iteration, const double time)
{
	timeStamps.push_back({iteration, time});
	std::cout << "Write snapshot at iteration " << iteration << ", time = " << time << std::endl;
	
	std::stringstream ss_filename;
	ss_filename << "snapshot-" << iteration << ".h5";
	const std::string filename = ss_filename.str();

	H5Writer h5writer(filename);

	twofluidSolver->writeSnapshot(h5writer);
	maxwellSolver->writeSnapshot(h5writer);

}

XdmfGrid AppmSolver::getSnapshotPrimalVertex(const int iteration) {
	XdmfGrid grid = primalMesh->getXdmfVertexGrid();

	// Attribute: Electric potential phi
	std::stringstream ss;
	ss << "snapshot-" << iteration << ".h5:/phi";
	XdmfAttribute potential(
		XdmfAttribute::Tags("electric potential", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	potential.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh->getNumberOfVertices()},
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			ss.str()
		)
	);
	grid.addChild(potential);

	return grid;
}

// TODO
XdmfGrid AppmSolver::getSnapshotPrimalEdge(const int iteration)
{	
	XdmfGrid grid = primalMesh->getXdmfEdgeGrid();

	// Attribute: Electric field E
	std::stringstream ss;
	ss << "snapshot-" << iteration << ".h5:/e";
	XdmfAttribute edge_voltage(
		XdmfAttribute::Tags("edge voltage", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	edge_voltage.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh->getNumberOfEdges()},
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			ss.str()
		)
	);
	grid.addChild(edge_voltage);

	return grid;
}

// TODO
XdmfGrid AppmSolver::getSnapshotPrimalFace(const int iteration)
{
	XdmfGrid grid = primalMesh->getXdmfFaceGrid();
	
	// Attribute: Magnetic flux B
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/b";
		XdmfAttribute mag_flux(
			XdmfAttribute::Tags("Magnetic Flux", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		mag_flux.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ primalMesh->getNumberOfFaces() },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(mag_flux);
	}

	return grid;
}

// TODO
XdmfGrid AppmSolver::getSnapshotDualEdge(const int iteration)
{
	XdmfGrid grid = dualMesh->getXdmfEdgeGrid();

	// Attribute: 
	{

	}
	return grid;
}

// TODO
XdmfGrid AppmSolver::getSnapshotDualFace(const int iteration)
{
	XdmfGrid grid = dualMesh->getXdmfFaceGrid();

	// Attribute: Electric current J
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5" << ":/j";
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Electric Current", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfFaces()},
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
		ss << "snapshot-" << iteration << ".h5:/electron-density";
		XdmfAttribute density(
			XdmfAttribute::Tags("electron number density", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
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
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/ion-density";
		XdmfAttribute density(
			XdmfAttribute::Tags("ion number density", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
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
		ss << "snapshot-" << iteration << ".h5:/electron-pressure";
		XdmfAttribute pressure(
			XdmfAttribute::Tags("electron pressure", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
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
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/ion-pressure";
		XdmfAttribute pressure(
			XdmfAttribute::Tags("ion pressure", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
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
		ss << "snapshot-" << iteration << ".h5:/electron-velocity";
		XdmfAttribute velocity(
			XdmfAttribute::Tags("electron velocity", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
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
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/ion-velocity";
		XdmfAttribute velocity(
			XdmfAttribute::Tags("ion velocity", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
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
		// Attribute: E-field
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/E";
		XdmfAttribute efield(
			XdmfAttribute::Tags("E-field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		efield.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfCells(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(efield);
	}
	// Attribute: B-field
	{
		std::stringstream ss;
		ss << "snapshot-" << iteration << ".h5:/B";
		XdmfAttribute bfield(
			XdmfAttribute::Tags("B-field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		bfield.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh->getNumberOfCells(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				ss.str()
			));
		grid.addChild(bfield);
	}
	return grid;
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
		if (tag == "lambda") {
			std::istringstream(line.substr(pos + 1)) >> this->lambda;
		}
		if (tag == "itersPerWrite") {
			std::istringstream(line.substr(pos + 1)) >> this->itersPerWrite;
		}
	}

	std::cout << std::endl;
	std::cout << "Appm Solver parameters:" << std::endl;
	std::cout << "======================="  << std::endl;
	std::cout << "maxIterations:  " << maxIterations << std::endl;
	std::cout << "maxTime:        " << maxTime << std::endl;
	std::cout << "lambda: " << lambda << std::endl;
	std::cout << "=======================" << std::endl;
}

void AppmSolver::verboseDiagnosis() const {
	
	double max_E = 0, max_B = 0;
	Eigen::MatrixXd E = maxwellSolver->getInterpolated_E();
	Eigen::MatrixXd B = maxwellSolver->getInterpolated_B();
	for (int i = 0; i < dualMesh->getNumberOfCells(); i++) {
		if (dualMesh->getCell(i)->getFluidType() == Cell::FluidType::Fluid) {
			max_E = E.row(i).norm() > max_E ? E.row(i).norm() : max_E;
			max_B = B.row(i).norm() > max_B ? B.row(i).norm() : max_B;
		}
	}
	std::cout << " ---------- max e component :" << maxwellSolver->e.array().abs().maxCoeff() << std::endl;
	std::cout << " ---------- max b component :" << maxwellSolver->b.array().abs().maxCoeff() << std::endl;
	std::cout << " ---------- max j component :" << maxwellSolver->j.array().abs().maxCoeff() << std::endl;
	std::cout << " ---------- max E-field :" << E.rowwise().norm().maxCoeff() << std::endl;
	std::cout << " ---------- max E-field in fluid:" << max_E << std::endl;
	std::cout << " ---------- max B-field :" << B.rowwise().norm().maxCoeff() << std::endl;
	std::cout << " ---------- max B-field in fluid:" << max_B << std::endl;
	std::cout << " ---------- max h component :" << maxwellSolver->hp.array().abs().maxCoeff() << std::endl;
	std::cout << " ---------- max d component :" << maxwellSolver->dp.array().abs().maxCoeff() << std::endl;
	std::cout << " ---------- max electron vel :" << twofluidSolver->electron_solver.U.middleCols(1,3).rowwise().norm().maxCoeff() << std::endl;
	std::cout << " ---------- max ion vel :" << twofluidSolver->ion_solver.U.middleCols(1,3).rowwise().norm().maxCoeff() << std::endl;
	std::cout << " ---------- electron density : [" << twofluidSolver->electron_solver.U.col(0).maxCoeff() << ", "
													<< twofluidSolver->electron_solver.U.col(0).minCoeff() << "]" << std::endl;
	std::cout << " ---------- ion density : [" << twofluidSolver->ion_solver.U.col(0).maxCoeff() << ", "
											   << twofluidSolver->ion_solver.U.col(0).minCoeff() << "]" << std::endl;

}
