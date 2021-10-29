#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "H5Writer.h"

#include <Eigen/Dense>

class MaxwellSolver
{
public:
	struct MaxwellParams {
		double lambdaSquare = 1.0;
	} parameters;

	MaxwellSolver();
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual);
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual, MaxwellParams & param);
	~MaxwellSolver();

	bool isMaxwellCurrentSource = false;


	virtual void updateMaxwellState(const double dt, const double time) = 0;

	void writeStates(H5Writer & writer) const;

	const Eigen::VectorXd & getBstate() const;

	void setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector);
	void setAzimuthalMagneticFluxField();
	void setTorusCurrent(const double x1, const double x2, const double z1, const double z2);

	Eigen::VectorXd electricPotentialTerminals(const double time);


protected:
	const PrimalMesh * primal = nullptr;
	const DualMesh * dual = nullptr;

	Mesh::MeshInfo primalMeshInfo;
	Mesh::MeshInfo dualMeshInfo;

	Eigen::VectorXd maxwellState;

	Eigen::VectorXd B_h, E_h, H_h, J_h;

	/// return the square matrix of size #primalEdge.      d = Mnu * e      (Yet in the paper d = Meps * e)
	Eigen::SparseMatrix<double> get_Mnu();
	/// return the square matrix of size #primalFaceInner. h = Meps^-1 * b  (Yet in the paper h = Mnu * e )
	Eigen::SparseMatrix<double> get_Meps();
	/// NOT implemented yet
	Eigen::SparseMatrix<double> get_Msigma();
	/// return the matrix of size (#primalVertex, #primalVertexBoundary) whose under part is an identity matrix.
	Eigen::SparseMatrix<double> get_bndryInclOp();
	/// return Q matrix
	Eigen::SparseMatrix<int> setupOperatorQ();
	/// return the square matrix of size #primalEdge       d = Meps * b     (consistent to the paper)
	Eigen::SparseMatrix<double> setupOperatorMeps();
	/// return the square matrix of size #primalFace       h = Mnu  * e     (consistent to the paper)
	Eigen::SparseMatrix<double> setupOperatorMnu();


private:
	void init();

	Eigen::SparseMatrix<double> Mnu;
	Eigen::SparseMatrix<double> Meps;
	Eigen::SparseMatrix<double> Msigma;
	Eigen::SparseMatrix<double> bndryInclOp;


};

