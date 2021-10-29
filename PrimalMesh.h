#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Mesh.h"
class PrimalMesh :
	public Mesh
{
public:
	class PrimalMeshParams {
	public:
		PrimalMeshParams();
		PrimalMeshParams(const std::string & filename);

		const int getRefinements() const;
		const int getAxialLayers() const;
		const int getOuterLayers() const;
		const double getElectrodeRadius() const;
		const double getZmax() const;

		friend std::ostream & operator<<(std::ostream & os, const PrimalMeshParams & obj);

	private:
		int nAxialLayers = 10;
		int nRefinements = 2;
		int nOuterLayers = 2;
		double electrodeRadius = 0.35;  /// NOTICE: The fluid radius is set one by default. The electro radius is different from fluid radius!
		double zmax = 1;

		void readParameters(const std::string & filename);
	};

	PrimalMesh();
	PrimalMesh(const std::string & meshPrefix);
	PrimalMesh(const PrimalMeshParams & p);
	~PrimalMesh();

	void init();

private:
	// Mesh parameters
	PrimalMeshParams params;

	/// Add the vertices, edges and faces of a hexagon centered at (0,0,zValue) (x-y plane)
	void init_hexagon(const double zValue);
	/// None
	void init_triangle();
	/// Build new f2v map; Remove old mesh; Bluid new mesh (x-y plane)
	void refineMesh(const int nRefinements);
	/// do 'outerMeshExtrude_prisms' by 'nLayers' times (x-y plane)
	void outerMeshExtrude(const int nLayers);
	/// Add two layers of triangles at the outer side of hexagon (x-y plane)
	void outerMeshExtrude_triangles();
	/// Add one layer of quatrilaterals at the out side of hexagon (x-y plane)
	void outerMeshExtrude_prisms();
	/// Get tensor product mesh on z-axis
	void extrudeMesh(const int nLayers, const double zmax);

	Eigen::Matrix3Xi refine_triangles();
	Eigen::Matrix3Xi refine_triangles_specialCorners();

	void test_quadFace();

	/// Sort vertices such that they have the sequence: INNER, BOUNDARY, TERMINAL. (type assigned)
	void sortVertices(const double electrodeRadius);
	/// Sort edges such that they have the sequence: Interior, InteriorToBoundary, Boundary. (type assigned)
	void sortEdges();
	/// Sort faces such that they have the sequence: inner, boundary. (type not assigned)
	void sortFaces();
	/// Sort cells such that they have the sequence: inner, outer. (type not assigned)
	void sortCells();
	
	void validateParameters();
	/// Check that all vertices have z-coordinate equal to z0.
	void check_zCoord(const double z0);
};

