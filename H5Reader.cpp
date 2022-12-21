#include "H5Reader.h"



H5Reader::H5Reader()
{
}

H5Reader::H5Reader(const std::string & filename)
{
	this->filename = filename;
	this->file = H5File(filename.c_str(), H5F_ACC_RDONLY);
}


H5Reader::~H5Reader()
{
}

int H5Reader::readVectorDataSize(const std::string & dataname) const {
	DataSet dataset = file.openDataSet(dataname.c_str());
	DataSpace dataspace = dataset.getSpace();
	const int rank = dataspace.getSimpleExtentNdims();
	assert(rank <= 1);
	hsize_t *dims = new hsize_t[rank];
	dataspace.getSimpleExtentDims(dims);
	int result = dims[0];
	delete[] dims;
	return result;
}

int H5Reader::readMatrixDataSize(const std::string & dataname) const {
	DataSet dataset = file.openDataSet(dataname.c_str());
	DataSpace dataspace = dataset.getSpace();
	const int rank = dataspace.getSimpleExtentNdims();
	assert(rank == 2);
	hsize_t *dims = new hsize_t[rank];
	dataspace.getSimpleExtentDims(dims);
	int result = dims[0];
	delete[] dims;
	return result;
}

Eigen::VectorXd H5Reader::readVectorData(const std::string & dataname) const {
	DataSet dataset = file.openDataSet(dataname.c_str());
	DataSpace dataspace = dataset.getSpace();
	const int rank = dataspace.getSimpleExtentNdims();
	assert(rank == 1);
	hsize_t *dims = new hsize_t[rank];
	dataspace.getSimpleExtentDims(dims);
	Eigen::VectorXd vec(dims[0]);
	dataset.read(vec.data(), PredType::NATIVE_DOUBLE);
	return vec;
}

Eigen::MatrixXd H5Reader::readMatrixData(const std::string & dataname) const {
	DataSet dataset = file.openDataSet(dataname.c_str());
	DataSpace dataspace = dataset.getSpace();
	const int rank = dataspace.getSimpleExtentNdims();
	assert(rank == 2);
	hsize_t *dims = new hsize_t[rank];
	dataspace.getSimpleExtentDims(dims);
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat((int)dims[0], (int)dims[1]);
	dataset.read(mat.data(), PredType::NATIVE_DOUBLE);
	return mat;
}