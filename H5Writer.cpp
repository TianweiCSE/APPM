#include "H5Writer.h"

H5Writer::H5Writer()
{
}

H5Writer::H5Writer(const std::string & filename)
{
	std::cout << "Create HDF5 file: " << filename << std::endl;
	const int length = 3;
	assert(filename.size() > 3);
	const int pos = filename.size() - 3;
	std::string substr = filename.substr(pos, length);
	assert(filename.compare(pos, length, ".h5") == 0);

	this->filename = filename;
	this->file = H5File(filename.c_str(), H5F_ACC_TRUNC);
}


H5Writer::~H5Writer()
{
}

void H5Writer::writeStdVector(const std::vector<int> & vector, const std::string & dataname){
	std::cout << "Write int data to " << this->filename << ": " << dataname << std::endl;
	const int rank = 1;
	hsize_t dims[rank];
	dims[0] = vector.size();
	DataSpace dataspace(rank, dims);
	const PredType type = PredType::NATIVE_INT;
	DataSet dataset = file.createDataSet(dataname.c_str(), type, dataspace);
	dataset.write(vector.data(), type);
}

void H5Writer::writeIntMatrix(const Eigen::MatrixXi & matrix, const std::string & dataname){
	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
	hsize_t dims[2];
	dims[0] = matrix.cols();
	dims[1] = matrix.rows();

	DataSpace dataspace(2, dims);
	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_INT, dataspace);
	dataset.write(matrix.transpose().data(), PredType::NATIVE_INT);
}

void H5Writer::writeDoubleMatrix(const Eigen::MatrixXd & matrix, const std::string & dataname){
	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
	hsize_t dims[2];
	dims[0] = matrix.cols();
	dims[1] = matrix.rows();

	DataSpace dataspace(2, dims);
	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_DOUBLE, dataspace);
	dataset.write(matrix.transpose().data(), PredType::NATIVE_DOUBLE);
}

void H5Writer::writeIntVector(const Eigen::VectorXi & vector, const std::string & dataname){
	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
	hsize_t dims[1];
	dims[0] = vector.size();

	DataSpace dataspace(1, dims);
	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_INT, dataspace);
	dataset.write(vector.data(), PredType::NATIVE_INT);
}

void H5Writer::writeDoubleVector(const Eigen::VectorXd & vector, const std::string & dataname){
	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
	hsize_t dims[1];
	dims[0] = vector.size();

	DataSpace dataspace(1, dims);
	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_DOUBLE, dataspace);
	dataset.write(vector.data(), PredType::NATIVE_DOUBLE);
}


