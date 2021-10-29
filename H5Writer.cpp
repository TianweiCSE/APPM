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

//void H5Writer::writeData(const Eigen::MatrixXd & matrix, const std::string & dataname)
//{
//	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
//	const int rank = 2;
//	hsize_t dims[rank];
//	dims[1] = matrix.rows();
//	dims[0] = matrix.cols();
//	DataSpace dataspace(rank, dims);
//	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_DOUBLE, dataspace);
//	dataset.write(matrix.transpose().data(), PredType::NATIVE_DOUBLE);
//}
//
//void H5Writer::writeData(const Eigen::VectorXi & vector, const std::string & dataname)
//{
//	std::cout << "Write data to " << this->filename << ": " << dataname << std::endl;
//	const int rank = 1;
//	hsize_t dims[rank];
//	dims[0] = vector.cols();
//	DataSpace dataspace(rank, dims);
//	DataSet dataset = file.createDataSet(dataname.c_str(), PredType::NATIVE_INT, dataspace);
//	dataset.write(vector.data(), PredType::NATIVE_INT);
//}
