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

const int H5Reader::readDataSize(const std::string & dataname) const
{
	DataSet dataset = file.openDataSet(dataname.c_str());
	DataSpace dataspace = dataset.getSpace();
	const int rank = dataspace.getSimpleExtentNdims();
	assert(rank <= 1);
	hsize_t *dims = new hsize_t(rank);
	dataspace.getSimpleExtentDims(dims);
	int result = dims[0];
	delete dims;
	return result;
}
