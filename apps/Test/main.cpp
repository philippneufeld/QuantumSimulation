// Philipp Neufeld, 2021-2022

#include <iostream>
#include <QSim/Util/DataFile2.h>

#include <H5Cpp.h>

#define TEST_WRITE

int main(int argc, const char* argv[])
{
    float data[4][4]; // buffer for data to write
    
#ifdef TEST_WRITE
    H5::H5File file("test.h5", H5F_ACC_TRUNC);
    
    // Create the data space for the dataset.
    hsize_t dims[2]; // dataset dimensions
    dims[0] = 4;
    dims[1] = 4;
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset = file.createDataSet("test", H5::PredType::NATIVE_FLOAT, dataspace);
    
    for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++)
            data[j][i] = i == j ? 1 : 0;// i * 6 + j + 1;
    dataset.write(data, H5::PredType::NATIVE_FLOAT);
#else
    H5::H5File file("test.h5", H5F_ACC_RDWR );
    H5::DataSet dataset = file.openDataSet("test");    
    dataset.read(data, H5::PredType::NATIVE_FLOAT);
#endif

    return 0;
}
