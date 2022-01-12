// Philipp Neufeld, 2021-2022

#include <iostream>
#include <QSim/Util/DataFile2.h>

#include <H5Cpp.h>

int main(int argc, const char* argv[])
{
    H5::H5File file("test.h5", H5F_ACC_TRUNC);
    
    // Create the data space for the dataset.
    hsize_t dims[2]; // dataset dimensions
    dims[0] = 4;
    dims[1] = 4;
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset = file.createDataSet("test", H5::PredType::STD_I32BE, dataspace);
    
    int data[4][4]; // buffer for data to write
    for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++)
            data[j][i] = i * 6 + j + 1;
    dataset.write(data, H5::PredType::NATIVE_INT);

    return 0;
}
