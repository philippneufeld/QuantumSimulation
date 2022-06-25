// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Util/DataFile.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char *argv[])
{
    
    DataFile file;
    file.Open("/home/pneufeld/test.h5", DataFile_TRUNCATE);
    auto root = file.OpenRootGroup();

    Matrix<double, 2, 3> test1;
    test1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    
    Vector4d test2;
    test2 << 1.0, 2.0, 3.0, 4.0;

    RowVector4d test3;
    test3 << 1.0, 2.0, 3.0, 4.0;

    Matrix<std::complex<double>, 2, 3> test4;
    test4 = std::complex<double>{1, -1} * test1; 

    std::cout << test1 << std::endl;
    std::cout << test2 << std::endl;
    std::cout << test3 << std::endl;
    std::cout << test4 << std::endl;
    std::cout << std::endl;

    auto dset = root.CreateDataset("test1", test1);
    dset.SetAttribute("test2", test2);
    root.SetAttribute("test3", test3);
    root.SetAttribute("test4", test4);
    root.CreateDataset("test4", test4);
    root.SetAttribute("test5", 5.0);

    // MULTIPLE WRITE TEST
    std::cout << (dset.Set(test3) ? "NOT CORRECT" : "CORRECT") << std::endl;
    std::cout << (dset.SetAttribute("test2", test3) ? "NOT CORRECT" : "CORRECT") << std::endl;
    std::cout << (dset.SetAttribute("test2", test2) ? "CORRECT" : "NOT CORRECT") << std::endl;
    
    std::cout << root.GetDataset("test1").template Get<decltype(test1)>() << std::endl;
    std::cout << root.GetDataset("test1").template GetAttribute<Vector4d>("test2") << std::endl;
    std::cout << root.template GetAttribute<RowVector4d>("test3") << std::endl;
    std::cout << root.template GetAttribute<decltype(test4)>("test4") << std::endl;
    std::cout << root.GetDataset("test4").template Get<decltype(test4)>() << std::endl;
    std::cout << root.template GetAttribute<double>("test5") << std::endl;
    
    return 0;
}
