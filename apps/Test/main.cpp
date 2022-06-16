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

    std::cout << test1 << std::endl;
    std::cout << test2 << std::endl;
    std::cout << test3 << std::endl;

    root.SetAttribute("test1", test1);
    root.SetAttribute("test2", test2);
    root.SetAttribute("test3", test3);
    root.SetAttribute("test4", 5.0);

    std::cout << root.template GetAttribute<decltype(test1)>("test1") << std::endl;
    std::cout << root.template GetAttribute<Vector4d>("test2") << std::endl;
    std::cout << root.template GetAttribute<RowVector4d>("test3") << std::endl;
    std::cout << root.template GetAttribute<double>("test4") << std::endl;
    
    // double data1[] = { 1.0, 2.0, 3.0, 4.0 };
    // root.CreateAttribute("attr1", {2, 2});
    // root.StoreAttribute("attr1", data1);
    // root.CreateAttribute("attr2", {1, });
    // root.StoreAttribute("attr2", data1);

    return 0;
}
