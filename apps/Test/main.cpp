// Philipp Neufeld, 2021-2022

#include <iostream>
#include <QSim/Util/ScopeGuard.h>
#include <QSim/Util/DataFile.h>
#include <QSim/Math/Matrix.h>
#include <QSim/Util/SimulationApp.h>

using namespace QSim;

int main(int argc, const char* argv[])
{
    DataFile file;
    file.Open("test.h5", DataFile_DEFAULT);
    DataFileGroup root = file.OpenRootGroup();

    auto test = CreateIdentity<double>(5);
    DataFileDataset testDataset;
    if (!root.DoesDatasetExist("test")) 
        testDataset = root.CreateDataset("test", {test.Rows(), test.Cols()});
    else
        testDataset = root.GetDataset("test");
    testDataset.Store(test.Data());

    DataFileGroup group1;
    if (!root.DoesSubgroupExist("group1")) 
        group1 = root.CreateSubgroup("group1");
    else
        group1 = root.GetSubgroup("group1");

    auto test2 = test * 2;
    DataFileDataset test2Dataset;
    if (!group1.DoesDatasetExist("test2")) 
        test2Dataset = group1.CreateDataset("test2", {test.Rows(), test.Cols()});
    else
        test2Dataset = group1.GetDataset("test2");
    test2Dataset.Store(test2.Data());
    
    double myattr[1] = { 127.0 };
    if (!group1.DoesAttributeExist("myattr"))
        group1.CreateAttribute("myattr", {1});
    group1.StoreAttribute("myattr", myattr);

    auto e1 = root.EnumerateSubgroups();
    auto e2 = root.EnumerateDatasets();

    return 0;
}
