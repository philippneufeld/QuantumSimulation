// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergDiatomic.h>
#include <QSim/Rydberg/DiatomicStarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#include "StorageThread.h"
#include <QSim/Util/PathUtil.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    constexpr double dE = 3.5 * EnergyInverseCm_v;
    
    NitricOxide molecule;
    RydbergDiatomicState_t state(38, 2, 2, 0, 0);
    DiatomicStarkMap starkMap(molecule, state, 30, 60, 3, dE);
    
    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

    auto eField = VectorXd::LinSpaced(1500, 0.0, 25.0); // V cm^-1

    // generate filename
    std::string path;
    path += GetHomeDirSubfolderPath("remote_home") + "/Masterarbeit/06_StarkMap/03_NO/";
    path += GenerateFilename("NOStarkMap") + ".h5";
    
    ThreadPool pool;
    StorageThread storageThread(path, state, starkMap.GetBasis(), dE, eField.size());
    ProgressBar progress(eField.size());

    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            auto [energies, states] = starkMap.GetEnergiesAndStates(eField[i] * 100);
            energies /= EnergyInverseCm_v;
            
            storageThread.AddData(i, eField[i], energies, states);
            progress.IncrementCount();
        });
    }

    progress.WaitUntilFinished();
    std::cout << "Finishing up..." << std::endl;
    storageThread.WaitUntilFinished();

    std::cout << "Finished. Data saved to " << path << std::endl;

    return 0;
}
