// Philipp Neufeld, 2021-2022

#include <iostream>
#include <algorithm>
#include <limits>
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

unsigned int GetNumberOfCalcThreads()
{
    std::string hostname = GetHostname();
    unsigned int logical = std::thread::hardware_concurrency();

    // use only half the cores on the calc* machines (hyperthreading)
    if (hostname.find("calc") != hostname.npos)
        return (logical / 2) - 1;
    else
        return logical;
}

int main(int argc, const char* argv[])
{
    // constexpr double dE = 3.5 * EnergyInverseCm_v;
    constexpr double dE = 0.05 * EnergyInverseCm_v;
    
    NitricOxide molecule;
    RydbergDiatomicState_t state(38, 2, 2, 0, 0);
    DiatomicStarkMap starkMap(molecule, state, 30, 60, 3, dE);

    // analyze basis
    std::vector<RydbergDiatomicState_t> basis = starkMap.GetBasis();
    int basisSize = basis.size();
    int minN = std::numeric_limits<int>::max();
    int maxN = std::numeric_limits<int>::min();
    int maxR = 0;
    int maxL = 0;
    for (const auto& state: basis)
    {
        auto [n, l, R, N, mN] = state;
        minN = std::min(minN, n);
        maxN = std::max(maxN, n);
        maxL = std::max(maxL, l);
        maxR = std::max(maxR, R);
    }
    std::cout << "Basis size: " << basisSize << std::endl;

    // generate character matrices
    MatrixXd nCharMatrix(maxN - minN + 1, basisSize);
    MatrixXd lCharMatrix(maxL + 1, basisSize);
    MatrixXd RCharMatrix(maxR + 1, basisSize);
    for (int i = 0; i < basisSize; i++)
    {
        auto [n, l, R, N, mN] = basis[i];
        nCharMatrix(n-minN, i) = 1.0;
        lCharMatrix(l, i) = 1.0;
        RCharMatrix(R, i) = 1.0;
    }

    // generate filename
    std::string path;
    path += GetHomeDirSubfolderPath("remote_home") + "/Masterarbeit/06_StarkMap/03_NO/";
    path += GenerateFilename("NOStarkMap") + ".h5";
    
    // calculate stark map
    auto eField = VectorXd::LinSpaced(1500, 0.0, 25.0); // V cm^-1

    // ThreadPool pool(GetNumberOfCalcThreads());
    ThreadPool pool(1);
    // StorageThread storageThread(path, state, starkMap.GetBasis(), dE, eField.size());
    ProgressBar progress(eField.size());

    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            auto [energies, states] = starkMap.GetEnergiesAndStates(eField[i] * 100);
            energies /= EnergyInverseCm_v;
            
            // get character of the states
            Matrix<double, Dynamic, 3> nlRChar(basisSize, 3);

            MatrixXd tmp1 = nCharMatrix * states;
            MatrixXd tmp2 = lCharMatrix * states;
            MatrixXd tmp3 = RCharMatrix * states;

            std::cout << tmp1 << std::endl << std::endl;
            std::cout << tmp2 << std::endl << std::endl;
            std::cout << tmp3 << std::endl << std::endl;

            // storageThread.AddData(i, eField[i], energies, states);
            progress.IncrementCount();
        });
    }

    progress.WaitUntilFinished();
    std::cout << "Finishing up..." << std::endl;
    // storageThread.WaitUntilFinished();

    std::cout << "Finished. Data saved to " << path << std::endl;

    return 0;
}
