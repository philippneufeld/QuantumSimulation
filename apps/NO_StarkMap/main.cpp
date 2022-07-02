// Philipp Neufeld, 2021-2022

#include <set>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <filesystem>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergDiatomic.h>
#include <QSim/Rydberg/DiatomicStarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>
#include <QSim/Util/Argparse.h>

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
    ArgumentParser argparse;

    std::string defaultFilename = GetHomeDirSubfolderPath("remote_home") 
        + "/Masterarbeit/06_StarkMap/03_NO/" + GenerateFilename("NOStarkMap") + ".h5";
    argparse.AddOptionDefault("file", "HDF5 file path", defaultFilename);
    argparse.AddOption("help", "Print this help string");

    auto args = argparse.Parse(argc, argv);
    
    if (args.IsError())
    {
        std::cout << args.GetError() << std::endl;
        return 1;
    }
    else if (args.IsOptionPresent("help"))
    {
        std::cout << argparse.GetHelpString() << std::endl;
        return 0;
    }

    std::string filename = args.GetOptionStringValue("file");

    constexpr double dE = 3.25 * EnergyInverseCm_v;
    
    NitricOxide molecule;
    RydbergDiatomicState_t state(38, 2, 2, 4, 0);
    DiatomicStarkMap starkMap(molecule, state, 20, 70, 4, dE);

    // analyze basis
    std::vector<RydbergDiatomicState_t> basis = starkMap.GetBasis();
    int basisSize = basis.size();
    int minn = std::numeric_limits<int>::max();
    int maxn = std::numeric_limits<int>::min();
    int maxR = std::numeric_limits<int>::min();
    int maxL = std::numeric_limits<int>::min();
    int minN = std::numeric_limits<int>::max();
    int maxN = std::numeric_limits<int>::min();
    for (const auto& state: basis)
    {
        auto [n, l, R, N, mN] = state;
        minn = std::min(minn, n);
        maxn = std::max(maxn, n);
        maxL = std::max(maxL, l);
        maxR = std::max(maxR, R);
        minN = std::min(minN, N);
        maxN = std::max(maxN, N);    
    }
    std::cout << "Basis size: " << basisSize << std::endl;

    // generate character matrices
    MatrixXd nCharMatrix = MatrixXd::Zero(maxn - minn + 1, basisSize);
    MatrixXd lCharMatrix = MatrixXd::Zero(maxL + 1, basisSize);
    MatrixXd RCharMatrix = MatrixXd::Zero(maxR + 1, basisSize);
    MatrixXd NCharMatrix = MatrixXd::Zero(maxN - minN + 1, basisSize);
    for (int i = 0; i < basisSize; i++)
    {
        auto [n, l, R, N, mN] = basis[i];
        nCharMatrix(n-minn, i) = 1.0;
        lCharMatrix(l, i) = 1.0;
        RCharMatrix(R, i) = 1.0;
        NCharMatrix(N-minN, i) = 1.0;
    }
    
    // calculate stark map
    auto eField = VectorXd::LinSpaced(1500, 0.0, 25.0); // V cm^-1
    // auto eField = VectorXd::LinSpaced(180, 0.0, 25.0); // V cm^-1

    StorageThread storageThread(filename, state, starkMap.GetBasis(), dE, eField.size());
    ThreadPool pool(GetNumberOfCalcThreads());
    ProgressBar progress(eField.size());

    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i]()
        {
            // calculate eigenenergies and eigenstates
            auto [energies, states] = starkMap.GetEnergiesAndStates(eField[i] * 100);
            energies /= EnergyInverseCm_v;
            
            // get character of the states
            Matrix<double, Dynamic, 4> character(basisSize, 4);
            auto statesSq = states.cwiseAbs2().eval();
            MatrixXd nMat = nCharMatrix * statesSq;
            MatrixXd lMat = lCharMatrix * statesSq;
            MatrixXd RMat = RCharMatrix * statesSq;
            MatrixXd NMat = NCharMatrix * statesSq;

            for (int j=0; j<basisSize; j++)
            {
                int nIdx=0, lIdx=0, RIdx=0, NIdx=0;
                nMat.col(j).maxCoeff(&nIdx);
                lMat.col(j).maxCoeff(&lIdx);
                RMat.col(j).maxCoeff(&RIdx);
                NMat.col(j).maxCoeff(&NIdx);
                
                character(j, 0) = minn + nIdx;
                character(j, 1) = lIdx;
                character(j, 2) = RIdx;
                character(j, 3) = minN + NIdx;
            }

            // Store data
            storageThread.StoreData(i, eField[i], energies, states, character);
            progress.IncrementCount();
        });
    }

    progress.WaitUntilFinished();
    std::cout << "Waiting for I/O to complete..." << std::endl;
    storageThread.WaitUntilFinished();

    // do post processing
    std::cout << "Do post processing..." << std::endl;

    MatrixXd prevStates;
    std::set<int> matched;
    std::vector<int> indices;
    ProgressBar progress2(eField.size());
    for(int i=0; i<eField.size(); i++)
    {
        auto data = storageThread.LoadData(i);
        auto [idx, eField, energies, states, character] = data;

        if (i == 0)
        {
            // prepare ordering auxilliary variables
            prevStates = states;
            indices.reserve(states.rows());
            for (int i=0; i<states.rows();i++)
                indices[i] = i;
            continue;
        }

        matched.clear();
        MatrixXd statesOrdered(states.rows(), states.cols());
        MatrixXd overlaps = (prevStates.transpose() * states).cwiseAbs();
        for(int j=0; j<states.rows(); j++)
        {
            auto overlap = overlaps.row(j);
            int idx = 0;
            int val = overlap.maxCoeff(&idx);

            // index is already taken? Find greates overlap state that is
            // not already taken
            if (matched.find(idx) != matched.end())
            {
                std::stable_sort(indices.begin(), indices.end(), 
                    [&](int i1, int i2){ return overlap[i1]<overlap[i2]; });

                int k = indices.size() - 1;
                while (matched.find(indices[--k]) == matched.end());
                idx = indices[k];
            }

            matched.insert(idx);
            statesOrdered.row(j) = states.row(idx);
        }
        
        storageThread.AdjustStates(i, statesOrdered);
        progress2.IncrementCount();
    }

    std::cout << "Data stored successfully (" << filename << ")" << std::endl;

    return 0;
}
