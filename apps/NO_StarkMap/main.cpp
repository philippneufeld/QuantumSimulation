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

    // constexpr double dE = 3.25 * EnergyInverseCm_v;
    constexpr double dE = 1.0 * EnergyInverseCm_v;
    
    NitricOxide molecule;
    RydbergDiatomicState_t state(38, 2, 2, 4, 0);
    DiatomicStarkMap starkMap(molecule, state, 20, 70, 2, dE);

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
    // auto eField = VectorXd::LinSpaced(1500, 0.0, 25.0); // V cm^-1
    auto eField = VectorXd::LinSpaced(500, 0.0, 2.5); // V cm^-1
    // auto eField = VectorXd::LinSpaced(20, 0.0, 0.025); // V cm^-1

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
    VectorXd prevEnergies, prevEnergies2, extrEnergies;

    std::set<int> matched;
    std::vector<int> indices;
    ProgressBar progress2(eField.size() - 1);
    for(int i=0; i<0*eField.size(); i++)
    {
        auto data = storageThread.LoadData(i);
        auto [idx, ef, energies, states, character] = data;

        if (i == 0)
        {
            // prepare ordering auxilliary variables
            prevEnergies = energies;
            prevEnergies2 = prevEnergies;
            prevStates = states;
            indices.resize(states.rows());
            for (int i=0; i<states.rows();i++)
                indices[i] = i;
            continue;
        }

        matched.clear();
        
        VectorXd energiesOrdered(states.rows());
        MatrixXd statesOrdered(states.rows(), states.cols());
        Matrix<double, Dynamic, 4> characterOrdered(states.rows(), 4);

        MatrixXd overlaps = (prevStates.transpose() * states).cwiseAbs().array().transpose().matrix();
        MatrixXd energyDiffs(overlaps.rows(), overlaps.cols());
        for (int i1=0; i1<energyDiffs.rows(); i1++)
        {
            for (int i2=0; i2<energyDiffs.rows(); i2++)
            {
                extrEnergies = prevEnergies2;
                if (i > 1)
                    extrEnergies += (prevEnergies - prevEnergies2) / (eField[i-1] - eField[i - 2]) * (eField[i] - eField[i-2]);
                energyDiffs(i1, i2) = std::abs(energies[i2] - extrEnergies[i1]);
            }
        }
        energyDiffs /= energyDiffs.maxCoeff();
        MatrixXd combinedCrits = overlaps - energyDiffs.cwiseSqrt();

        std::vector<double> overlapMax(states.rows());
        std::vector<double> combinedMax(states.rows());
        std::vector<int> overlapIndices(states.rows());
        std::vector<int> stateIndices(states.rows());
        for(int j=0; j<states.rows(); j++)
        {
            stateIndices[j]= j;
            overlapMax[j] = overlaps.row(j).maxCoeff(&overlapIndices[j]);
            combinedMax[j] = combinedCrits.row(j).maxCoeff();
        }

        std::stable_sort(stateIndices.begin(), stateIndices.end(), 
            [&](int i1, int i2){ return overlapMax[i1]>overlapMax[i2]; });
        
        constexpr double threshold = 0.8; // 0.70710678118;
        auto it = stateIndices.begin();
        for(; it < stateIndices.end() && overlapMax[*it] > threshold; it++);
        std::stable_sort(it, stateIndices.end(), 
            [&](int i1, int i2){ return combinedMax[i1]>combinedMax[i2]; });

        for(int j: stateIndices)
        {
            auto overlap = overlaps.row(j);
            auto combinedCrit = combinedCrits.row(j);
            idx = overlapIndices[j];

            // index is already taken? Find greates overlap state that is
            // not already taken
            if (overlapMax[j] <= threshold || matched.find(idx) != matched.end())
            {
                std::stable_sort(indices.begin(), indices.end(), 
                    [&](int i1, int i2){ return combinedCrit[i1]<combinedCrit[i2]; });

                int k = indices.size();
                while (matched.find(indices[--k]) != matched.end() && k >= 0);
                idx = indices[k];
            }

            matched.insert(idx);
            energiesOrdered(j) = energies(idx);
            statesOrdered.row(j) = states.row(idx);
            characterOrdered.row(j) = character.row(idx);
        }

        // if (i>1) energiesOrdered = extrEnergies;

        prevEnergies2 = prevEnergies;
        prevEnergies = energiesOrdered;
        prevStates = statesOrdered;
        
        storageThread.ChangeData(i, ef, energiesOrdered, statesOrdered, characterOrdered);
        progress2.IncrementCount();
    }

    std::cout << "Data stored successfully (" << filename << ")" << std::endl;

    return 0;
}
