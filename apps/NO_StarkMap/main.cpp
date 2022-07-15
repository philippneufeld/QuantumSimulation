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

#include "IOThread.h"
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



class NOStarkMapApp
{
public:
    NOStarkMapApp(const std::string path)
        : m_path(path), m_ioThread(path) { }

    void RunCalculation(const VectorXd& eField)
    {
        // constexpr double dE = 3.25 * EnergyInverseCm_v;
        constexpr double dE = 8.0 * EnergyInverseCm_v;

        // initialize calculation
        NitricOxide molecule;
        RydbergDiatomicState_t state(33, 2, 2, 0, 0);
        DiatomicStarkMap starkMap(molecule, state, 20, 70, 2, dE);

        // process basis
        std::vector<RydbergDiatomicState_t> basis = starkMap.GetBasis();
        AnalyzeBasis(basis);

        std::cout << "Basis size: " << basis.size() << std::endl;

        // start i/o thread
        m_ioThread.Start(state, basis, dE);

        // do the actual calculation
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
                auto character = DetermineStateCharacter(states);

                // Store data
                m_ioThread.StoreData(i, eField[i], energies, states, character);
                progress.IncrementCount();
            });
        }

        progress.WaitUntilFinished();
        std::cout << "Waiting for I/O to complete..." << std::endl;
        m_ioThread.WaitUntilFinished();

        // do post processing
        std::cout << "Overlapping states..." << std::endl;
        MatchStates(eField);

        std::cout << "Data stored successfully (" << m_path << ")" << std::endl;
    }

protected:

    void AnalyzeBasis(const std::vector<RydbergDiatomicState_t>& basis)
    {
        int basisSize = basis.size();
        m_minn = std::numeric_limits<int>::max();
        int maxn = std::numeric_limits<int>::min();
        int maxR = std::numeric_limits<int>::min();
        int maxL = std::numeric_limits<int>::min();
        int maxN = std::numeric_limits<int>::min();
        for (const auto& state: basis)
        {
            auto [n, l, R, N, mN] = state;
            m_minn = std::min(m_minn, n);
            maxn = std::max(maxn, n);
            maxL = std::max(maxL, l);
            maxR = std::max(maxR, R);
            maxN = std::max(maxN, N);    
        }

        // generate character matrices
        m_nCharMatrix = MatrixXd::Zero(maxn - m_minn + 1, basisSize);
        m_lCharMatrix = MatrixXd::Zero(maxL + 1, basisSize);
        m_RCharMatrix = MatrixXd::Zero(maxR + 1, basisSize);
        m_NCharMatrix = MatrixXd::Zero(maxN + 1, basisSize);
        for (int i = 0; i < basisSize; i++)
        {
            auto [n, l, R, N, mN] = basis[i];
            m_nCharMatrix(n-m_minn, i) = 1.0;
            m_lCharMatrix(l, i) = 1.0;
            m_RCharMatrix(R, i) = 1.0;
            m_NCharMatrix(N, i) = 1.0;
        }
    }

    Matrix<double, Dynamic, 4> DetermineStateCharacter(const MatrixXd& states)
    {
        Matrix<double, Dynamic, 4> character(states.rows(), 4);
        auto statesSq = states.cwiseAbs2().eval();
        MatrixXd nMat = m_nCharMatrix * statesSq;
        MatrixXd lMat = m_lCharMatrix * statesSq;
        MatrixXd RMat = m_RCharMatrix * statesSq;
        MatrixXd NMat = m_NCharMatrix * statesSq;

        for (int j=0; j<states.rows(); j++)
        {
            int nIdx=0, lIdx=0, RIdx=0, NIdx=0;
            nMat.col(j).maxCoeff(&nIdx);
            lMat.col(j).maxCoeff(&lIdx);
            RMat.col(j).maxCoeff(&RIdx);
            NMat.col(j).maxCoeff(&NIdx);
            
            character(j, 0) = m_minn + nIdx;
            character(j, 1) = lIdx;
            character(j, 2) = RIdx;
            character(j, 3) = NIdx;
        }

        return character;
    }
    
    void MatchStates(const VectorXd& eField)
    {
        MatrixXd prevStates;
        VectorXd prevEnergies, prevEnergies2, extrEnergies;

        std::set<int> matched;
        std::vector<int> indices;
        ProgressBar progress(eField.size() - 1);

        auto nextDataPromise = m_ioThread.LoadData(0);
        for(int i=0; i<eField.size(); i++)
        {
            auto [idx, ef, energies, states, character] = nextDataPromise.get();
            if (i < eField.size() - 1)
                nextDataPromise = m_ioThread.LoadData(i + 1);

            // prepare ordering auxilliary variables
            if (i == 0)
            {
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

            // calculate heuristic matching criteria
            MatrixXd overlaps = (states.transpose() * prevStates).cwiseAbs();
            MatrixXd energyDiffs(overlaps.rows(), overlaps.cols());
            for (int i1=0; i1<energyDiffs.rows(); i1++)
            {
                for (int i2=0; i2<energyDiffs.rows(); i2++)
                {
                    extrEnergies = prevEnergies2;
                    if (i > 1)
                        extrEnergies += (prevEnergies - prevEnergies2) / (eField[i-1] - eField[i - 2]) * (eField[i] - eField[i-2]);
                    energyDiffs(i1, i2) = std::abs(energies[i1] - extrEnergies[i2]);
                }
            }
            energyDiffs /= energyDiffs.maxCoeff();
            MatrixXd heuristics = overlaps - energyDiffs.cwiseSqrt();

            // determine order in which states are processed
            std::vector<double> overlapMax(states.rows());
            std::vector<double> heuristicMax(states.rows());
            std::vector<int> overlapIndices(states.rows());
            std::vector<int> stateIndices(states.rows());
            for(int j=0; j<states.rows(); j++)
            {
                stateIndices[j]= j;
                overlapMax[j] = overlaps.col(j).maxCoeff(&overlapIndices[j]);
                heuristicMax[j] = heuristics.col(j).maxCoeff();
            }

            std::stable_sort(stateIndices.begin(), stateIndices.end(), 
                [&](int i1, int i2){ return overlapMax[i1]>overlapMax[i2]; });
            
            constexpr double threshold = 1 + 0.70710678118;
            auto it = stateIndices.begin();
            for(; it < stateIndices.end() && overlapMax[*it] > threshold; it++);
            std::stable_sort(it, stateIndices.end(), 
                [&](int i1, int i2){ return heuristicMax[i1]>heuristicMax[i2]; });

            // do the actual matching procedure
            for(int j: stateIndices)
            {
                auto overlap = overlaps.col(j);
                auto heuristic = heuristics.col(j);
                idx = overlapIndices[j];

                // index is already taken? Find greates overlap state that is
                // not already taken
                if (overlapMax[j] <= threshold || matched.find(idx) != matched.end())
                {
                    std::stable_sort(indices.begin(), indices.end(), 
                        [&](int i1, int i2){ return heuristic[i1]<heuristic[i2]; });
                    int k = indices.size();
                    while (matched.find(indices[--k]) != matched.end() && k >= 0);
                    idx = indices[k];
                }

                matched.insert(idx);
                energiesOrdered[j] = energies[idx];
                statesOrdered.col(j) = states.col(idx);
                characterOrdered.row(j) = character.row(idx);
            }

            // store processed data back to file
            m_ioThread.StoreData(i, ef, energiesOrdered, statesOrdered, characterOrdered);

            // shift auxilliary variables
            prevEnergies2 = prevEnergies;
            prevEnergies = energiesOrdered;
            prevStates = statesOrdered;
            
            progress.IncrementCount();
        }
    }
    

private:
    std::string m_path;
    IOThread m_ioThread;

    // state characterization
    int m_minn;
    MatrixXd m_nCharMatrix;
    MatrixXd m_lCharMatrix;
    MatrixXd m_RCharMatrix;
    MatrixXd m_NCharMatrix;
};


int main(int argc, const char* argv[])
{
    // parsing command line arguments
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

    // Run calculation
    NOStarkMapApp app(args.GetOptionStringValue("file"));
    VectorXd eField = VectorXd::LinSpaced(200, 0.0, 16.0); // V cm^-1
    app.RunCalculation(eField);

    return 0;
}
