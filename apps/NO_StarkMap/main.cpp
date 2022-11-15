// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <set>

#include <QSim/Constants.h>
#include <QSim/Rydberg/DiatomicStarkMap.h>
#include <QSim/Rydberg/RydbergDiatomic.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/Argparse.h>
#include <QSim/Util/ProgressBar.h>

#include "IOThread.h"
#include <QSim/Util/PathUtil.h>

using namespace QSim;
using namespace Eigen;

unsigned int GetNumberOfCalcThreads() {
    std::string hostname = GetHostname();
    unsigned int logical = std::thread::hardware_concurrency();

    // use only half the cores on the calc* machines (hyperthreading)
    if (hostname.find("calc") != hostname.npos)
        return (logical / 2) - 1;
    else
        return logical;
}

class NOStarkMapApp {
  public:
    NOStarkMapApp(const std::string path) : m_path(path), m_ioThread(path) {
    }

    void RunCalculation(const VectorXd &eFields, double energy, double dE,
                        int Rmin, int Rmax, int mN, int nMin = 1,
                        int nMax = 100) {
        // initialize calculation
        DiatomicStarkMap starkMap(NitricOxide{}, nMin, nMax, Rmin, Rmax, mN,
                                  energy, dE);

        // process basis
        std::vector<RydbergDiatomicState_t> basis = starkMap.GetBasis();
        AnalyzeBasis(basis);

        std::cout << "Basis size: " << basis.size() << std::endl;

        // start i/o thread
        m_ioThread.Start(basis, energy, dE);

        // do the actual calculation
        ThreadPool pool(GetNumberOfCalcThreads());
        ProgressBar progress(eFields.size());
        for (int i = 0; i < eFields.size(); i++) {
            pool.Submit([&, i = i]() {
                // calculate eigenenergies and eigenstates
                auto [energies, states] =
                    starkMap.GetEnergiesAndStates(eFields[i]);

                // get character of the states
                auto character = DetermineStateCharacter(states);

                // Store data
                m_ioThread.StoreData(i, eFields[i], energies, states,
                                     character);
                progress.IncrementCount();
            });
        }

        progress.WaitUntilFinished();
        std::cout << "Waiting for I/O to complete..." << std::endl;
        m_ioThread.WaitUntilFinished();

        // do post processing
        std::cout << "Overlapping states..." << std::endl;
        MatchStates(eFields, basis.size());

        std::cout << "Data stored successfully (" << m_path << ")" << std::endl;
    }

  protected:
    void AnalyzeBasis(const std::vector<RydbergDiatomicState_t> &basis) {
        int basisSize = basis.size();
        m_minn = std::numeric_limits<int>::max();
        int maxn = std::numeric_limits<int>::min();
        int maxR = std::numeric_limits<int>::min();
        int maxL = std::numeric_limits<int>::min();
        int maxN = std::numeric_limits<int>::min();
        for (const auto &state : basis) {
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
        for (int i = 0; i < basisSize; i++) {
            auto [n, l, R, N, mN] = basis[i];
            m_nCharMatrix(n - m_minn, i) = 1.0;
            m_lCharMatrix(l, i) = 1.0;
            m_RCharMatrix(R, i) = 1.0;
            m_NCharMatrix(N, i) = 1.0;
        }
    }

    Matrix<double, Dynamic, 4> DetermineStateCharacter(const MatrixXd &states) {
        Matrix<double, Dynamic, 4> character(states.rows(), 4);
        auto statesSq = states.cwiseAbs2().eval();
        MatrixXd nMat = m_nCharMatrix * statesSq;
        MatrixXd lMat = m_lCharMatrix * statesSq;
        MatrixXd RMat = m_RCharMatrix * statesSq;
        MatrixXd NMat = m_NCharMatrix * statesSq;

        for (int j = 0; j < states.rows(); j++) {
            int nIdx = 0, lIdx = 0, RIdx = 0, NIdx = 0;
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

    void MatchStates(const VectorXd &eFields, int basisSize) {
        MatrixXd prevStates;
        VectorXd prevEnergies, prevEnergies2, extrEnergies;

        ProgressBar progress(eFields.size() - 1);

        // Allocate memory
        std::vector<double> overlapMax(basisSize);
        std::vector<double> heuristicMax(basisSize);
        std::vector<int> overlapIndices(basisSize);
        std::vector<int> stateIndices(basisSize);

        VectorXd energiesOrdered(basisSize);
        MatrixXd statesOrdered(basisSize, basisSize);
        Matrix<double, Dynamic, 4> characterOrdered(basisSize, 4);

        MatrixXd overlaps(basisSize, basisSize);
        MatrixXd energyDiffs(basisSize, basisSize);
        MatrixXd heuristics;

        std::vector<int> indices(basisSize);
        for (int i = 0; i < indices.size(); i++)
            indices[i] = i;

        auto nextDataPromise = m_ioThread.LoadData(0);
        for (int i = 0; i < eFields.size(); i++) {
            auto [idx, ef, energies, states, character] = nextDataPromise.get();
            if (i < eFields.size() - 1)
                nextDataPromise = m_ioThread.LoadData(i + 1);

            // prepare ordering auxilliary variables
            if (i == 0) {
                prevEnergies = energies;
                prevEnergies2 = prevEnergies;
                prevStates = states;
                continue;
            }

            // calculate heuristic matching criteria
            overlaps = (states.transpose() * prevStates).cwiseAbs();
            energyDiffs(overlaps.rows(), overlaps.cols());
            for (int i1 = 0; i1 < energyDiffs.rows(); i1++) {
                for (int i2 = 0; i2 < energyDiffs.rows(); i2++) {
                    extrEnergies = prevEnergies2;
                    if (i > 1)
                        extrEnergies += (prevEnergies - prevEnergies2) /
                                        (eFields[i - 1] - eFields[i - 2]) *
                                        (eFields[i] - eFields[i - 2]);
                    energyDiffs(i1, i2) =
                        std::abs(energies[i1] - extrEnergies[i2]);
                }
            }
            energyDiffs /= energyDiffs.maxCoeff();
            heuristics = overlaps - energyDiffs.cwiseSqrt();

            // determine order in which states are processed
            for (int j = 0; j < states.rows(); j++) {
                stateIndices[j] = j;
                overlapMax[j] = overlaps.col(j).maxCoeff(&overlapIndices[j]);
                heuristicMax[j] = heuristics.col(j).maxCoeff();
            }

            constexpr double threshold = 0.70710678118;
            auto it = stateIndices.begin();
            std::sort(it, stateIndices.end(), [&](int i1, int i2) {
                return overlapMax[i1] > overlapMax[i2];
            });
            for (; it < stateIndices.end() && overlapMax[*it] > threshold; it++)
                ;
            std::sort(it, stateIndices.end(), [&](int i1, int i2) {
                return heuristicMax[i1] > heuristicMax[i2];
            });

            // do the actual matching procedure
            std::set<int> matched;
            for (int j : stateIndices) {
                auto overlap = overlaps.col(j);
                auto heuristic = heuristics.col(j);
                idx = overlapIndices[j];

                // index is already taken? Find greates overlap state that is
                // not already taken
                if (overlapMax[j] <= threshold ||
                    matched.find(idx) != matched.end()) {
                    std::sort(indices.begin(), indices.end(),
                              [&](int i1, int i2) {
                                  return heuristic[i1] < heuristic[i2];
                              });
                    int k = indices.size();
                    while (matched.find(indices[--k]) != matched.end() &&
                           k >= 0)
                        ;
                    idx = indices[k];
                }

                matched.insert(idx);
                energiesOrdered[j] = energies[idx];
                statesOrdered.col(j) = states.col(idx);
                characterOrdered.row(j) = character.row(idx);
            }

            // store processed data back to file
            m_ioThread.StoreData(i, ef, energiesOrdered, statesOrdered,
                                 characterOrdered);

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

int main(int argc, const char *argv[]) {
    // parsing command line arguments
    ArgumentParser argparse;

    std::string filename = GetDefaultAppDataDir("NO_StarkMap") + '/' +
                           GenerateFilename("NOStarkMap") + ".h5";
    argparse.AddOptionDefault("file", "HDF5 file path", filename);
    argparse.AddOption("help", "Print this help string");
    auto args = argparse.Parse(argc, argv);

    if (args.IsError()) {
        std::cout << args.GetError() << std::endl;
        return 1;
    } else if (args.IsOptionPresent("help")) {
        std::cout << argparse.GetHelpString() << std::endl;
        return 0;
    }
    filename = args.GetOptionStringValue("file");

    constexpr double rcm = EnergyInverseCm_v;
    constexpr double GHz = EnergyGHz_v;

    // parameters
    int n = 75;
    int R = 3;
    int dR = 2;
    double dE = 80 * rcm;

    // find state
    NitricOxide molecule;
    double energy = molecule.GetEnergy(RydbergDiatomicState_t(n, 10, R, 10, 0));
    std::cout << "State energy: " << energy / EnergyGHz_v << "GHz" << std::endl;

    // Run calculation
    NOStarkMapApp app(filename);
    // VectorXd eFields = VectorXd::LinSpaced(256, 0.0, 7.0); // V cm^-1
    VectorXd eFields = VectorXd::LinSpaced(1, 0.0, 0.0); // V cm^-1
    app.RunCalculation(100.0 * eFields, energy, dE, std::max(0, R - dR), R + dR,
                       0, 1, 100);

    return 0;
}
