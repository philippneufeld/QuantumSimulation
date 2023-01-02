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
#include <QSim/Execution/ServerPool.h>

#include <QSim/Util/Argparse.h>
#include <QSim/Util/ProgressBar.h>
#include <QSim/Util/PathUtil.h>

#include "IOThread.h"

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

class NOStarkMapAppWorker : public ServerPoolWorker
{
public:
    NOStarkMapAppWorker() : ServerPoolWorker(GetNumberOfCalcThreads()) { }

    virtual DataPackagePayload DoWork(DataPackagePayload data) override
    {
        std::cout << "DoWork()" << std::endl;

        std::uint32_t n = *reinterpret_cast<std::uint32_t*>(data.GetData());
        Eigen::MatrixXd hamiltonian = Eigen::Map<Eigen::MatrixXd>(
            reinterpret_cast<double*>(data.GetData() + 4), n, n);

        auto [energies, states] = DiatomicStarkMap::GetEnergiesAndStatesFromHamiltonian(hamiltonian);

        DataPackagePayload result(4 + (n+n*n)*sizeof(double));
        *reinterpret_cast<std::uint32_t*>(result.GetData()) = n;
        double* dres = reinterpret_cast<double*>(result.GetData() + 4);
        Eigen::Map<Eigen::VectorXd>(dres, n) = energies;
        Eigen::Map<Eigen::MatrixXd>(dres+n, n, n) = states;

        return result;
    }

};

class NOStarkMapApp : public ServerPool
{
  public:
    NOStarkMapApp(const std::string path) : m_path(path), m_ioThread(path) {}

    void RunCalculation(const VectorXd &eFields, double energy, double dE, const std::vector<int>& Rs, int mN, int nMax) 
    {
        if (GetWorkerCount() == 0)
        {
            std::cout << "No worker servers connected!" << std::endl;
            return;
        }

        // initialize calculation and process basis
        std::cout << "Searching for appropriate basis set..." << std::endl;
        DiatomicStarkMap starkMap(std::make_unique<NitricOxide>(), Rs, nMax, mN, energy, dE);
        std::vector<RydbergDiatomicState_t> basis = starkMap.GetBasis();
        AnalyzeBasis(basis);
        std::cout << "Found basis set of size: " << basis.size() << std::endl;

        std::cout << "Calculating hamiltonian..." << std::endl;
        {
            ThreadPool tp;
            starkMap.PrepareCalculation(tp);
        }
        std::cout << "Calculation of hamiltonian finished. Starting diagonalization..." << std::endl;

        // start i/o thread
        m_ioThread.Start(basis, energy, dE);

        // do the actual calculation
        ProgressBar progress(eFields.size());
        std::thread serverThread([&](){ Run(); });
        for (int i = 0; i < eFields.size(); i++) 
        {
            auto [task, fut] = Submit([&, ef=eFields[i]]()
            {
                std::cout << "Submit " << std::to_string(ef) << std::endl;

                auto n = basis.size();
                DataPackagePayload data(4+n*n*sizeof(double));
                *reinterpret_cast<std::uint32_t*>(data.GetData()) = n;
                Eigen::Map<Eigen::MatrixXd> rawHamiltonian(
                    reinterpret_cast<double*>(data.GetData() + 4), n, n);
                rawHamiltonian = starkMap.GetHamiltonian(ef);
                return data;
            });
            m_tasks[task] = std::make_tuple(i, eFields[i], std::move(fut), &progress);
        }

        progress.WaitUntilFinished();

        // stop thread
        WaitUntilFinished();
        Stop();
        serverThread.join();

        std::cout << "Waiting for I/O to complete..." << std::endl;
        m_ioThread.WaitUntilFinished();

        // do post processing
        std::cout << "Overlapping states..." << std::endl;
        MatchStates(eFields, basis.size());

        std::cout << "Data stored successfully (" << m_path << ")" << std::endl;
    }

    virtual void OnTaskCompleted(UUIDv4 id) override
    {
        auto it = m_tasks.find(id);
        std::cout << "Task completed " << (it != m_tasks.end()) << std::endl;
        if (it != m_tasks.end())
        {
            auto [i, ef, fut, prog] = std::move(it->second);
            DataPackagePayload data = std::move(fut.get());
            std::uint32_t n = *reinterpret_cast<std::uint32_t*>(data.GetData());
            double* dptr = reinterpret_cast<double*>(data.GetData() + 4);
            Eigen::VectorXd energies = Eigen::Map<Eigen::VectorXd>(dptr, n);
            Eigen::MatrixXd states = Eigen::Map<Eigen::MatrixXd>(dptr+n, n, n);

            auto character = DetermineStateCharacter(states);

            m_ioThread.StoreData(i, ef, energies, states, character);
            m_tasks.erase(it);

            prog->IncrementCount();
        }
    }

  protected:
    void AnalyzeBasis(const std::vector<RydbergDiatomicState_t> &basis) 
    {
        int basisSize = basis.size();
        m_minn = std::numeric_limits<int>::max();
        int maxn = std::numeric_limits<int>::min();
        int maxR = std::numeric_limits<int>::min();
        int maxL = std::numeric_limits<int>::min();
        int maxN = std::numeric_limits<int>::min();
        for (const auto &state : basis) 
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
            m_nCharMatrix(n - m_minn, i) = 1.0;
            m_lCharMatrix(l, i) = 1.0;
            m_RCharMatrix(R, i) = 1.0;
            m_NCharMatrix(N, i) = 1.0;
        }
    }

    Matrix<double, Dynamic, 4> DetermineStateCharacter(const MatrixXd &states) 
    {
        Matrix<double, Dynamic, 4> character(states.rows(), 4);
        auto statesSq = states.cwiseAbs2().eval();
        MatrixXd nMat = m_nCharMatrix * statesSq;
        MatrixXd lMat = m_lCharMatrix * statesSq;
        MatrixXd RMat = m_RCharMatrix * statesSq;
        MatrixXd NMat = m_NCharMatrix * statesSq;

        for (int j = 0; j < states.rows(); j++) 
        {
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

    void MatchStates(const VectorXd &eFields, int basisSize) 
    {
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
        for (int i = 0; i < eFields.size(); i++) 
        {
            auto [idx, ef, energies, states, character] = nextDataPromise.get();
            if (i < eFields.size() - 1)
                nextDataPromise = m_ioThread.LoadData(i + 1);

            // prepare ordering auxilliary variables
            if (i == 0) 
            {
                prevEnergies = energies;
                prevEnergies2 = prevEnergies;
                prevStates = states;
                continue;
            }

            // calculate heuristic matching criteria
            overlaps = (states.transpose() * prevStates).cwiseAbs();
            // energyDiffs(overlaps.rows(), overlaps.cols());
            for (int i1 = 0; i1 < energyDiffs.rows(); i1++) 
            {
                for (int i2 = 0; i2 < energyDiffs.rows(); i2++) 
                {
                    extrEnergies = prevEnergies2;
                    if (i > 1)
                        extrEnergies += (prevEnergies - prevEnergies2) / (eFields[i - 1] - eFields[i - 2]) *
                                        (eFields[i] - eFields[i - 2]);
                    energyDiffs(i1, i2) = std::abs(energies[i1] - extrEnergies[i2]);
                }
            }
            energyDiffs /= energyDiffs.maxCoeff();
            heuristics = overlaps - energyDiffs.cwiseSqrt();

            // determine order in which states are processed
            for (int j = 0; j < states.rows(); j++) 
            {
                stateIndices[j] = j;
                overlapMax[j] = overlaps.col(j).maxCoeff(&overlapIndices[j]);
                heuristicMax[j] = heuristics.col(j).maxCoeff();
            }

            constexpr double threshold = 0.70710678118;
            auto it = stateIndices.begin();
            std::sort(it, stateIndices.end(), [&](int i1, int i2) { return overlapMax[i1] > overlapMax[i2]; });
            for (; it < stateIndices.end() && overlapMax[*it] > threshold; it++);
            std::sort(it, stateIndices.end(), [&](int i1, int i2) { return heuristicMax[i1] > heuristicMax[i2]; });

            // do the actual matching procedure
            std::set<int> matched;
            for (int j : stateIndices) 
            {
                auto overlap = overlaps.col(j);
                auto heuristic = heuristics.col(j);
                idx = overlapIndices[j];

                // index is already taken? Find greates overlap state that is
                // not already taken
                if (overlapMax[j] <= threshold || matched.find(idx) != matched.end()) {
                    std::sort(indices.begin(), indices.end(),
                              [&](int i1, int i2) { return heuristic[i1] < heuristic[i2]; });
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
    std::map<UUIDv4, std::tuple<
        int, double, 
        std::future<DataPackagePayload>, 
        ProgressBar*>> m_tasks;

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

    std::string filename = GetDefaultAppDataDir("NO_StarkMap", false) + '/' + GenerateFilename("NOStarkMap") + ".h5";
    argparse.AddOption("worker", "Start worker thread");
    argparse.AddOptionDefault("file", "HDF5 file path", filename);
    argparse.AddOptionDefault("n", "Principal quantum number", "42");
    argparse.AddOptionDefault("l", "Orbital angular momentum quantum number", "3");
    argparse.AddOptionDefault("R", "Rotational quantum number", "2");
    argparse.AddOptionDefault("N", "Total angular momentum quantum number", "3");
    argparse.AddOptionDefault("mN", "Projection total angular momentum quantum number", "0");
    argparse.AddOptionDefault("nMax", "Principal quantum number basis maximum", "100");
    argparse.AddOptionDefault("Rs", "Rotational quantum number basis", "0,1,2,3");
    argparse.AddOptionDefault("Fmin", "Rotational quantum number (V/cm)", "0.0");
    argparse.AddOptionDefault("Fmax", "Rotational quantum number (V/cm)", "7.0");
    argparse.AddOptionDefault("Fsteps", "Rotational quantum number", "256");
    argparse.AddOption("dErcm", "Energy range cm^-1");
    argparse.AddOption("dEGHz", "Energy range GHz");
    argparse.AddOption("help", "Print this help string");
    auto args = argparse.Parse(argc, argv);

    if (args.IsError()) {
        std::cout << args.GetError() << std::endl;
        return 1;
    } 
    else if (args.IsOptionPresent("dErcm") && args.IsOptionPresent("dEGHz"))
    {
        std::cout << "dErcm and dEGHz option are mutually exclusive" << std::endl;
        return 1;
    }
    else if (args.IsOptionPresent("help")) {
        std::cout << argparse.GetHelpString() << std::endl;
        return 0;
    }


    short port = 8000;
    if (args.IsOptionPresent("worker"))
    {
        NOStarkMapAppWorker worker;
        return worker.Run(port);
    }
    else
    {
        filename = args.GetOptionStringValue("file");

        std::string dir = std::filesystem::path(filename).parent_path().string();
        std::filesystem::create_directories(dir);

        // retrieve parameter
        int n = args.GetOptionValue<int>("n");
        int l = args.GetOptionValue<int>("l");
        int R = args.GetOptionValue<int>("R");
        int N = args.GetOptionValue<int>("N");
        int mN = args.GetOptionValue<int>("mN");
        int nMax = args.GetOptionValue<int>("nMax");
        std::vector<int> Rs = args.GetOptionValue<std::vector<int>>("Rs");
        double Fmin = args.GetOptionValue<double>("Fmin");
        double Fmax = args.GetOptionValue<double>("Fmax");
        int Fsteps = args.GetOptionValue<int>("Fsteps");

        double dE = 10 * EnergyInverseCm_v;
        if (args.IsOptionPresent("dErcm"))
            dE = args.GetOptionValue<double>("dErcm") * EnergyInverseCm_v;
        else if (args.IsOptionPresent("dEGHz"))
            dE = args.GetOptionValue<double>("dEGHz") * EnergyGHz_v;

        // validate quantum numbers
        auto state = RydbergDiatomicState_t(n, l, R, N, mN);
        if (n < 0 || l < 0 || l >= n || R < 0 || N < std::abs(R-l) || N > std::abs(R+l) || std::abs(mN) > N)
        {
            std::cout << "Invalid quantum numbers" << std::endl;
            return 0;
        }
        
        // find state
        NitricOxide molecule;
        double energy = molecule.GetEnergy(state);
        std::cout << "State energy: " << energy / EnergyGHz_v << "GHz / " << energy / EnergyInverseCm_v << "cm^{-1}" << std::endl;

        // Run calculation
        NOStarkMapApp app(filename);

        app.ConnectWorkerHostname("ludwigsburg", port);
        app.ConnectWorkerHostname("calca", port);
        app.ConnectWorkerHostname("calcb", port);
        app.ConnectWorkerHostname("calcc", port);
        app.ConnectWorkerHostname("calcv", port);
        app.ConnectWorkerHostname("calcr", port);
        app.ConnectWorkerHostname("sost", port);

        VectorXd eFields = VectorXd::LinSpaced(Fsteps, Fmin, Fmax); // V cm^-1
        app.RunCalculation(100.0 * eFields, energy, dE, Rs, mN, nMax);
    }

    return 0;
}
