// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Util/DataFile.h>

using namespace QSim;
using namespace Eigen;

std::pair<VectorXd, VectorXd> CalculateTrajectory(
    double f0, double detuning, double rabi, double gamma, double tmax)
{
    TNLevelSystem<2> system;
    system.SetLevel(0, 0);
    system.SetLevel(1, f0);
    system.AddCoupling(0, 1, rabi);
    system.AddDecay(1, 0, gamma);

    // get properties of the system and the lasers
    VectorXd laserFreqs = system.GetCouplingResonanceFreqs().array() + detuning;

    // simulate trajectory
    double dt = 1.0e-2 / std::max(rabi, gamma);
    auto rhos = system.GetTrajectory(laserFreqs, system.CreateGroundState(), 0.0, tmax, dt);
    auto ts = system.GetTrajectoryTimeaxis(0.0, dt, rhos.size());
    
    VectorXd pops(rhos.size());
    for (int i=0; i<rhos.size(); i++) 
        pops[i] = std::real(rhos[i](1, 1));

    return std::make_pair(ts, pops);
}

int main(int argc, const char* argv[])
{
    double rabi = 10.0e6;

    std::map<std::string, std::pair<double, double>> sims;
    sims["R0G0"] = {0*rabi, 0.0};
    sims["R1G0"] = {1*rabi, 0.0};
    sims["R2G0"] = {2*rabi, 0.0};
    sims["R0G.1"] = {0*rabi, 0.1*rabi};
    sims["R1G.1"] = {1*rabi, 0.1*rabi};
    sims["R2G.1"] = {2*rabi, 0.1*rabi};

    DataFile file;
    file.Open("RabiOscillations.h5", DataFile_TRUNCATE);
    auto root = file.OpenRootGroup();

    for (auto it=sims.begin(); it!=sims.end(); it++)
    {
        auto [det, gamma] = it->second;
        auto grp = root.CreateSubgroup(it->first);

        auto [ts, pops] = CalculateTrajectory(1e15, det, rabi, gamma, 25.0 / rabi);
        grp.CreateDataset("t", (ts * rabi * TwoPi_v).eval());
        grp.CreateDataset("pop", pops);
    }
        
    file.Close();

    return 0;
}
