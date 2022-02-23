// Philipp Neufeld, 2021-2022

#include <iostream>
#include <random>
#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/NLevelSystemQM.h>
#include <QSim/NLevel/NLevelSystemSC.h>
#include <QSim/NLevel/Doppler.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
    double intProbe = NLevelLaser::RabiToIntensity(dip, 30.5e6);

    double modPeriod = 1e-6;
    ModulatedNLevelLaser laser({0, 1});
    laser.SetIntensity(intProbe);
    laser.SetModulationFunc([=](double t){ t = std::fmod(t, modPeriod); return t <= modPeriod / 2 ? 1.0 : 0.0; });

    TNLevelSystemSC<DynamicDim_v, true> system2(2);
    system2.SetLevel(0, 0.0);
    system2.SetLevel(1, SpeedOfLight_v / 780.241e-9);
    system2.SetDecay(1, 0, 6.065e6);
    system2.SetDipoleElement(0, 1, dip);
    system2.AddLaser(laser);

    double dt = 2e-10;
    auto rho0 = system2.CreateGroundState();
    
    TDopplerIntegrator<QuadMidpointPolicy> doppler;
    doppler.SetIntegrationSteps(300);
    doppler.SetMass(1.44316060e-25);
    doppler.SetTemperature(0);

    auto traj = doppler.Integrate([&](double vel)
    { 
        auto rhos = system2.GetTrajectory(VectorXd::Zero(1), rho0, vel, 0.0, 2e-6, dt);
        VectorXd pops(rhos.size());
        for (int i=0; i<rhos.size(); i++) 
            pops[i] = std::real(rhos[i](1, 1));
        return pops;
    });

    auto ts = system2.GetTrajectoryTimeaxis(0.0, dt, traj.size());


#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(ts.data(), traj.data(), ts.size());
    matplotlib.RunGUILoop();
#endif

    return 0;
}
