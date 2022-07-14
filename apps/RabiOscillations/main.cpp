// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    double modPeriod = 1e-6;
    auto rabiFunc = [=](double t){ t = std::fmod(t, modPeriod); return 30.5e6 * (t <= modPeriod / 2 ? 1.0 : 0.0); }; 

    TNLevelSystem<DynamicDim_v, true> system(2);
    system.SetLevel(0, 0.0);
    system.SetLevel(1, SpeedOfLight_v / 780.241e-9);
    system.SetDecay(1, 0, 6.065e6);
    system.AddCoupling(0, 1, rabiFunc);

    double dt = 2e-10;
    auto rho0 = system.CreateGroundState();
    
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);
    doppler.SetTemperature(0);

    // get properties of the system and the lasers
    VectorXd resonanceFreqs = system.GetCouplingResonanceFreqs();
    VectorXd dirs = VectorXd::Ones(resonanceFreqs.size());

    auto traj = doppler.Integrate([&](double vel)
    { 
        VectorXd laserFreqs = doppler.ShiftFrequencies(resonanceFreqs, dirs, vel);
        auto rhos = system.GetTrajectory(laserFreqs, rho0, 0.0, 2e-6, dt);
        VectorXd pops(rhos.size());
        for (int i=0; i<rhos.size(); i++) 
            pops[i] = std::real(rhos[i](1, 1));
        return pops;
    });

    auto ts = system.GetTrajectoryTimeaxis(0.0, dt, traj.size());

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(ts.data(), traj.data(), ts.size());
    matplotlib.RunGUILoop();
#endif

    return 0;
}
