// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <fstream>

#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    // calculate parameters
    constexpr double decay = 6.065e6;
    constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
    double intProbe = NLevelLaser::RabiToIntensity(dip, 3.5e6);
    double intPump = NLevelLaser::RabiToIntensity(dip, 10.0e6);
    constexpr double freq = SpeedOfLight_v / 780.241e-9;

    // dt << Rabi^-1, Doppler^-1, detuning^-1
    constexpr double decayTime = 1.0 / decay;
    constexpr double tint = 10.0 * decayTime;
    constexpr double dt = decayTime / 250.0;
    
    // create system
    TNLevelSystemSC<2> system;
    system.SetLevel(0, 0.0);
    system.SetLevel(1, freq);
    system.SetDipoleElement(0, 1, dip);
    system.AddLaser(NLevelLaser({0, 1}, intProbe, 1.0));
    system.AddLaser(NLevelLaser({0, 1}, intPump, -1.0));
    system.SetDecay(1, 0, decay);
    
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);

    // Generate detuning axis
    VectorXd detunings = VectorXd::LinSpaced(501, -1e9, 1e9);
    // VectorXd detunings = VectorXd::LinSpaced(51, -10.5e7, 10.5e7);
    VectorXd absCoeffs(detunings.size());

    // get properties of the system and the lasers
    auto transitions = system.GetTransitionFreqs();
    auto dirs = system.GetLaserDirs();

    // start calculation
    ThreadPool pool; 
    ProgressBar progress(detunings.size());
    auto rho0 = system.CreateGroundState();
    for (std::size_t i = 0; i < detunings.size(); i++)
    {
        pool.Submit([&, i=i](){ 
            absCoeffs[i] = doppler.Integrate([&](double vel)
            {
                VectorXd laserFreqs = transitions + Vector2d(0.0, detunings[i]);
                laserFreqs = doppler.ShiftFrequencies(laserFreqs, dirs, vel);

                auto rho = system.GetDensityMatrixAv(laserFreqs, rho0, 0.0, tint, 0.25*tint, dt);
                return std::real(rho(1, 1));
            });
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

    std::ofstream file;
    std::string home = getenv("HOME");
    file.open(home + "/data.txt");
    for (int i = 0; i < detunings.size(); i++)
    {
        file << detunings[i] << " " << absCoeffs[i] << std::endl;
    }
    file.close();


    // plot data
#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(detunings.data(), absCoeffs.data(), detunings.size());
    matplotlib.RunGUILoop();
#endif

    return 0;
}
