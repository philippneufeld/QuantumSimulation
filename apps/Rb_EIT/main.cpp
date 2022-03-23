// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

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
    constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
    double intProbe = NLevelLaser::RabiToIntensity(dip, 3.5e6);
    double intPump = NLevelLaser::RabiToIntensity(dip, 10.0e6);

    // Create system
    TNLevelSystemQM<3> system;
    system.SetLevel(0, -4.271e9);
    system.SetLevel(1, 2.563e9);
    system.SetLevel(2, SpeedOfLight_v / 780.241e-9);
    system.SetDipoleElement(0, 2, dip);
    system.SetDipoleElement(1, 2, dip);
    system.AddLaser(NLevelLaser({0, 2}, intProbe, 1.0));
    system.AddLaser(NLevelLaser({1, 2}, intPump, 1.0));
    system.SetDecay(2, 0, 3.0/8.0 * 6.065e6);
    system.SetDecay(2, 1, 5.0/8.0 * 6.065e6);
    
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);

    // Generate detuning axis
    VectorXd detunings = VectorXd::LinSpaced(501, -125e6, 125e6);
    VectorXd absCoeffs(detunings.size());

    // get properties of the system and the lasers
    auto transitions = system.GetTransitionFreqs();
    auto dirs = system.GetLaserDirs();

    // start calculation
    ThreadPool pool; 
    ProgressBar progress(detunings.size());
    for (std::size_t i = 0; i < detunings.size(); i++)
    {
        pool.Submit([&, i=i](){ 
            absCoeffs[i] = doppler.Integrate([&](double vel)
            {
                VectorXd laserFreqs = transitions + Vector2d(detunings[i], 0.0);
                laserFreqs = doppler.ShiftFrequencies(laserFreqs, dirs, vel);

                auto rho = system.GetDensityMatrixSS(laserFreqs);
                return std::imag(rho(0, 2));
            });
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

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
