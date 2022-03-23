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
    constexpr double freq = SpeedOfLight_v / 780.241e-9;

    // create system
    TNLevelSystemQM<2> system;
    system.SetLevel(0, 0.0);
    system.SetLevel(1, freq);
    system.SetDipoleElement(0, 1, dip);
    system.AddLaser(NLevelLaser({0, 1}, intProbe, 1.0));
    system.SetDecay(1, 0, 6.065e6);
    
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);

    // Generate detuning axis
    VectorXd detunings = VectorXd::LinSpaced(501, -1e9, 1e9);
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
                VectorXd laserFreqs = transitions + Matrix<double, 1, 1>(detunings[i]);
                laserFreqs = doppler.ShiftFrequencies(laserFreqs, dirs, vel);

                auto rho = system.GetDensityMatrixSS(laserFreqs);
                return std::imag(rho(0, 1));
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
