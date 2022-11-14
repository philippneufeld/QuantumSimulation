// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    constexpr double rabi = 3.5e6;
    constexpr double freq = SpeedOfLight_v / 780.241e-9;

    // create system
    TNLevelSystem<2> system;
    system.SetLevel(0, 0.0);
    system.SetLevel(1, freq);
    system.AddCoupling(0, 1, rabi);
    system.SetDecay(1, 0, 6.065e6);
    
    // create doppler integartor
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);
    VectorXd resFreqs = system.GetCouplingResonanceFreqs();
    VectorXd laserDirs = VectorXd::Ones(1);

    // Generate detuning axis
    VectorXd detunings = VectorXd::LinSpaced(501, -1e9, 1e9);
    VectorXd absCoeffs(detunings.size());

    // start calculation
    ThreadPool pool; 
    ProgressBar progress(detunings.size());
    for (std::size_t i = 0; i < detunings.size(); i++)
    {
        pool.Submit([&, i=i](){ 
            absCoeffs[i] = doppler.Integrate([&](double vel)
            {
                VectorXd laserFreqs = resFreqs + Matrix<double, 1, 1>(detunings[i]);
                laserFreqs = doppler.ShiftFrequencies(laserFreqs, laserDirs, vel);
                auto rho = system.GetDensityMatrixSS(laserFreqs);
                return std::imag(rho(0, 1));
            });
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

    // TODO: Save data

    return 0;
}
