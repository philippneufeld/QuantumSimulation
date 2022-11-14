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
    constexpr double rabiProbe = 3.5e6;
    constexpr double rabiPump = 10.0e6;

    // Create system
    TNLevelSystem<3> system;
    system.SetLevel(0, -4.271e9);
    system.SetLevel(1, 2.563e9);
    system.SetLevel(2, SpeedOfLight_v / 780.241e-9);
    system.AddCoupling(0, 2, rabiProbe);
    system.AddCoupling(1, 2, rabiPump);
    system.SetDecay(2, 0, 3.0/8.0 * 6.065e6);
    system.SetDecay(2, 1, 5.0/8.0 * 6.065e6);
    
    // create doppler integartor
    TDopplerIntegrator<> doppler;
    doppler.SetMass(1.44316060e-25);
    VectorXd resFreqs = system.GetCouplingResonanceFreqs();
    VectorXd laserDirs = VectorXd::Ones(2);

    // Generate detuning axis
    VectorXd detunings = VectorXd::LinSpaced(501, -125e6, 125e6);
    VectorXd absCoeffs(detunings.size());

    // start calculation
    ThreadPool pool; 
    ProgressBar progress(detunings.size());
    for (std::size_t i = 0; i < detunings.size(); i++)
    {
        pool.Submit([&, i=i](){ 
            absCoeffs[i] = doppler.Integrate([&](double vel)
            {
                VectorXd laserFreqs = resFreqs + Vector2d(detunings[i], 0.0);
                laserFreqs = doppler.ShiftFrequencies(laserFreqs, laserDirs, vel);

                auto rho = system.GetDensityMatrixSS(laserFreqs);
                return std::imag(rho(0, 2));
            });
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

    // TODO: Save data

    return 0;
}
