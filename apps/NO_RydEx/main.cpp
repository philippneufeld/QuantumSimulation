// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/Util/SimulationApp.h>
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

class NOGasSensor
{
    constexpr static double xRadius = 115e-12; // bond length NO ground state
    constexpr static double rydRadius = 1e-9; // estimation of rydberg radius
public:

    NOGasSensor() 
        : m_beamRadius(1e-3), 
        m_mass(30.0061 * AtomicMassUnit_v), 
        m_temperature(300), 
        m_pressure(10.0) // 1e-3 mbar = 100Pa
    {
        m_mass = 30.0061 * AtomicMassUnit_v;
        
        // setup system levels
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, m_system.GetLevel(0) + SpeedOfLight_v / 226.97e-9);
        m_system.SetLevel(2, m_system.GetLevel(1) + SpeedOfLight_v / 540.0e-9);
        m_system.SetLevel(3, m_system.GetLevel(2) + SpeedOfLight_v / 834.92e-9);
        m_system.SetLevel(4, 9.27 * ElementaryCharge_v / PlanckConstant_v);

        // add dipole matrix elements
        m_system.SetDipoleElement(0, 1, 0.1595 * Debye_v); // https://doi.org/10.1093/mnras/stx1211
        m_system.SetDipoleElement(1, 2, 1e-2 * Debye_v);
        m_system.SetDipoleElement(2, 3, 1e-1 * Debye_v);

        // add coupling lasers
        double uvInt = NLevelLaser::PowerToIntensity(0.05, m_beamRadius); // 50mW
        double greenInt = NLevelLaser::PowerToIntensity(1.0, m_beamRadius); // 1W
        double redInt = NLevelLaser::PowerToIntensity(0.5, m_beamRadius); // 500mW
        m_system.AddLaser(NLevelLaser({0, 1}, uvInt, 1.0));
        m_system.AddLaser(NLevelLaser({1, 2}, greenInt, -1.0));
        m_system.AddLaser(NLevelLaser({2, 3}, redInt, -1.0));

        // set level decays
        m_system.AddDecay(1, 0, 13.8e6); // https://doi.org/10.1063/1.454958
        m_system.AddDecay(2, 1, 1.0e6);
        m_system.AddDecay(3, 2, 0.5e6);

        // add transit broadening effect
        double meanVel = std::sqrt(2*BoltzmannConstant_v*300 / m_mass);
        double decayTransit = meanVel / m_beamRadius;
        m_system.AddDecay(1, 0, decayTransit);
        m_system.AddDecay(2, 0, decayTransit);
        m_system.AddDecay(3, 0, decayTransit);
        m_system.AddDecay(4, 0, decayTransit);

        // add collisional ionization term
        // collisions between rydberg and ground state
        double rateRX = GetCollisionRate(rydRadius + xRadius);
        m_system.AddDecay(3, 4, rateRX);

        // add recombination rate (assume to be the same as ionization term for now)
        m_system.AddDecay(4, 0, rateRX);

        // initialize helper variables
        m_transitions = m_system.GetTransitionFreqs();
        m_laserDirs = m_system.GetLaserDirs();
    }

    double GetCollisionRate(double rEff) const
    {
        double kT = BoltzmannConstant_v * m_temperature;
        double n = m_pressure / kT; // number density
        double sigma = Pi_v * std::pow(rEff, 2);
        double mtp = 1.0 / (std::sqrt(2) * sigma * n); // mean free path
        double vExp = std::sqrt(8*kT / (Pi_v*m_mass)); // expected velocity
        return vExp / mtp;
    }

    Matrix<double, 5, 1> GetPopulations(double uvDet, double greenDet, double redDet) const
    {
        // make a local copy of the system (makes concurrent calls to this function possible)
        auto sys = m_system;

        // calculate ionization rate (rate0 + rate1*rho_ion)
        double rate0 = GetCollisionRate(rydRadius + xRadius); // constant term
        double rate1 = GetCollisionRate(rydRadius);

        // starting point of the calculation
        auto rho0 = sys.CreateGroundState();
        
        TDopplerIntegrator<> doppler(m_mass, m_temperature);
        auto dopplerCalc = [&](double vel)
        {
            // calculate laser frequencies
            Vector3d detunings(uvDet, greenDet, redDet);
            Vector3d laserFreqs = m_transitions + detunings;
            
            // adjust laser frequencies for doppler shift
            laserFreqs = doppler.ShiftFrequencies(laserFreqs, m_laserDirs, vel);
            
            // solve Lindblad master equation for the steady state
            auto rho = rho0;
            auto rhoPrev = rho;
            for (int i = 0; i < 10; i++)
            {
                double rhoIon = std::real(rho(4, 4));
                sys.SetDecay(3, 4, rate0 + rate1*rhoIon);
                
                rhoPrev = rho;
                rho = sys.GetDensityMatrixSS(laserFreqs);

                if (rho.isApprox(rhoPrev))
                    break;
            }

            // return populations of the levels
            return real(rho.diagonal().array()).matrix().eval();
        };

        return doppler.Integrate(dopplerCalc);
    }

private:
    TNLevelSystemQM<5> m_system;

    // environmental parameters
    double m_temperature;
    double m_mass;
    double m_pressure;

    // laser variables
    double m_beamRadius;
    Vector3d m_transitions;
    Vector3d m_laserDirs;
};


int main(int argc, const char* argv[])
{
    NOGasSensor gasSensor;

    VectorXd detunings = VectorXd::LinSpaced(501, -1e9, 1e9);
    Matrix<double, 5, Dynamic> populations(5, detunings.size());

    // start parallelized calculation
    ThreadPool pool; 
    ProgressBar progress(detunings.size());
    for (std::size_t i = 0; i < detunings.size(); i++)
    {
        pool.Submit([&, i=i]{ 
            populations.col(i) = gasSensor.GetPopulations(0.0, 0.0, detunings[i]);
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();


#ifdef QSIM_PYTHON3
        PythonMatplotlib matplotlib;
        std::string levelNames[] = {"X", "A", "H", "R", "Ion"};

        // for (std::size_t i = 0; i < 5; i++)
        std::size_t i = 4;
        {
            VectorXd pops = populations.row(i);
            auto figure = matplotlib.CreateFigure();
            auto ax = figure.AddSubplot();
            ax.SetXLabel("Red laser detuning [GHz]");
            ax.SetYLabel(levelNames[i] + " population");
            ax.Plot((detunings / 1e9).eval().data(), pops.data(), detunings.size());
        }

        matplotlib.RunGUILoop();
#endif

}
