// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/Util/SimulationApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>
#include <QSim/Util/PathUtil.h>

#include <mutex>


#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <iostream>

using namespace QSim;
using namespace Eigen;

class NOGasSensorTD
{
    constexpr static double xRadius = 115e-12; // bond length NO ground state
    constexpr static double rydRadius = 1e-9; // estimation of rydberg radius
public:

    NOGasSensorTD() 
        : m_beamRadius(1e-3), 
        m_mass(30.0061 * AtomicMassUnit_v), 
        m_temperature(300), 
        m_pressure(100.0) // 1e-3 mbar = 100Pa
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
        m_system.SetDipoleElement(1, 2, 1e-1 * Debye_v);
        m_system.SetDipoleElement(2, 3, 5e-1 * Debye_v);

        // add coupling lasers
        using TLaser = ModulatedNLevelLaser;
        double uvInt = TLaser::PowerToIntensity(0.05, m_beamRadius); // 50mW
        double greenInt = TLaser::PowerToIntensity(1.0, m_beamRadius); // 1W
        double redInt = TLaser::PowerToIntensity(1.0, m_beamRadius); // 1W
        m_system.AddLaser(TLaser({0, 1}, uvInt, 1.0));
        m_system.AddLaser(TLaser({1, 2}, greenInt, -1.0));
        m_system.AddLaser(TLaser({2, 3}, redInt, -1.0));

        // set level decays
        m_system.AddDecay(1, 0, 13.8e6); // https://doi.org/10.1063/1.454958
        m_system.AddDecay(2, 1, 10.0e6);
        m_system.AddDecay(3, 2, 10.0e6);

        // add transit broadening effect
        double meanVel = std::sqrt(2*BoltzmannConstant_v*300 / m_mass);
        double decayTransit = meanVel / m_beamRadius;
        m_system.AddDecay(1, 0, decayTransit);
        m_system.AddDecay(2, 0, decayTransit);
        m_system.AddDecay(3, 0, decayTransit);
        m_system.AddDecay(4, 0, decayTransit);

        // add collisional ionization term
        // collisions between rydberg and ground state
        double sigmaRX = Pi_v * (rydRadius + xRadius) * (rydRadius + xRadius);
        double kT = BoltzmannConstant_v * m_temperature;
        double n = m_pressure / kT; // number density
        double rate = n * std::sqrt(8*kT / Pi_v) * sigmaRX * std::sqrt(2 / m_mass);
        m_system.AddDecay(3, 4, rate);

        // add recombination rate (assume to be the same as ionization term for now)
        m_system.AddDecay(4, 0, rate);

        // initialize helper variables
        m_transitions = m_system.GetTransitionFreqs();
        m_laserDirs = m_system.GetLaserDirs();
    }

    auto GetPopulationsTrajectory(
        double uvDet, double greenDet, double redDet, 
        double t, double dt, double modFreq) const
    {
        // make a local copy of the system (makes concurrent calls to this function possible)
        auto sys = m_system;

        // modulate green laser
        modFreq = std::abs(modFreq);
        if (modFreq == 0)
        {
            sys.GetLaser(1).SetModulationFunc([=](double t){ return 1.0; });
        }
        else
        {
            double modPeriod = 1.0 / modFreq;
            auto modFunc = [=](double t){ t = std::fmod(t, modPeriod); return t <= modPeriod / 2 ? 1.0 : 0.0; };
            sys.GetLaser(1).SetModulationFunc(modFunc);
        }     

        // geometric collisional cross-sections
        double sigmaRe = Pi_v * rydRadius * rydRadius;
        
        // calculate ionization rate (rate0 + rate1*rho_ion)
        double kT = BoltzmannConstant_v * m_temperature;
        double n = m_pressure / kT; // number density
        double rate0 = sys.GetDecay(3, 4); // constant term
        double rate1 = n * std::sqrt(8*kT / Pi_v) * sigmaRe / std::sqrt(ElectronMass_v);

        // starting point of the calculation
        auto rho0 = sys.CreateGroundState();

        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil(t / dt));
        dt = t / steps;

        // limit runtime by sacrificing some precision in the doppler integration
        TDopplerIntegrator<> doppler(m_mass, m_temperature);
        doppler.SetIntegrationRTol(1e-4);
        doppler.SetIntegrationWidth(1.0);

        auto dopplerCalc = [&](double vel)
        {
            VectorXd ionPops(steps);

            // calculate laser frequencies
            Vector3d detunings(uvDet, greenDet, redDet);
            Vector3d laserFreqs = m_transitions + detunings;
            
            // adjust laser frequencies for doppler shift
            laserFreqs = doppler.ShiftFrequencies(laserFreqs, m_laserDirs, vel);
            const auto auxData = sys.GetHamiltonianAux(laserFreqs);

            // define function to be integrated (non-linearity included)
            auto func = [&](double x, const decltype(rho0)& rho) 
            {
                // non-linear dependece on ion population
                double rhoIon = std::real(rho(4, 4));
                sys.SetDecay(3, 4, rate0 + rate1*rhoIon);

                return sys.GetDensityOpDerivativeFast(auxData, rho, x); 
            };

            // integrate trajectory
            // TODEIntegrator<ODERK4Policy> integrator;
            TODEIntegrator<ODEAd54DPPolicy> integrator;
            auto rho = rho0;
            double dtEff = dt; // dt that is controlled by the adaptive stepsize control
            for (int i = 0; i < steps; i++)
            {
                rho = integrator.IntegrateTo(func, rho, i*dt, (i+1)*dt, dtEff);
                ionPops(i) = std::real(rho(4, 4));
            }
            
            // return populations of the levels
            return ionPops;
        };

        // auto trajectory = dopplerCalc(100.0);
        auto trajectory = doppler.Integrate(dopplerCalc);

        VectorXd timeAxis = sys.GetTrajectoryTimeaxis(0.0, dt, trajectory.size());
        return std::make_pair(timeAxis, trajectory);
    }

private:
    TNLevelSystemSC<5, true> m_system;

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
    NOGasSensorTD gasSensor;

    double fmin = 1e6;
    double fmax = 1e8;
    VectorXd freqs = ArrayXd::LinSpaced(50, std::log(fmin), std::log(fmax)).exp();
    VectorXd populations = VectorXd::Zero(freqs.size());

    double dt = 1e-9;
    double tmax = 5.0 / fmin;

    // data storage
    // generate filename (first store locally and then move to desired location)
    std::string filename = GenerateFilename("NORydExTD") + ".h5";
    std::string dstPath = GetHomeDirSubfolderPath("remote_home") + "/Masterarbeit/07_TimeDependence/01_Ion_modG/" + filename;

    {
        DataFile file;
        file.Open(filename, DataFile_MUST_NOT_EXIST);
        auto root = file.OpenRootGroup();
        auto groupNameLen = std::to_string(freqs.size()).size();
        std::mutex file_mutex;

        ThreadPool threadPool;
        ProgressBar progressBar(freqs.size());

        for (int i = 0; i < freqs.size(); i++)
        {
            threadPool.Submit([&, i=i]()
            {
                auto [ts, pops] = gasSensor.GetPopulationsTrajectory(0.0, 0.0, 0.0, tmax, dt, freqs[i]);
                populations[i] = pops.sum() / tmax;

                // store trajectory
                std::unique_lock lock(file_mutex);
                std::string groupName = std::to_string(i);
                groupName.insert(groupName.begin(), groupNameLen - groupName.size(), '0');
                auto group = root.CreateSubgroup(groupName);
                group.CreateAttribute("frequency", { 1 });
                group.StoreAttribute("frequency", &freqs[i]);
                auto tStorage = group.CreateDataset("t", { static_cast<std::size_t>(ts.size()) });
                auto pStorage = group.CreateDataset("populations", { static_cast<std::size_t>(pops.size()) });
                tStorage.StoreMatrix(ts);
                pStorage.StoreMatrix(pops);

                progressBar.IncrementCount();
            });
        }

        progressBar.WaitUntilFinished();

        // store overall
        auto fStorage = root.CreateDataset("frequencies", { static_cast<std::size_t>(freqs.size()) });
        auto pStorage = root.CreateDataset("populations", { static_cast<std::size_t>(populations.size()) });
        fStorage.Store(freqs.data());
        pStorage.Store(populations.data());

    }

    if (MoveFile(filename, dstPath))
        std::cout << "Finished. Data saved to " << dstPath << std::endl;
    else
        std::cout << "Finished. Data saved locally to " << filename << std::endl;
    
    return 0;
}
