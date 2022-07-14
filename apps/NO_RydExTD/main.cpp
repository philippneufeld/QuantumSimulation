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

// #define DOPPLER_ADAPTIVE

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
        m_pressure(100.0) // 1 mbar = 100Pa
    {
        m_mass = 30.0061 * AtomicMassUnit_v;
        
        // setup system levels
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, m_system.GetLevel(0) + SpeedOfLight_v / 226.97e-9);
        m_system.SetLevel(2, m_system.GetLevel(1) + SpeedOfLight_v / 540.0e-9);
        m_system.SetLevel(3, m_system.GetLevel(2) + SpeedOfLight_v / 834.92e-9);
        m_system.SetLevel(4, 9.27 * ElementaryCharge_v / PlanckConstant_v);

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

        std::cout << "Transit decay rate: " << decayTransit / 1e6 << "MHz" << std::endl;

        // add collisional ionization term
        // collisions between rydberg and ground state
        double rateRX = GetCollisionRate(rydRadius + xRadius);
        m_system.AddDecay(3, 4, rateRX);

        // add recombination rate (assume to be the same as ionization term for now)
        m_system.AddDecay(4, 0, rateRX);

        std::cout << "Rydberg collision rate: " << rateRX / 1e6 << "MHz" << std::endl;
        std::cout << "GS collision rate: " << GetCollisionRate(2*xRadius) / 1e6 << "MHz" << std::endl;

        // initialize helper variables
        m_laserDirs = Vector3d{1.0, -1.0, -1.0};
    }

    static void SetCouplings(TNLevelSystem<5, true>& sys, double greenModFreq)
    {
        sys.ClearCouplings();

        sys.AddCoupling(0, 1, [](double) { return 3.5e6; });

        double greenRabi = 1.0e6;
        if (greenModFreq == 0)
            sys.AddCoupling(1, 2, [=](double t){ return 1.0; });
        else
        {
            double modPeriod = 1.0 / greenModFreq;
            sys.AddCoupling(1, 2, [=](double t){ t = std::fmod(t, modPeriod); return t <= modPeriod / 2 ? 1.0 : 0.0; });
        }

        sys.AddCoupling(2, 3, [](double) { return 2.0e6; });
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

    auto GetPopulationsTrajectory(
        double uvDet, double greenDet, double redDet, 
        double t, double dt, double modFreq) const
    {
        // make a local copy of the system (makes concurrent calls to this function possible)
        auto sys = m_system;

        // modulate green laser
        SetCouplings(sys, std::abs(modFreq));
        Vector3d resonanceFreqs = sys.GetCouplingResonanceFreqs();

        // calculate ionization rate (rate0 + rate1*rho_ion)
        double rate0 = GetCollisionRate(rydRadius + xRadius); // constant term
        double rate1 = GetCollisionRate(rydRadius);

        // starting point of the calculation
        auto rho0 = sys.CreateGroundState();
        using Rho_t = decltype(rho0);

        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil(t / dt));
        dt = t / steps;

        // limit runtime by sacrificing some precision in the doppler integration
#ifndef DOPPLER_ADAPTIVE
        TDopplerIntegrator<QuadSimpsonPolicy> doppler(m_mass, m_temperature, 1001);
#else
        TDopplerIntegrator<> doppler(m_mass, m_temperature);
        doppler.SetIntegrationRTol(1e-3);
#endif
        doppler.SetIntegrationWidth(1.0);

        auto dopplerCalc = [&](double vel)
        {
            VectorXd ionPops(steps);

            // calculate laser frequencies
            Vector3d detunings(uvDet, greenDet, redDet);
            Vector3d laserFreqs = resonanceFreqs + detunings;
            
            // adjust laser frequencies for doppler shift
            laserFreqs = doppler.ShiftFrequencies(laserFreqs, m_laserDirs, vel);

            // define function to be integrated (non-linearity included)
            auto func = [&](double x, const Rho_t& rho) 
            {
                // non-linear dependece on ion population
                double rhoIon = std::real(rho(4, 4));
                sys.SetDecay(3, 4, rate0 + rate1*rhoIon);

                return sys.GetDensityOpDerivative(laserFreqs, rho, x); 
            };

            // integrate trajectory
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

        // auto trajectory = dopplerCalc(0.0);
        auto trajectory = doppler.Integrate(dopplerCalc);

        VectorXd timeAxis = sys.GetTrajectoryTimeaxis(0.0, dt, trajectory.size());
        return std::make_pair(timeAxis, trajectory);
    }

private:
    TNLevelSystem<5, true> m_system;

    // environmental parameters
    double m_temperature;
    double m_mass;
    double m_pressure;

    // laser variables
    double m_beamRadius;
    Vector3d m_laserDirs;
};


int main(int argc, const char* argv[])
{
    NOGasSensorTD gasSensor;

    double fmin = 7500;
    double fmax = 1e9;
    VectorXd freqs = ArrayXd::LinSpaced(300, std::log(fmin), std::log(fmax)).exp();
    VectorXd populations = VectorXd::Zero(freqs.size());

    double dt = 1e-9;
    double tmin = 1e-4;

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
	    	double tsim = std::max(50.0 / freqs[i], tmin);
                auto [ts, pops] = gasSensor.GetPopulationsTrajectory(0.0, 0.0, 0.0, tsim, dt, freqs[i]);

                populations[i] = pops.sum() / tsim;

                // store trajectory
                std::unique_lock lock(file_mutex);
                std::string groupName = std::to_string(i);
                groupName.insert(groupName.begin(), groupNameLen - groupName.size(), '0');
                auto group = root.CreateSubgroup(groupName);

                group.SetAttribute("frequency", freqs[i]);
                group.CreateDataset("t", ts);
                group.CreateDataset("populations", pops);

                progressBar.IncrementCount();
            });
        }

        progressBar.WaitUntilFinished();

        // store overall
        root.CreateDataset("frequencies", freqs);
        root.CreateDataset("populations", populations);
    }

    if (MoveFile(filename, dstPath))
        std::cout << "Finished. Data saved to " << dstPath << std::endl;
    else
        std::cout << "Finished. Data saved locally to " << filename << std::endl;
    
    return 0;
}
