// Philipp Neufeld, 2021-2022

#include <iostream>
#include <mutex>

#include <Eigen/Dense>

#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>
#include <QSim/Util/PathUtil.h>
#include <QSim/Util/DataFile.h>

using namespace QSim;
using namespace Eigen;

template<bool AM=false>
class TNOGasSensor
{
public:
    using System_t = TNLevelSystem<5, AM>;
    using Rabi_t = typename System_t::Rabi_t;
    using MatOp_t = typename System_t::MatOp_t;

    static constexpr double s_mass = 30.01 * AtomicMassUnit_v;
    constexpr static double s_xRadius = 115e-12; // bond length NO ground state
    constexpr static double s_rydRadius = BohrRadius_v * 25*25; // estimation of rydberg radius at n=25

    TNOGasSensor()
        : m_laserDirs{1.0, -1.0, -1.0},
        m_beamRadius(1e-3), 
        m_temperature(300), 
        m_pressure(100.0), // 1 mbar = 100Pa
        m_dopplerSteps(301)
    {
        // setup system levels
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, m_system.GetLevel(0) + SpeedOfLight_v / 226.97e-9);
        m_system.SetLevel(2, m_system.GetLevel(1) + SpeedOfLight_v / 540.0e-9);
        m_system.SetLevel(3, m_system.GetLevel(2) + SpeedOfLight_v / 834.92e-9);
        m_system.SetLevel(4, 9.27 * ElementaryCharge_v / PlanckConstant_v);

        // add couplings
        m_system.AddCoupling(0, 1, CreateConstRabi(0));
        m_system.AddCoupling(1, 2, CreateConstRabi(0));
        m_system.AddCoupling(2, 3, CreateConstRabi(0));

        // add natural decays
        m_system.AddDecay(1, 0, 13.8e6); // https://doi.org/10.1063/1.454958
        m_system.AddDecay(2, 1, 1.0e6);
        m_system.AddDecay(3, 2, 0.5e6);

        // add transit broadening effect
        double meanVel = std::sqrt(8*BoltzmannConstant_v*m_temperature / (Pi_v * s_mass));
        double decayTransit = meanVel / m_beamRadius;
        m_system.AddDecay(1, 0, decayTransit);
        m_system.AddDecay(2, 0, decayTransit);
        m_system.AddDecay(3, 0, decayTransit);
        m_system.AddDecay(4, 0, decayTransit);

        // calculate collision rates
        m_collRateRX = GetCollisionRate(s_rydRadius + s_xRadius, s_mass); // rydberg-ground state
        m_collRateRE = GetCollisionRate(s_rydRadius, ElectronMass_v); // rydberg-electron
        m_collRateIE = GetCollisionRate(s_xRadius, ElectronMass_v); // ion-electron
    }

    void SetRabiUV(Rabi_t rabi) { m_system.GetCoupling(0) = rabi; }
    void SetRabiGreen(Rabi_t rabi) { m_system.GetCoupling(1) = rabi; }
    void SetRabiRed(Rabi_t rabi) { m_system.GetCoupling(2) = rabi; }

    auto GetPopulationsTrajectory(
        double uvDet, double greenDet, double redDet, double t, double dt) const
    {
        // make a local copy of the system (makes concurrent calls to this function possible)
        System_t sys = m_system;

        Vector3d resonanceFreqs = sys.GetCouplingResonanceFreqs();

        // get rates
        double prevRI = sys.GetDecay(3, 4);
        double prevIX = sys.GetDecay(4, 0);

        // starting point of the calculation
        MatOp_t rho0 = sys.CreateGroundState();

        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil(t / dt));
        dt = t / steps;

        TDopplerIntegrator<QuadMidpointPolicy> doppler(s_mass, m_temperature, m_dopplerSteps);
        doppler.SetIntegrationWidth(1.25);

        auto dopplerCalc = [&](double vel)
        {
            VectorXd ionPops(steps + 1);

            // calculate laser frequencies
            Vector3d detunings(uvDet, greenDet, redDet);
            Vector3d laserFreqs = resonanceFreqs + detunings;
            
            // adjust laser frequencies for doppler shift
            laserFreqs = doppler.ShiftFrequencies(laserFreqs, m_laserDirs, vel);

            // define function to be integrated (non-linearity included)
            auto func = [&](double x, const MatOp_t& rho) 
            {
                // non-linear dependece on ion population
                double rhoIon = std::real(rho(4, 4));
                sys.SetDecay(3, 4, prevRI + GetIonizationRate(rhoIon));
                sys.SetDecay(4, 0, prevIX + GetRecombinationRate(rhoIon));

                return sys.GetDensityOpDerivative(laserFreqs, rho, x); 
            };

            MatOp_t rho = rho0;
            ionPops(0) = std::real(rho(4, 4));

            // integrate trajectory
            TODEIntegrator<ODEAd54CKPolicy> integrator;
            integrator.SetPrecision(1e-13);
            double dtEff = dt; // dt that is controlled by the adaptive stepsize control
            for (int i = 0; i < steps; i++)
            {
                rho = integrator.IntegrateTo(func, rho, i*dt, (i+1)*dt, dtEff);
                ionPops(i + 1) = std::real(rho(4, 4));
            }
            
            // return populations of the levels
            return ionPops;
        };

        auto trajectory = doppler.Integrate(dopplerCalc);
        VectorXd timeAxis = sys.GetTrajectoryTimeaxis(0.0, dt, trajectory.size());
        
        return std::make_pair(timeAxis, trajectory);
    }

    double GetIonizationRate(double ionPop) const
    {
        return m_collRateRX + m_collRateRE * ionPop;
    }

    double GetRecombinationRate(double ionPop) const
    {
        return m_collRateIE * ionPop;
    }

protected:
    double GetCollisionRate(double rEff, double mColl) const
    {
        // pV=nkT
        double kT = BoltzmannConstant_v * m_temperature;
        double n = m_pressure / kT; // number density
        double sigma = Pi_v * std::pow(rEff, 2);
        double mtp = 1.0 / (std::sqrt(2) * sigma * n); // mean free path
        double vExp = std::sqrt(8*kT / (Pi_v*mColl)); // expected velocity
        return vExp / mtp;
    }

    template<bool Dummy=AM>
    static std::enable_if_t<!Dummy, Rabi_t> CreateConstRabi(double rabi) { return rabi; }
    template<bool Dummy=AM>
    static std::enable_if_t<Dummy, Rabi_t> CreateConstRabi(double rabi) { return [=](double) { return rabi; }; }

private:
    System_t m_system;
    const std::size_t m_dopplerSteps;

    double m_collRateRX;
    double m_collRateRE;
    double m_collRateIE;

    // environmental parameters
    double m_temperature;
    double m_pressure;

    // laser variables
    double m_beamRadius;
    const Vector3d m_laserDirs;
};


std::pair<VectorXd, VectorXd> GenerateTrace(TNOGasSensor<true> gasSensor, double rabiAH, double freq, double tmin, double dt, int osc)
{
    double tsim = std::max(osc / freq, std::ceil(tmin*freq)/freq);
    double period = 1.0 / freq;
    dt = std::min(0.025*period, dt);
    gasSensor.SetRabiGreen([=](double t){ t = std::fmod(t, period); return t <= period / 2 ? rabiAH : 0.0; });
    return gasSensor.GetPopulationsTrajectory(0, 0, 0, tsim, dt);
}

std::string RunSimulation(double rabiXA, double rabiAH, double rabiHR, double fmin, 
    double fmax, std::size_t n, double tmin, double dt, int osc)
{
    // initialize gas sensor
    TNOGasSensor<true> gasSensor;
    gasSensor.SetRabiUV([=](double){ return rabiXA; });
    gasSensor.SetRabiGreen([=](double t){ return rabiAH; });
    gasSensor.SetRabiRed([=](double){ return rabiHR; });

    // generate filename
    std::string filename = GenerateFilename("NORydExTD") + ".h5";

    // create file locally
    DataFile file;
    file.Open(filename, DataFile_MUST_NOT_EXIST);
    auto root = file.OpenRootGroup();
    auto groupNameLen = std::to_string(n).size();
    std::mutex file_mutex;

    ThreadPool threadPool;
    ProgressBar progressBar(n);

    // create frequency axis
    VectorXd freqs = ArrayXd::LinSpaced(n, std::log(fmin), std::log(fmax)).exp();

    for (int i = 0; i < freqs.size(); i++)
    {
        threadPool.Submit([&, i=i]()
        {
            auto [ts, pops] = GenerateTrace(gasSensor, rabiAH, freqs[i], tmin, dt, osc);

            // create group
            std::unique_lock lock(file_mutex);
            std::string groupName = std::to_string(i);
            groupName.insert(groupName.begin(), groupNameLen - groupName.size(), '0');
            auto group = root.CreateSubgroup(groupName);

            // Store data
            group.SetAttribute("frequency", freqs[i]);
            group.CreateDataset("t", ts);
            group.CreateDataset("populations", pops);

            progressBar.IncrementCount();
        });
    }

    progressBar.WaitUntilFinished();
    
    return filename;
}

int main(int argc, const char* argv[])
{
    double rabiXA = 5.0e6;
    double rabiAH = 1.0e6;
    double rabiHR = 1.0e6;

    double fmin = 5e4;
    double fmax = 1e9;
    
    double dt = 1e-9;
    double tmin = 1e-5;

    std::string filename = RunSimulation(rabiXA, rabiAH, rabiHR, fmin, fmax, 300, tmin, dt, 2);

    std::string dstPath = GetDefaultAppDataDir("NO_RydExTD") + '/' + filename;
    if (MoveFile(filename, dstPath))
        std::cout << "Finished. Data saved to " << dstPath << std::endl;
    else
        std::cout << "Finished. Data saved locally to " << filename << std::endl;
    
    return 0;

}
