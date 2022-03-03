// Philipp Neufeld, 2021-2022

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
public:

    


private:
    TNLevelSystemQM<5> m_system;
};



class CNORydEx : public SimulationApp
{
public:

    CNORydEx()
    {
        m_scanLaser = 2; // Red laser

        // 
        constexpr double beamRadius = 1e-3;
        constexpr double mass = 30.0061 * AtomicMassUnit_v;
        double meanv = std::sqrt(2*BoltzmannConstant_v*300/mass);

        // define parameters
        constexpr double lvlX = 0;
        constexpr double lvlA = lvlX + SpeedOfLight_v / 226.97e-9;
        constexpr double lvlH = lvlA + SpeedOfLight_v / 540e-9;
        constexpr double lvlR = lvlH + SpeedOfLight_v / 834.92e-9;
        constexpr double lvlI = 9.27 * ElementaryCharge_v / PlanckConstant_v;

        // dipole matrix elements
        constexpr double dipXA = 0.1595 * Debye_v; // https://doi.org/10.1093/mnras/stx1211
        constexpr double dipAH = 1e-1 * Debye_v;
        constexpr double dipHR = 5e-1 * Debye_v;

        // decay rates
        constexpr double decayAX = 13.8e6; // https://doi.org/10.1063/1.454958
        constexpr double decayHA = 10.0e6;
        constexpr double decayRH = 100.0e6;
        constexpr double decayIonizaion = 6e5;
        double decayTransit = meanv / beamRadius;

        // laser intensities
        double uvInt = NLevelLaser::PowerToIntensity(0.05, beamRadius); // 50mW
        double greenInt = NLevelLaser::PowerToIntensity(1.0, beamRadius); // 1W
        double redInt = NLevelLaser::PowerToIntensity(1.0, beamRadius); // 1W

        double rabiUV = NLevelLaser::IntensityToRabi(dipXA, uvInt);
        double rabiGreen = NLevelLaser::IntensityToRabi(dipAH, greenInt);
        double rabiRed = NLevelLaser::IntensityToRabi(dipHR, redInt);

        m_system.SetLevel(0, lvlX);
        m_system.SetLevel(1, lvlA);
        m_system.SetLevel(2, lvlH);
        m_system.SetLevel(3, lvlR);
        m_system.SetLevel(4, lvlI);

        m_system.SetDipoleElement(0, 1, dipXA);
        m_system.SetDipoleElement(1, 2, dipAH);
        m_system.SetDipoleElement(2, 3, dipHR);
        m_system.AddLaser(NLevelLaser({0, 1}, uvInt, 1.0));
        m_system.AddLaser(NLevelLaser({1, 2}, greenInt, -1.0));
        m_system.AddLaser(NLevelLaser({2, 3}, redInt, -1.0));
        
        m_system.SetDecay(1, 0, decayAX + decayTransit);
        m_system.SetDecay(2, 1, decayHA);
        m_system.SetDecay(3, 2, decayRH);
        m_system.SetDecay(3, 4, decayIonizaion);
        m_system.SetDecay(2, 0, decayTransit);
        m_system.SetDecay(3, 0, decayTransit);
        m_system.SetDecay(4, 0, decayTransit);

        m_doppler.SetMass(30.0061 * AtomicMassUnit_v);
    }

    virtual void Init(DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        Matrix<double, 3, cnt> detunings;
        detunings.setZero();
        detunings.row(m_scanLaser) = RowVectorXd::LinSpaced(cnt, -500.0e6, 500.0e6);

        simdata.CreateDataset("Detunings", { 3, cnt }).StoreMatrix(detunings);
        simdata.CreateDataset("Population", { 5, cnt });
    }

    virtual void Continue(DataFileGroup& simdata)  override
    {
        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();
        MatrixXd population(5, detunings.cols());
        
        // get properties of the system and the lasers
        auto transitions = m_system.GetTransitionFreqs();
        auto dirs = m_system.GetLaserDirs();

        ThreadPool pool; 
        ProgressBar progress(detunings.cols());

        // start calculation
        for (std::size_t i = 0; i < detunings.cols(); i++)
        {
            pool.Submit([&, i=i](){

                auto natural = [&](double vel)
                { 
                    VectorXd laserFreqs = transitions + detunings.col(i);
                    laserFreqs = m_doppler.ShiftFrequencies(laserFreqs, dirs, vel);
                    
                    auto rho = m_system.GetDensityMatrixSS(laserFreqs);
                    return real(rho.diagonal().array()).matrix().eval();
                };

                population.col(i) = m_doppler.Integrate(natural);
                progress.IncrementCount();
            });
        } 
        progress.WaitUntilFinished();

        simdata.GetDataset("Population").StoreMatrix(population);
        SetFinished(simdata);
    }

    virtual void Plot(DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        VectorXd x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(m_scanLaser) / 1e6;
        MatrixXd populations = simdata.GetDataset("Population").LoadMatrix();

        std::string laserNames[] = {"UV", "Green", "Red"};
        std::string levelNames[] = {"X", "A", "H", "R", "Ion"};

        PythonMatplotlib matplotlib;

        // for (std::size_t i = 0; i < 5; i++)
        auto i = 4;
        {
            VectorXd pop = populations.row(i);
            auto figure = matplotlib.CreateFigure();
            auto ax = figure.AddSubplot();
            ax.SetXLabel(laserNames[m_scanLaser] + " detuning [MHz]");
            ax.SetYLabel(levelNames[i] + " population");
            ax.Plot(x_axis.data(), pop.data(), x_axis.size());
        }

        matplotlib.RunGUILoop();
#endif
    }

private:
    unsigned int m_scanLaser;
    TNLevelSystemQM<5> m_system;
    TDopplerIntegrator<> m_doppler;
};

int main(int argc, const char* argv[])
{
    CNORydEx app;
    return app.Run(argc, argv);
}
