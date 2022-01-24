// Philipp Neufeld, 2021-2022

#include <QSim/Util/SimulationApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Executor/ThreadPool.h>
#include <QSim/Util/CLIProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

class CNORydEx : public QSim::SimulationApp
{
    constexpr static std::size_t DIM = 5;
public:

    CNORydEx()
    {
        m_scanLaser = "Red";
        m_desiredLevel = "Ion";

        // define parameters
        constexpr double lvlX = 0;
        constexpr double lvlX2 = 10e6;
        constexpr double lvlA = lvlX + QSim::SpeedOfLight_v / 226.97e-9;
        constexpr double lvlH = lvlA + QSim::SpeedOfLight_v / 540e-9;
        constexpr double lvlR = lvlH + QSim::SpeedOfLight_v / 834.92e-9;
        constexpr double lvlI = 9.27 * QSim::ElementaryCharge_v / QSim::PlanckConstant_v;

        // dipole matrix elements
        constexpr double dipXA = 0.1595 * QSim::Debye_v; // https://doi.org/10.1093/mnras/stx1211
        constexpr double dipAH = 2e-3 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dipHR = 1e-3 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;

        // decay rates
        constexpr double decayAX = 13.8e6; // https://doi.org/10.1063/1.454958
        constexpr double decayHA = 1.0e6;
        constexpr double decayRH = 1.0e5;
        constexpr double decayIonizaion = 1e4;
        constexpr double decayTransit = 1e4;

        // laser intensities
        constexpr double uvInt = QSim::GetIntensityFromPower(0.02, 1e-3); // 20mW
        constexpr double greenInt = QSim::GetIntensityFromPower(1.0, 1e-3); // 1W
        constexpr double redInt = QSim::GetIntensityFromPower(1.0, 1e-3); // 1W

        std::array<std::string, DIM> lvlNames = {"X", "A", "H", "R", "Ion"};
        std::array<double, DIM> levels = {lvlX, lvlA, lvlH, lvlR, lvlI};
        m_system = QSim::TNLevelSystemQM<DIM>(lvlNames, levels);
        m_system.SetDipoleElementByName("X", "A", dipXA);
        m_system.SetDipoleElementByName("X", "A", dipXA);
        m_system.SetDipoleElementByName("A", "H", dipAH);
        m_system.SetDipoleElementByName("H", "R", dipHR);
        m_system.AddLaserByName("UV", "X", "A", uvInt, false);
        m_system.AddLaserByName("Green", "A", "H", greenInt, true);
        m_system.AddLaserByName("Red", "H", "R", redInt, true);
        
        m_system.SetDecayByName("A", "X", decayAX + decayTransit);
        m_system.SetDecayByName("H", "A", decayHA);
        m_system.SetDecayByName("R", "H", decayRH);
        m_system.SetDecayByName("R", "Ion", decayIonizaion);
        m_system.SetDecayByName("H", "X", decayTransit);
        m_system.SetDecayByName("R", "X", decayTransit);
        m_system.SetDecayByName("Ion", "X", decayTransit);

        m_doppler.SetMass(30.0061 * QSim::AtomicMassUnit_v);
        m_doppler.SetATol(1e-10);
        // m_doppler.SetRTol(1e-8);
        // m_doppler.SetIntegrationWidth(m_scanLaser != "UV" ? 0.35 : 3.5); // peak is narrow
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TDynamicMatrix<double> detunings(m_system.GetLaserCount(), cnt);
        detunings.SetRow(QSim::CreateLinspaceRow(-5e7, 5e7, cnt), 
            m_system.GetLaserIdxByName(m_scanLaser));
        simdata.CreateDataset("Detunings", { detunings.Rows(), cnt }).Store(detunings.Data());

        // Generate result dataset
        simdata.CreateDataset("Population", { cnt });
    }

    virtual void Continue(QSim::DataFileGroup& simdata)  override
    {
        QSim::ThreadPoolExecutor pool;
        
        QSim::CLIProgBarInt progress;
        auto func = [&](auto dets)
        {
            auto res = m_doppler.Integrate([&](double vel)
            { 
                auto rho = m_system.GetDensityMatrixSS(dets, vel);
                return rho.GetPopulation(m_desiredLevel); 
            });
            progress.IncrementCount();
            return res;
        };

        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").Load2DMatrix();

        progress.SetTotal(detunings.Cols());
        progress.Start();
        
        auto absCoeffs = QSim::CreateZeros(detunings.Cols());
        pool.Map(func, absCoeffs, detunings.GetColIterBegin(), detunings.GetColIterEnd());

        simdata.GetDataset("Population").Store(absCoeffs.Data());
        SetFinished(simdata);
    }

    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto detunings = simdata.GetDataset("Detunings").Load2DMatrix();
        auto x_axis = detunings.GetRow(m_system.GetLaserIdxByName(m_scanLaser));
        auto y_axis = simdata.GetDataset("Population").Load1DRowVector();
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.SetXLabel(m_scanLaser + " detuning [MHz]");
        ax.SetYLabel(m_desiredLevel + " population");
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    std::string m_scanLaser;
    std::string m_desiredLevel;
    QSim::TNLevelSystemQM<DIM> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CNORydEx app;
    return app.Run(argc, argv);
}
