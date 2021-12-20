// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Util/ThreadPool.h>

class CNORydEx : public QSim::CalcApp
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
        constexpr double lvlI = 0; // 9.27 * QSim::ElementaryCharge_v / QSim::PlanckConstant_v;

        // dipole matrix elements
        constexpr double dipXA = 2e-3 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dipAH = 1e-3 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dipHR = 5e-4 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;

        // decay rates
        constexpr double decayAX = 10.0e6;
        constexpr double decayHA = 1.0e6;
        constexpr double decayRH = 1.0e5;
        constexpr double decayIonizaion = 1e4;
        constexpr double decayTransit = 1e4;

        // laser intensities
        constexpr double uvInt = QSim::GetIntensityFromPower(0.03, 1e-3); // 30mW
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
        // m_doppler.SetATol(1e-10);
        // m_doppler.SetRTol(1e-8);
        m_doppler.SetIntegrationWidth(m_scanLaser != "UV" ? 0.35 : 3.5); // peak is narrow
    }

    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;
        
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TDynamicMatrix<double> detunings(m_system.GetLaserCount(), cnt);
        detunings.SetRow(QSim::CreateLinspaceRow(-3e7, 3e7, cnt), 
            m_system.GetLaserIdxByName(m_scanLaser));

        auto func = [&](auto dets)
        {
            return m_doppler.Integrate([&](double vel)
            { 
                auto rho = m_system.GetDensityMatrixSS(dets, vel);
                return rho.GetPopulation(m_desiredLevel); 
            });
        };

        QSim::TDynamicRowVector<double> absCoeffs = pool.Map(
            func, detunings.GetColIterBegin(), detunings.GetColIterEnd());
        
        this->StoreMatrix("Detunings", detunings);
        this->StoreMatrix("Population " + m_desiredLevel, absCoeffs);
    }

    virtual void Plot() override
    {
        auto detunings = this->LoadMatrix("Detunings");
        auto x_axis = detunings.GetRow(m_system.GetLaserIdxByName(m_scanLaser));
        auto y_axis = this->LoadMatrix("Population " + m_desiredLevel);
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.SetXLabel(m_scanLaser + " detuning [MHz]");
        ax.SetYLabel(m_desiredLevel + " population");
        ax.Plot((x_axis / 1e6).Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
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
