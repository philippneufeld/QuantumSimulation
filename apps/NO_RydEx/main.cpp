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

class CNORydEx : public QSim::SimulationApp
{
public:

    CNORydEx()
    {
        m_scanLaser = 2; // Red laser
        m_desiredLevel = 4; // Ion level

        // define parameters
        constexpr double lvlX = 0;
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

        m_system.SetLevel(0, lvlX);
        m_system.SetLevel(1, lvlA);
        m_system.SetLevel(2, lvlH);
        m_system.SetLevel(3, lvlR);
        m_system.SetLevel(4, lvlI);

        m_system.SetDipoleElement(0, 1, dipXA);
        m_system.SetDipoleElement(1, 2, dipAH);
        m_system.SetDipoleElement(2, 3, dipHR);
        m_system.AddLaser(0, 1, uvInt, false);
        m_system.AddLaser(1, 2, greenInt, true);
        m_system.AddLaser(2, 3, redInt, true);
        
        m_system.SetDecay(1, 0, decayAX + decayTransit);
        m_system.SetDecay(2, 1, decayHA);
        m_system.SetDecay(3, 2, decayRH);
        m_system.SetDecay(3, 4, decayIonizaion);
        m_system.SetDecay(2, 0, decayTransit);
        m_system.SetDecay(3, 0, decayTransit);
        m_system.SetDecay(4, 0, decayTransit);

        m_doppler.SetMass(30.0061 * QSim::AtomicMassUnit_v);
        m_doppler.SetATol(1e-10);
        // m_doppler.SetRTol(1e-8);
        // m_doppler.SetIntegrationWidth(m_scanLaser != "UV" ? 0.35 : 3.5); // peak is narrow
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        Eigen::Matrix<double, 3, cnt> detunings;
        detunings.setZero();
        detunings.row(m_scanLaser) = Eigen::RowVectorXd::LinSpaced(cnt, -50.0e6, 50.0e6);

        simdata.CreateDataset("Detunings", { 3, cnt }).StoreMatrix(detunings);
        simdata.CreateDataset("Population", { cnt });
    }

    virtual void Continue(QSim::DataFileGroup& simdata)  override
    {
        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();
        Eigen::VectorXd population(detunings.cols());
        
        QSim::ThreadPool pool; 
        QSim::ProgressBar progress(detunings.cols());

        // start calculation
        for (std::size_t i = 0; i < detunings.cols(); i++)
        {
            pool.Submit([&, i=i](){ 
                population[i] = m_doppler.Integrate([&](double vel)
                { 
                    auto rho = m_system.GetDensityMatrixSS(detunings.col(i), vel);
                    return std::real(rho(m_desiredLevel, m_desiredLevel));
                });
                progress.IncrementCount();
            });
        } 
        progress.WaitUntilFinished();

        simdata.GetDataset("Population").StoreMatrix(population);
        SetFinished(simdata);
    }

    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(m_scanLaser).eval();
        auto y_axis = simdata.GetDataset("Population").LoadMatrix();
        
        std::string laserNames[] = {"UV", "Green", "Red"};
        std::string levelNames[] = {"X", "A", "H", "R", "Ion"};

        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.SetXLabel(laserNames[m_scanLaser] + " detuning [MHz]");
        ax.SetYLabel(levelNames[m_desiredLevel] + " population");
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    unsigned int m_scanLaser;
    unsigned int m_desiredLevel;
    QSim::TNLevelSystemQM<5> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CNORydEx app;
    return app.Run(argc, argv);
}
