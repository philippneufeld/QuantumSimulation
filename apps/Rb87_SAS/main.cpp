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

class CRb87SASApp : public QSim::SimulationApp
{
    constexpr static double decay = 6.065e6;
public:
    CRb87SASApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10e6);       
        constexpr double freq = QSim::SpeedOfLight_v / 780.241e-9;

        // Create system
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, freq);
        m_system.SetDipoleElement(0, 1, dip);
        m_system.AddLaser(0, 1, intProbe, false);
        m_system.AddLaser(0, 1, intPump, true);
        m_system.SetDecay(1, 0, 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr std::size_t cnt = 501;
        Eigen::Matrix<double, 2, cnt> detunings;
        detunings.row(0) = Eigen::RowVectorXd::LinSpaced(cnt, -1e9, 1e9);
        detunings.row(1) = detunings.row(0);

        simdata.CreateDataset("Detunings", { 2, cnt }).StoreMatrix(detunings);
        simdata.CreateDataset("Populations", { cnt });
    }

    virtual void Continue(QSim::DataFileGroup& simdata)  override
    {
        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();
 
        QSim::ThreadPoolExecutor pool; 
        QSim::CLIProgBar progress(detunings.cols());

        // dt << Rabi^-1, Doppler^-1, detuning^-1
        double dt = 1e-10;
        constexpr double tint = 1.0 / decay;

        auto rho0 = m_system.CreateGroundState();
        auto func = [&](auto dets)
        {
            auto res = m_doppler.Integrate([&](double vel)
            { 
                auto rho = m_system.GetDensityMatrixAv(dets, rho0, vel, 0.0, tint, 0.25*tint, dt);
                return std::real(rho(1, 1));
            });
            progress.IncrementCount();
            return res;
        };

        // start calculation
        Eigen::VectorXd populations(detunings.cols());
        auto genDetuning = [&detunings](){static int i=0; return detunings.col(i++).eval(); };
        
        progress.Start();
        pool.MapNonBlockingG(func, populations, genDetuning, detunings.cols());
        progress.WaitUntilFinished();
        
        simdata.GetDataset("Populations").StoreMatrix(populations);
        SetFinished(simdata);
    }

    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(0).eval();
        auto y_axis = simdata.GetDataset("Populations").LoadMatrix();

        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    QSim::TNLevelSystemSC<2> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87SASApp app;
    return app.Run(argc, argv);
}
