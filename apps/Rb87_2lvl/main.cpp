// Philipp Neufeld, 2021-2022

#include <QSim/Util/SimulationApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem2.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Executor/ThreadPool.h>
#include <QSim/Util/CLIProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

class CRb87TwoLvlApp : public QSim::SimulationApp
{
public:

    CRb87TwoLvlApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double freq = QSim::SpeedOfLight_v / 780.241e-9;

        // Create system
        std::array<double, 2> levels = {0, freq};
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, freq);
        m_system.SetDipoleElement(0, 1, dip);
        m_system.AddLaser(0, 1, intProbe, false);
        m_system.SetDecay(1, 0, 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr std::size_t cnt = 501;
        Eigen::VectorXd detunings = Eigen::VectorXd::LinSpaced(cnt, -1e9, 1e9);
        simdata.CreateDataset("Detunings", { cnt }).Store(detunings.data());

        // Generate result dataset
        simdata.CreateDataset("AbsCoeffs", { cnt });
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
                return std::imag(rho(0, 1));
            });
            progress.IncrementCount();
            return res;
        };

        // Load detuning axis
        auto detunings1 = simdata.GetDataset("Detunings").Load1DRowVector();
        Eigen::RowVectorXd detunings = Eigen::Map<const Eigen::RowVectorXd>(detunings1.Data(), detunings1.Size());

        // start progress bar
        progress.SetTotal(detunings.cols());
        progress.Start();
        
        // start calculation
        Eigen::VectorXd absCoeffs(detunings.cols());
        auto genDetuning = [&detunings](){static int i=0; return detunings.col(i++).eval(); };
        
        pool.MapNonBlockingG(func, absCoeffs, genDetuning, detunings.cols());
        progress.WaitUntilFinished();

        simdata.GetDataset("AbsCoeffs").Store(absCoeffs.data());
        SetFinished(simdata);
    }

    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis1 = simdata.GetDataset("Detunings").Load1DRowVector();
        auto y_axis1 = simdata.GetDataset("AbsCoeffs").Load1DRowVector();

        Eigen::VectorXd x_axis = Eigen::Map<const Eigen::VectorXd>(x_axis1.Data(), x_axis1.Size());
        Eigen::VectorXd y_axis = Eigen::Map<const Eigen::VectorXd>(y_axis1.Data(), y_axis1.Size());

        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    QSim::TNLevelSystemQM2<2> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87TwoLvlApp app;
    return app.Run(argc, argv);
}
