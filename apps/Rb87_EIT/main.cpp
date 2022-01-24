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

class CRb87EITApp : public QSim::SimulationApp
{
public:

    CRb87EITApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10.0e6);

        // Create system
        m_system.SetLevel(0, -4.271e9);
        m_system.SetLevel(1, 2.563e9);
        m_system.SetLevel(2, QSim::SpeedOfLight_v / 780.241e-9);
        m_system.SetDipoleElement(0, 2, dip);
        m_system.SetDipoleElement(1, 2, dip);
        m_system.AddLaser(0, 2, intProbe, false);
        m_system.AddLaser(1, 2, intPump, false);
        m_system.SetDecay(2, 0, 3.0/8.0 * 6.065e6);
        m_system.SetDecay(2, 1, 5.0/8.0 * 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        Eigen::Matrix<double, 2, cnt> detunings;
        detunings.setZero();
        detunings.row(0) = Eigen::RowVectorXd::LinSpaced(cnt, -100.0e6, 100.0e6);
        
        simdata.CreateDataset("Detunings", { 2, cnt }).StoreMatrix(detunings);
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
                return std::imag(rho(0, 2)); 
            });
            progress.IncrementCount();
            return res;
        };

        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();

        // start progress bar
        progress.SetTotal(detunings.cols());
        progress.Start();
        
        // start calculation
        Eigen::VectorXd absCoeffs(detunings.cols());
        auto genDetuning = [&detunings](){static int i=0; return detunings.col(i++).eval(); };
        
        pool.MapNonBlockingG(func, absCoeffs, genDetuning, detunings.cols());
        progress.WaitUntilFinished();
        
        simdata.GetDataset("AbsCoeffs").StoreMatrix(absCoeffs);
        SetFinished(simdata);
    }


    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(0);
        auto y_axis = simdata.GetDataset("AbsCoeffs").LoadMatrix();
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    QSim::TNLevelSystemQM2<3> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87EITApp app;
    return app.Run(argc, argv);
}
