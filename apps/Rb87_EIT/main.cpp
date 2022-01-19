// Philipp Neufeld, 2021-2022

#include <QSim/Util/SimulationApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Executor/Executor.h>
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
        std::array<std::string, 3> lvlNames = {"S1_2_F1", "S1_2_F2", "P3_2"};
        std::array<double, 3> levels = {-4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9};
        m_system = QSim::TNLevelSystemQM<3>(lvlNames, levels);
        m_system.SetDipoleElementByName("S1_2_F1", "P3_2", dip);
        m_system.SetDipoleElementByName("S1_2_F2", "P3_2", dip);
        m_system.AddLaserByName("Probe", "S1_2_F1", "P3_2", intProbe, false);
        m_system.AddLaserByName("Pump", "S1_2_F2", "P3_2", intPump, false);
        m_system.SetDecayByName("P3_2", "S1_2_F1", 3.0/8.0 * 6.065e6);
        m_system.SetDecayByName("P3_2", "S1_2_F2", 5.0/8.0 * 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TStaticMatrix<double, 2, cnt> detunings;
        detunings.SetRow(QSim::CreateLinspaceRow(-100.0e6, 100.0e6, cnt), 
            m_system.GetLaserIdxByName("Probe"));
        simdata.CreateDataset("Detunings", { 2, cnt }).Store(detunings.Data());

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
                return rho.GetAbsCoeff("S1_2_F1", "P3_2"); 
            });
            progress.IncrementCount();
            return res;
        };

        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").Load2DMatrix();

        // start progress bar
        progress.SetTotal(detunings.Cols());
        progress.Start();
        
        // start calculation
        auto absCoeffs = QSim::CreateZeros(detunings.Cols());
        pool.Map(func, absCoeffs, detunings.GetColIterBegin(), detunings.GetColIterEnd());
        
        simdata.GetDataset("AbsCoeffs").Store(absCoeffs.Data());
        SetFinished(simdata);
    }


    virtual void Plot(QSim::DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto detunings = simdata.GetDataset("Detunings").Load2DMatrix();
        auto x_axis = detunings.GetRow(m_system.GetLaserIdxByName("Probe"));
        auto y_axis = simdata.GetDataset("AbsCoeffs").Load1DRowVector();
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    QSim::TNLevelSystemQM<3> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87EITApp app;
    return app.Run(argc, argv);
}
