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

class CRb87MWApp : public QSim::SimulationApp
{
public:

    CRb87MWApp()
    {
        // define parameters
        constexpr double lvl5S = 0;
        constexpr double lvl5P = lvl5S + QSim::SpeedOfLight_v / 780e-9;
        constexpr double lvl53D = lvl5P + QSim::SpeedOfLight_v / 480e-9;
        constexpr double lvl54P = lvl53D - 14.3e9;

        constexpr double dip5S5P = 5.178 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dip5P53D = 2.394e-2 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dip54P53D = 3611 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;

        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip5S5P, 1.5e6);
        constexpr double intCoupling = QSim::GetIntensityFromRabiFrequency(dip5P53D, 0.25e6);
        constexpr double intMicrowave = QSim::GetIntensityFromRabiFrequency(dip54P53D, 10e6);

        // Create system (some decay rates from https://atomcalc.jqc.org.uk)
        std::array<std::string, 4> lvlNames = {"5S", "5P", "54P", "53D"};
        std::array<double, 4> levels = {lvl5S, lvl5P, lvl54P, lvl53D};
        m_system = QSim::TNLevelSystemQM<4>(lvlNames, levels);
        m_system.SetDipoleElementByName("5S", "5P", dip5S5P);
        m_system.SetDipoleElementByName("5P", "53D", dip5P53D);
        m_system.SetDipoleElementByName("54P", "53D", dip54P53D);
        m_system.AddLaserByName("Probe", "5S", "5P", intProbe, false);
        m_system.AddLaserByName("Coupling", "5P", "53D", intCoupling, false);
        m_system.AddLaserByName("Microwaves", "54P", "53D", intMicrowave, false);
        m_system.SetDecayByName("5P", "5S", 6.065e6);
        m_system.SetDecayByName("53D", "5P", 6.053e2);
        m_system.SetDecayByName("54P", "53D", 7.964e1);
        m_system.SetDecayByName("54P", "5S", 3.583e1);

        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(QSim::DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TDynamicMatrix<double> detunings(m_system.GetLaserCount(), cnt);
        detunings.SetRow(QSim::CreateLinspaceRow(-3e7, 3e7, cnt), 
            m_system.GetLaserIdxByName("Probe"));
        simdata.CreateDataset("Detunings", { detunings.Rows(), cnt }).Store(detunings.Data());

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
                return rho.GetAbsCoeff("5S", "5P"); 
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
    QSim::TNLevelSystemQM<4> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87MWApp app;
    return app.Run(argc, argv);
}
