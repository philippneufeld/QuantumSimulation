// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Util/CLIProgressBar.h>

class CRb87TwoLvlApp : public QSim::CalcApp
{
public:

    CRb87TwoLvlApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double freq = QSim::SpeedOfLight_v / 780.241e-9;

        // Create system
        std::array<std::string, 2> lvlNames = {"S1_2", "P3_2"};
        std::array<double, 2> levels = {0, freq};
        m_system = QSim::TNLevelSystemQM<2>(lvlNames, levels);
        m_system.SetDipoleElementByName("S1_2", "P3_2", dip);
        m_system.AddLaserByName("Probe", "S1_2", "P3_2", intProbe, false);
        m_system.SetDecayByName("P3_2", "S1_2", 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;
        
        QSim::CLIProgBarInt progress;
        auto func = [&](auto dets)
        {
            auto res = m_doppler.Integrate([&](double vel)
            { 
                auto rho = m_system.GetDensityMatrixSS(dets, vel);
                return rho.GetAbsCoeff("S1_2", "P3_2"); 
            });
            progress.IncrementCount();
            return res;
        };

        // Generate detuning axis
        auto detunings = QSim::CreateLinspaceRow(-1e9, 1e9, 501);
        
        // start progress bar
        progress.SetTotal(detunings.Cols());
        progress.Start();
        
        // start calculation
        QSim::TDynamicRowVector<double> absCoeffs = pool.Map(
            func, detunings.GetColIterBegin(), detunings.GetColIterEnd());
        
        this->StoreMatrix("Detunings", detunings);
        this->StoreMatrix("AbsCoeffs S1/2->P3/2", absCoeffs);
    }

    virtual void Plot() override
    {
        auto x_axis = this->LoadMatrix("Detunings");
        auto y_axis = this->LoadMatrix("AbsCoeffs S1/2->P3/2");
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
    }
private:
    QSim::TNLevelSystemQM<2> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87TwoLvlApp app;
    return app.Run(argc, argv);
}
