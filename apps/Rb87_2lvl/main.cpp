// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>

class CRb87TwoLvlApp : public QSim::CalcApp
{
public:
    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;
        
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);

        // Create system
        std::array<std::string, 2> lvlNames = {"S1_2", "P3_2"};
        std::array<double, 2> levels = {0, QSim::SpeedOfLight_v / 780.241e-9};
        QSim::TNLevelSystemQM<2> system(lvlNames, levels);
        system.SetDipoleElementByName("S1_2", "P3_2", dip);
        system.AddLaserByName("Probe", "S1_2", "P3_2", intProbe, false);
        system.SetDecayByName("P3_2", "S1_2", 6.065e6);
        system.SetMass(1.44316060e-25);

        auto detunings = QSim::CreateLinspaceRow(-1e9, 1e9, 501);

        auto func = [&](auto dets)
        { 
            auto rho = system.GetDensityMatrixSS(dets);
            return rho.GetAbsCoeff("S1_2", "P3_2");
        }; 
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
};

int main(int argc, const char* argv[])
{
    CRb87TwoLvlApp app;
    return app.Run(argc, argv);
}
