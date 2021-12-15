// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>

class CRb87EITApp : public QSim::CalcApp
{
public:
    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;
        
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10.0e6);

        // Create system
        std::array<std::string, 3> lvlNames = {"S1_2_F1", "S1_2_F2", "P3_2"};
        std::array<double, 3> levels = {-4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9};
        QSim::TNLevelSystemQM<3> system(lvlNames, levels);
        system.SetDipoleElementByName("S1_2_F1", "P3_2", dip);
        system.SetDipoleElementByName("S1_2_F2", "P3_2", dip);
        system.AddLaserByName("Probe", "S1_2_F1", "P3_2", intProbe, false);
        system.AddLaserByName("Pump", "S1_2_F2", "P3_2", intPump, false);
        system.SetDecayByName("P3_2", "S1_2_F1", 3.0/8.0 * 6.065e6);
        system.SetDecayByName("P3_2", "S1_2_F2", 5.0/8.0 * 6.065e6);
        system.SetMass(1.44316060e-25);

        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TStaticMatrix<double, 2, cnt> detunings;
        detunings.SetRow(QSim::CreateLinspaceRow(-100.0e6, 100.0e6, cnt), 
            system.GetLaserIdxByName("Probe"));

        auto func = [&](auto dets)
        {
            auto rho = system.GetDensityMatrixSS(dets);
            return rho.GetAbsCoeff("S1_2_F1", "P3_2");
        }; 
        QSim::TDynamicRowVector<double> absCoeffs = pool.Map(
            func, detunings.GetColIterBegin(), detunings.GetColIterEnd());
        
        this->StoreMatrix("Detunings", detunings);
        this->StoreMatrix("AbsCoeffs S1/2 F1->P3/2", absCoeffs);
    }

    virtual void Plot() override
    {
        auto x_axis = this->LoadMatrix("Detunings").GetRow(0);
        auto y_axis = this->LoadMatrix("AbsCoeffs S1/2 F1->P3/2");
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
    }
};

int main(int argc, const char* argv[])
{
    CRb87EITApp app;
    return app.Run(argc, argv);
}
