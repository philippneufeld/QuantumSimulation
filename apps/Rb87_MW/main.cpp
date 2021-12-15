// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>

class CRb87MWApp : public QSim::CalcApp
{
public:
    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;
        
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
        QSim::TNLevelSystemQM<4> system(lvlNames, levels);
        system.SetMass(1.44316060e-25);
        system.SetDipoleElementByName("5S", "5P", dip5S5P);
        system.SetDipoleElementByName("5P", "53D", dip5P53D);
        system.SetDipoleElementByName("54P", "53D", dip54P53D);
        system.AddLaserByName("Probe", "5S", "5P", intProbe, false);
        system.AddLaserByName("Coupling", "5P", "53D", intCoupling, false);
        system.AddLaserByName("Microwaves", "54P", "53D", intMicrowave, false);
        system.SetDecayByName("5P", "5S", 6.065e6);
        system.SetDecayByName("53D", "5P", 6.053e2);
        system.SetDecayByName("54P", "53D", 7.964e1);
        system.SetDecayByName("54P", "5S", 3.583e1);
        
        system.SetDopplerIntegrationSteps(15000);

        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TDynamicMatrix<double> detunings(system.GetLaserCount(), cnt);
        detunings.SetRow(QSim::CreateLinspaceRow(-3e7, 3e7, cnt), 
            system.GetLaserIdxByName("Probe"));

        auto func = [&](auto dets)
        {
            auto rho = system.GetDensityMatrixSS(dets);
            return rho.GetAbsCoeff("5S", "5P");
        }; 
        QSim::TDynamicRowVector<double> absCoeffs = pool.Map(
            func, detunings.GetColIterBegin(), detunings.GetColIterEnd());
        
        this->StoreMatrix("Detunings", detunings);
        this->StoreMatrix("AbsCoeffs 5S->5P", absCoeffs);
    }

    virtual void Plot() override
    {
        auto x_axis = this->LoadMatrix("Detunings").GetRow(0);
        auto y_axis = this->LoadMatrix("AbsCoeffs 5S->5P");
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
    }
};

int main(int argc, const char* argv[])
{
    CRb87MWApp app;
    return app.Run(argc, argv);
}
