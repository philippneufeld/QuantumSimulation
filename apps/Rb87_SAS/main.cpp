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

using namespace QSim;
using namespace Eigen;

class CRb87SASApp : public SimulationApp
{
    constexpr static double decay = 6.065e6;
public:
    CRb87SASApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
        double intProbe = NLevelLaser::RabiToIntensity(dip, 3.5e6);
        double intPump = NLevelLaser::RabiToIntensity(dip, 10e6);       
        constexpr double freq = SpeedOfLight_v / 780.241e-9;

        // Create system
        m_system.SetLevel(0, 0.0);
        m_system.SetLevel(1, freq);
        m_system.SetDipoleElement(0, 1, dip);
        m_system.AddLaser(0, 1, intProbe, false);
        m_system.AddLaser(0, 1, intPump, true);
        m_system.SetDecay(1, 0, 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr std::size_t cnt = 501;
        Matrix<double, 2, cnt> detunings;
        detunings.row(0) = RowVectorXd::LinSpaced(cnt, -1e9, 1e9);
        detunings.row(1) = detunings.row(0);

        simdata.CreateDataset("Detunings", { 2, cnt }).StoreMatrix(detunings);
        simdata.CreateDataset("Populations", { cnt });
    }

    virtual void Continue(DataFileGroup& simdata)  override
    {
        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();
        VectorXd populations(detunings.cols());

        ThreadPool pool; 
        ProgressBar progress(detunings.cols());

        // dt << Rabi^-1, Doppler^-1, detuning^-1
        double dt = 1e-10;
        constexpr double tint = 1.0 / decay;

        // start calculation
        auto rho0 = m_system.CreateGroundState();
        for (std::size_t i = 0; i < detunings.cols(); i++)
        {
            pool.Submit([&, i=i](){ 
                populations[i] = m_doppler.Integrate([&](double vel)
                { 
                    auto rho = m_system.GetDensityMatrixAv(detunings.col(i), rho0, vel, 0.0, tint, 0.25*tint, dt);
                return std::real(rho(1, 1));
                });
                progress.IncrementCount();
            });
        } 
        progress.WaitUntilFinished();
        
        simdata.GetDataset("Populations").StoreMatrix(populations);
        SetFinished(simdata);
    }

    virtual void Plot(DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(0).eval();
        auto y_axis = simdata.GetDataset("Populations").LoadMatrix();

        PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    TNLevelSystemSC<2> m_system;
    TDopplerIntegrator<> m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87SASApp app;
    return app.Run(argc, argv);
}
