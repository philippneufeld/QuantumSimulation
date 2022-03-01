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

#include <iostream>

using namespace QSim;
using namespace Eigen;

class CRb87EITApp : public SimulationApp
{
public:

    CRb87EITApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
        double intProbe = NLevelLaser::RabiToIntensity(dip, 3.5e6);
        double intPump = NLevelLaser::RabiToIntensity(dip, 10.0e6);

        // Create system
        m_system.SetLevel(0, -4.271e9);
        m_system.SetLevel(1, 2.563e9);
        m_system.SetLevel(2, SpeedOfLight_v / 780.241e-9);
        m_system.SetDipoleElement(0, 2, dip);
        m_system.SetDipoleElement(1, 2, dip);
        m_system.AddLaser(NLevelLaser({0, 2}, intProbe, 1.0));
        m_system.AddLaser(NLevelLaser({1, 2}, intPump, 1.0));
        m_system.SetDecay(2, 0, 3.0/8.0 * 6.065e6);
        m_system.SetDecay(2, 1, 5.0/8.0 * 6.065e6);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void Init(DataFileGroup& simdata) override
    {
        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        Matrix<double, 2, cnt> detunings;
        detunings.setZero();
        detunings.row(0) = RowVectorXd::LinSpaced(cnt, -100.0e6, 100.0e6);
        
        simdata.CreateDataset("Detunings", { 2, cnt }).StoreMatrix(detunings);
        simdata.CreateDataset("AbsCoeffs", { cnt });
    }

    virtual void Continue(DataFileGroup& simdata)  override
    {
        // Load detuning axis
        auto detunings = simdata.GetDataset("Detunings").LoadMatrix();
        VectorXd absCoeffs(detunings.cols());
        
        ThreadPool pool; 
        ProgressBar progress(detunings.cols());

        // start calculation
        for (std::size_t i = 0; i < detunings.cols(); i++)
        {
            pool.Submit([&, i=i](){ 
                auto natural = [&](double vel)
                { 
                    auto rho = m_system.GetDensityMatrixSS(detunings.col(i), vel);
                    return std::imag(rho(0, 2));
                };

                absCoeffs[i] = m_doppler.Integrate(natural);
                progress.IncrementCount();
            });
        }
        progress.WaitUntilFinished();
        
        simdata.GetDataset("AbsCoeffs").StoreMatrix(absCoeffs);
        SetFinished(simdata);
    }


    virtual void Plot(DataFileGroup& simdata) override
    {
#ifdef QSIM_PYTHON3
        auto x_axis = simdata.GetDataset("Detunings").LoadMatrix().row(0).eval();
        auto y_axis = simdata.GetDataset("AbsCoeffs").LoadMatrix();
        
        PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
#endif
    }

private:
    TNLevelSystemQM<3> m_system;
    TDopplerIntegrator<> m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87EITApp app;
    return app.Run(argc, argv);
}
