// Philipp Neufeld, 2021-2022

#include <iostream>
#include <random>
#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/NLevelSystemQM.h>
#include <QSim/NLevel/NLevelSystemSC.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Execution/SingleThreaded.h>
#include <QSim/Util/ProgressBar.h>

#include <QSim/Math/Quadrature.h>
#include <QSim/Math/LevenbergMarquard.h>

using namespace QSim;
using namespace Eigen;

double test1(double x, const Vector4d& beta)
{
    return ((beta[0]*x + beta[1])*x + beta[2])*x + beta[3];
}

double randf()
{
    return (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * 2;
}

#include <utility>
int main(int argc, const char* argv[])
{
    Vector4d beta = Vector4d{0.5, 3, -2, 1};

    constexpr int n = 100;
    VectorXd x = VectorXd::LinSpaced(n, -5, 5);
    VectorXd y(n);
    for (int i = 0; i < n; i++)
        y[i] = test1(x[i], beta) + 5*randf();
    Vector4d beta0 = Vector4d{0, 0, 0, 0};

    LevenbergMarquard fit;
    Vector4d betafit = fit.CurveFitV(test1, x, y, beta0, (Vector4d::Ones()*1e-3).eval());

    VectorXd xfit = VectorXd::LinSpaced(500, -5, 5);
    VectorXd yfit(xfit.size());
    for (int i = 0; i < xfit.size(); i++)
        yfit[i] = test1(xfit[i], betafit);

    std::cout << betafit << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.eval().data(), y.data(), x.size(), "Data", "x");
    ax.Plot(xfit.data(), yfit.data(), xfit.size(), "Fit");
    matplotlib.RunGUILoop();
#endif

    return 0;
    
    /*constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
    constexpr double intProbe = GetIntensityFromRabiFrequency(dip, 30.5e6);

    TNLevelSystemSC<DynamicDim_v> system2(2);
    system2.SetLevel(0, 0.0);
    system2.SetLevel(1, SpeedOfLight_v / 780.241e-9);
    system2.SetDecay(1, 0, 6.065e6);
    system2.SetDipoleElement(0, 1, dip);
    system2.AddLaser(0, 1, intProbe, false);

    auto rho0 = system2.CreateGroundState();
    auto[ts, rhos] = system2.GetTrajectory(Eigen::VectorXd::Zero(1), rho0, 0.0, 0.0, 4e-7, 1e-9);

    std::vector<double> pops(rhos.size());
    for (std::size_t i = 0; i<rhos.size(); i++)
        pops[i] = std::real(rhos[i](1, 1));

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(ts.data(), pops.data(), pops.size());
    matplotlib.RunGUILoop();
#endif

    return 0;*/
}
