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

#include <QSim/Util/Functor.h>

#include <QSim/Math/Quadrature.h>
#include <QSim/Math/LevenbergMarquardt.h>
#include <QSim/Math/Gamma.h>

using namespace QSim;
using namespace Eigen;

double test1(double x, const Vector4d& beta)
{
    return (((beta[0]/1e8)*x + beta[1])*x + beta[2])*x + beta[3];
}

double randf()
{
    return (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * 2;
}

#include <utility>
int main(int argc, const char* argv[])
{

    /*std::cout << std::tgamma(5) << std::endl;
    std::cout << GammaFunction::GammaP(2, 4) << std::endl;
    std::cout << GammaFunction::InvGammaP(2, GammaFunction::GammaP(2, 4)) << std::endl;

    VectorXd x = VectorXd::LinSpaced(10000, 0.0, 500.0);
    VectorXd y1 = x.unaryExpr([](double x){ return GammaFunction::GammaP(10, x); });
    VectorXd y2 = x.unaryExpr([](double x){ return GammaFunction::GammaP(3, x); });
    VectorXd y3 = x.unaryExpr([](double x){ return GammaFunction::GammaP(1, x); });
    VectorXd y4 = x.unaryExpr([](double x){ return GammaFunction::GammaP(0.5, x); });
    VectorXd y5 = x.unaryExpr([](double x){ return GammaFunction::GammaP(100, x); });
    VectorXd y6 = x.unaryExpr([](double x){ return GammaFunction::GammaP(400, x); });
    
    VectorXd y = VectorXd::LinSpaced(1000, 0.0, 0.999);
    VectorXd x1 = y.unaryExpr([](double y){ return GammaFunction::InvGammaP(200, y); });
    VectorXd x2 = y.unaryExpr([](double y){ return GammaFunction::InvGammaP(0.75, y); });
    

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.eval().data(), y1.data(), x.size(), "a=10");
    ax.Plot(x.eval().data(), y2.data(), x.size(), "a=3");
    ax.Plot(x.eval().data(), y3.data(), x.size(), "a=1");
    ax.Plot(x.eval().data(), y4.data(), x.size(), "a=0.5");
    ax.Plot(x.eval().data(), y5.data(), x.size(), "a=100");
    ax.Plot(x.eval().data(), y6.data(), x.size(), "a=400");
    ax.Plot(x1.eval().data(), y.data(), y.size(), "a=200");
    ax.Plot(x2.eval().data(), y.data(), y.size(), "a=0.75");
    fig.Legend();
    matplotlib.RunGUILoop();
#endif

    return 0;*/

    Vector4d beta = Vector4d{0.5*1e8, 3, -2, 1};

    constexpr int n = 100;
    double yerrMag = 5;
    VectorXd x = VectorXd::LinSpaced(n, -5, 5);
    VectorXd y = x.unaryExpr([&](double x){ return test1(x, beta) + yerrMag*randf(); });
    VectorXd yerrs = VectorXd::Ones(y.size())*yerrMag/1.5;
    Vector4d beta0 = Vector4d{0, 0, 0, 0};

    LevenbergMarquardt fit;
    auto [params, cov, good] = fit.CurveFitV(test1, x, y, yerrs, beta0);
    Vector4d errs = cov.diagonal().cwiseSqrt();

    VectorXd xfit = VectorXd::LinSpaced(500, -5, 5);
    VectorXd yfit = xfit.unaryExpr([p=params](double x){ return test1(x, p); });
    VectorXd yfit2 = xfit.unaryExpr([p=params+errs](double x){ return test1(x, p); });
    VectorXd yfit3 = xfit.unaryExpr([p=params-errs](double x){ return test1(x, p); });
    
    std::cout << params << std::endl << std::endl;
    std::cout << errs.cwiseQuotient(params) * 100 << std::endl << std::endl;
    std::cout << good << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.eval().data(), y.data(), x.size(), "Data", "x");
    ax.Plot(xfit.data(), yfit.data(), xfit.size(), "Fit");
    ax.Plot(xfit.data(), yfit2.data(), xfit.size(), "Fit upper");
    ax.Plot(xfit.data(), yfit3.data(), xfit.size(), "Fit lower");
    fig.Legend();
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
