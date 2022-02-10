// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Executor/ThreadPool.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/NLevelSystemQM.h>
#include <QSim/NLevel/NLevelSystemSC.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Math/Functor.h>
#include <QSim/Math/Differentiation.h>

using namespace QSim;

struct Test1
{
    double operator()(double x) { return std::sin(x); }
};

struct Test2
{
    double operator()(Eigen::Vector2d x) { return std::sin(x[0]) + std::cos(x[1]); }
};

struct Generic {};

template<typename Func>
void foo(Func& func) { std::cout << func(1.0) << std::endl; }
template<typename Func>
void bar(TFunctor<Func, double, std::tuple<double>>& func) { std::cout << func(1.0) << std::endl; }

int main(int argc, const char* argv[])
{
    Test2 test;
    auto f1 = CreateFunctor(test);

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(500, -5.0, 5.0);
    Eigen::VectorXd y = x;
    Eigen::VectorXd dy = x;
    
    TDiff1O2<Eigen::Vector2d> diff;

    for (int i = 0; i < x.size(); i++)
    {
        y[i] = f1(Eigen::Vector2d{x[i], 0.0});
        dy[i] = diff.Differentiate(f1, Eigen::Vector2d{x[i], 0.0}, Eigen::Vector2d::Unit(0) * 1e-3);
    }

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size());
    ax.Plot(x.data(), dy.data(), x.size());
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
