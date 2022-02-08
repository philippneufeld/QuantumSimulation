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

using namespace QSim;

struct Functor1
{
    void operator()(double) {}
    void test(double) {}
};

struct Generic {};

int main(int argc, const char* argv[])
{

    constexpr bool b4 = TIsCallable_test_v<Functor1, double>;
    constexpr bool b5 = TIsCallable_test_v<Generic, double>;

    constexpr bool b6 = TIsCallable_test2_v<Functor1, double>;
    constexpr bool b7 = TIsCallable_test2_v<Generic, double>;

    constexpr bool b1 = std::is_invocable_v<Functor1, double>;
    constexpr bool b2 = std::is_invocable_v<Functor1>;
    constexpr bool b3 = std::is_invocable_v<Generic, double>;

    constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
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

    return 0;
}
