// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/NLevelSystemQM.h>
#include <QSim/NLevel/NLevelSystemSC.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Util/Functor.h>
#include <QSim/Math/Differentiation.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Execution/SingleThreaded.h>
#include <QSim/Util/ProgressBar.h>

#include <QSim/Execution/Progress.h>

using namespace QSim;


double test1(Eigen::Vector2d x) 
{ 
    return std::sin(x[0]) + std::cos(x[1]); 
}

Eigen::Vector2d test2(Eigen::Vector2d x)
{
    return { 2*x[0] + x[1], x[0]*x[1] -3*x[0] + x[1]*x[1] };
}

Eigen::VectorXd test3(Eigen::Vector2d x) 
{ 
    return Eigen::Matrix<double, 1, 1>{std::sin(x[0]) + std::cos(x[1])}; 
}

double test4(Eigen::VectorXd x)
{
    return 2*x[0]*x[0]+3*x[0]+5;
}

double test5(double x)
{
    return 2*x*x+3*x+5;
}



// template<typename Ty>
// using TMatrixEvalType_t = typename TMatrixEvalType<Ty>::type;

// template<typename Ty, typename=std::void_t<decltype(std::declval<std::decay_t<ty2>>().eval())>



int main(int argc, const char* argv[])
{
    auto func = test5;

    ThreadPool pool;
    
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(500, -5.0, 5.0);
    Eigen::VectorXd y = x;
    Eigen::VectorXd dy = x;
    
    Internal::TJacobianHelper diff;
    auto J1 = diff.Jacobian(test2, Eigen::Vector2d{0.0, 0.0}, Eigen::VectorXd::Ones(2) * 1e-3);
    auto J2 = diff.Jacobian(test5, 0.0, /*Eigen::VectorXd::Ones(1) * */ 1e-3);

    std::cout << J1 << std::endl;
    std::cout << J2 << std::endl;

    std::cout << Eigen::Matrix<double, 1, 1>{5.0} / 2.0;

    return 0;

    /*auto progress = pool.CreateProgressTracker(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        pool.Submit([&, i=i]()
        {
            y[i] = func(Eigen::Vector2d{x[i], 0.0});
            dy[i] = diff.Differentiate(func, Eigen::Vector2d{x[i], 0.0}, Eigen::Vector2d::Unit(0) * 1e-3);
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size());
    ax.Plot(x.data(), dy.data(), x.size());
    matplotlib.RunGUILoop();
#endif

    return 0;*/

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
