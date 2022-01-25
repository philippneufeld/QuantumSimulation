// Philipp Neufeld, 2021-2022

#include <chrono>

#include <QSim/Math/CurveFit.h>
#include <QSim/Executor/ThreadPool.h>
#include <QSim/Util/CLIProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace std::chrono_literals;

double func(double x, double a, double g, double x0)
{
    auto y = x - x0;
    return a / (y*y + g*g/4);
}

float randf()
{
    return 2 * (static_cast<double>(rand()) / RAND_MAX) - 1;
}

int main(int argc, const char* argv[])
{
    QSim::DefaultExecutor executor;
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(100, -5.0, 10.0);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(x.size());
    executor.Map([](double x){ return func(x, 2, 3, 1.5) + 0.1 * randf(); }, 
        y, x.data(), x.data() + x.size());

    double a = 1;
    double g = 1;
    double x0 = 1;
    
    QSim::CurveFit(executor, func, x, y, a, g, x0);

    Eigen::VectorXd xfit = Eigen::VectorXd::LinSpaced(1000, -5.0, 10.0);
    Eigen::VectorXd yfit = Eigen::VectorXd::Zero(xfit.size());
    executor.Map([=](double x) { return func(x, a, g, x0); }, 
        yfit, xfit.data(), xfit.data() + xfit.size());

#ifdef QSIM_PYTHON3
    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size(), "Data", ".");
    ax.Plot(xfit.data(), yfit.data(), xfit.size(), "Fit", "-");
    matplotlib.RunGUILoop();
#else
    std::cout << "Plotting disabled" << std::endl;
#endif

    return 0;
}
