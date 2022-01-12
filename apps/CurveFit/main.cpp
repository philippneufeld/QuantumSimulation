// Philipp Neufeld, 2021-2022

#include <chrono>

#include <QSim/Math/CurveFit.h>
#include <QSim/Executor/Executor.h>
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
    auto x = QSim::CreateLinspace(-5.0, 10.0, 100);
    auto y = QSim::CreateZerosLike(x);
    executor.Map([](double x){ return func(x, 2, 3, 1.5) + 0.1 * randf(); }, 
        y, x.Data(), x.Data() + x.Size());

    double a = 1;
    double g = 1;
    double x0 = 1;
    
    QSim::CurveFit(executor, func, x, y, a, g, x0);

    auto xfit = QSim::CreateLinspace(-5.0, 10.0, 500);
    auto yfit = QSim::CreateZerosLike(xfit);
    executor.Map([=](double x) { return func(x, a, g, x0); }, 
        yfit, xfit.Data(), xfit.Data() + xfit.Size());

#ifdef QSIM_PYTHON3
    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(x.Data(), y.Data(), x.Size(), "Data", ".");
    ax.Plot(xfit.Data(), yfit.Data(), xfit.Size(), "Fit", "-");
    matplotlib.RunGUILoop();
#else
    std::cout << "Plotting disabled" << std::endl;
#endif

    return 0;
}
