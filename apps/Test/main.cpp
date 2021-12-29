// Philipp Neufeld, 2021-2022

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Executor/Executor.h>
#include <QSim/Python/Plotting.h>
#include <QSim/Util/CLIProgressBar.h>

#include <chrono>

#include <QSim/Math/CurveFit.h>

using namespace std::chrono_literals;

double func(double x, double a, double g, double x0)
{
    std::this_thread::sleep_for(1ms);
    auto y = x - x0;
    return a / (y*y + g*g/4);
}

int main(int argc, const char* argv[])
{
    QSim::DefaultExecutor executor;
    auto x = QSim::CreateLinspace(-5.0, 10.0, 100);
    auto y = QSim::CreateZerosLike(x);
    executor.Map([](double x){ return func(x, 2, 3, 1.5) + rand() / 2.5e10; }, 
        y, x.Data(), x.Data() + x.Size());

    double a = 1;
    double g = 1;
    double x0 = 1;
    
    QSim::CurveFit(executor, func, x, y, a, g, x0);

    auto xfit = QSim::CreateLinspace(-5.0, 10.0, 500);
    auto yfit = QSim::CreateZerosLike(xfit);
    executor.Map([=](double x) { return func(x, a, g, x0); }, 
        yfit, xfit.Data(), xfit.Data() + xfit.Size());

    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(x.Data(), y.Data(), x.Size(), "Data", ".");
    ax.Plot(xfit.Data(), yfit.Data(), xfit.Size(), "Fit", "-");
    matplotlib.RunGUILoop();

    return 0;
}
