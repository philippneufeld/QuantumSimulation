// Philipp Neufeld, 2021-2022

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>
#include <QSim/Util/CLIProgressBar.h>

#include <chrono>

#include <QSim/Math/CurveFit.h>

double func(double x, double a, double g, double x0)
{
    auto y = x - x0;
    return a / (y*y + g*g/4);
}

int main(int argc, const char* argv[])
{
    QSim::ThreadPool pool;
    auto x = QSim::CreateLinspace(-5.0, 10.0, 40);
    auto y = QSim::CreateZerosLike(x);

    for (std::size_t i = 0; i < x.Size(); i++)
        y[i] = func(x[i], 2, 3, 1.5) + rand() / 1e10;
    
    QSim::TStaticColVector<double, 3> step;
    step[0] = 1e-3; step[1] = 1e-3; step[2] = 1e-3;
    
    double a = 1;
    double g = 1;
    double x0 = 1;
    
    QSim::CurveFit(pool, func, x, y, a, g, x0);

    auto xfit = QSim::CreateLinspace(-5.0, 10.0, 500);
    auto yfit = QSim::CreateZerosLike(xfit);

    for (std::size_t i = 0; i < xfit.Size(); i++)
        yfit[i] = func(xfit[i], a, g, x0);
    
    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(x.Data(), y.Data(), x.Size(), "Data", ".");
    ax.Plot(xfit.Data(), yfit.Data(), xfit.Size(), "Fit", "-");
    matplotlib.RunGUILoop();

    return 0;
}
