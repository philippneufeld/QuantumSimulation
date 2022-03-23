// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Math/Ode.h>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergAtom.h>
#include <QSim/Rydberg/StarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int fevs = 0;

double func(double x, double y)
{
    fevs++;
    return -x*y*y;
}


int main(int argc, const char* argv[])
{    
    std::vector<double> xs;
    std::vector<double> ys;

    xs.push_back(0.0);
    ys.push_back(1.0);
    
    double dx = 0.01;
    
    TODEStepper<ODEAd54DPPolicy> stepper;

    for(; xs.back() <= 200.0;)
    {
        auto [dy, err] = stepper.StepWithErrorEst(func, ys.back(), xs.back(), dx);
        if (std::abs(err) < 1e-7)
        {
            xs.push_back(xs.back() + dx);
            ys.push_back(ys.back() + dy);

            if (std::abs(err) < 1e-8)
            {
                dx *= 1.5;
            }
        }
        else
        {
            dx /= 1.5;
        }
    }

    std::cout << "Steps: " << fevs << std::endl;

    std::vector<double> ysErr;
    for (int i = 0; i < xs.size(); i++)
        ysErr.push_back(std::abs(ys[i] - (2.0/(xs[i]*xs[i]+2))));

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    // ax.Plot(xs.data(), ys.data(), xs.size(), "", "-C0");
    ax.Plot(xs.data(), ysErr.data(), xs.size(), "", "-C1");

    matplotlib.RunGUILoop();
#endif
    
    return 0;
}
