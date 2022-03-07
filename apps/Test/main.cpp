// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/HydrogenicSystem.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    HydrogenicSystem hyd;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    
    for (int n=1; n<=10; n++)
    {
        for (int l=0; l < 1; l++)
        {
            int outerA0 = 2*(n+10)*n;
            int steps = 100*outerA0;
            auto [r, psi] = hyd.GetRadialWF(n, l, BohrRadius_v * outerA0 / steps, BohrRadius_v*outerA0, n*100);

            VectorXd prob = 4*Pi_v*(r.array()*r.array()*psi.array()*psi.array());
            ax.Plot(r.data(), prob.data(), r.size());
        }
    }

    matplotlib.RunGUILoop();
#endif

    return 0;
}
