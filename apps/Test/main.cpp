// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/HydrogenicSystem.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Math/Wigner.h>

using namespace QSim;
using namespace Eigen;

double error(double exact, double val)
{
    return std::abs((exact-val)/exact);
}

double sq(double x) { return x == 0 ? 0 : x*x*x/std::abs(x); }

int main(int argc, const char* argv[])
{
    HydrogenicSystem hyd;


    std::cout << sq(ClebshGordan(1, 1, 2, 1, -1, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 2, 0, 0, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 2, -1, 1, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 1, 1, -1, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 1, 0, 0, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 1, -1, 1, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 0, 1, -1, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 0, 0, 0, 0)) << std::endl;
    std::cout << sq(ClebshGordan(1, 1, 0, -1, 1, 0)) << std::endl;
    
    return 0;


#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    
    /*auto dip1 = hyd.GetRadialMatrixElementLinear(1, 0, 1, 0, ax);
    std::cout << error(dip1, 1.5*BohrRadius_v) << std::endl;

    auto dip2 = hyd.GetRadialMatrixElementLinear(3, 0, 3, 0, ax);
    std::cout << error(dip2, 27.0/2*BohrRadius_v) << std::endl;*/

    double exact = std::pow(2.0, 5.5)/81*BohrRadius_v;
    auto dip3 = hyd.GetRadialMatrixElement(2, 0, 1, 0);
    std::cout << error(exact, dip3) << std::endl;

    for (int n=1; n<=10; n++)
    {
        int l = 0;
        int outerA0 = 2*(n+10)*n;
        int steps = 100*outerA0;
        auto [r, psi] = hyd.GetRadialWF(n, l, BohrRadius_v * outerA0 / steps, BohrRadius_v*outerA0, n*100);

        VectorXd prob = 4*Pi_v*(r.array()*r.array()*psi.array()*psi.array());
        ax.Plot(r.data(), prob.data(), r.size());
    }

    matplotlib.RunGUILoop();
#endif

    return 0;
}
