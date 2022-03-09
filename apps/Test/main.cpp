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

int main(int argc, const char* argv[])
{
    HydrogenicSystem hyd;

    double ang = 1 / std::sqrt(3);
    double rad = 128.0/81.0*std::sqrt(2.0/3.0)*BohrRadius_v;
    double exact1 = std::sqrt(2)*128/243*BohrRadius_v;
    std::cout << error(ang, hyd.GetDipAngularME(0,0,1,0)) << std::endl;
    std::cout << error(rad, hyd.GetDipRadialME(1,0,2,1)) << std::endl;
    std::cout << error(exact1, hyd.GetDipoleME(1,0,0,2,1,0)) << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    
    int outerA0 = 100;
    int steps = 100*outerA0;
    auto [r1, psi1] = hyd.GetRadialWF(1, 0, BohrRadius_v * outerA0 / steps, BohrRadius_v*outerA0, 1000);
    auto [r2, psi2] = hyd.GetRadialWF(2, 1, BohrRadius_v * outerA0 / steps, BohrRadius_v*outerA0, 1000);

    VectorXd theo1 = r1.unaryExpr([](double r){ return 2/std::pow(BohrRadius_v, 1.5)*std::exp(-r/BohrRadius_v); });
    VectorXd theo2 = r2.unaryExpr([](double r){ return 0.5/std::sqrt(6)/std::pow(BohrRadius_v, 2.5)*r*std::exp(-r/(2*BohrRadius_v)); });
    ax.Plot(r1.data(), psi1.data(), r1.size());
    ax.Plot(r1.data(), theo1.data(), r1.size());
    ax.Plot(r2.data(), psi2.data(), r2.size());
    ax.Plot(r2.data(), theo2.data(), r2.size());


    /*for (int n=1; n<=2; n++)
    {
        for (int l=0; l<n; l++)
        {   
            int outerA0 = 2*(n+10)*n;
            int steps = 100*outerA0;
            auto [r, psi] = hyd.GetRadialWF(n, l, BohrRadius_v * outerA0 / steps, BohrRadius_v*outerA0, n*100);

            VectorXd prob = 4*Pi_v*(r.array()*r.array()*psi.array()*psi.array());
            ax.Plot(r.data(), prob.data(), r.size());
        }
    }*/

    matplotlib.RunGUILoop();
#endif

    return 0;
}
