// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Math/Quad.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    
    auto x = VectorXd::LinSpaced(8, 0.0, 1.0);
    auto y = x.unaryExpr([](double t){return t*t*t*t;});

    TQuadrature<QuadSimpsonPolicy> quad;

    // std::cout << y.eval()(Eigen::seq(0, y.size() - 1, 10)) << std::endl;

    std::cout << quad.Integrate(y, x[1] - x[0]) << std::endl;
    std::cout << quad.Integrate([](double t){return t*t*t*t;}, 0.0, 1.0, 8) << std::endl;

/*#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size());
    matplotlib.RunGUILoop();
#endif*/

    return 0;
}
