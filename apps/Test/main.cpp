// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    
    // auto x = VectorXd::LinSpaced(11, 0.0, 10.0);
    auto y = VectorXd::NullaryExpr(11, [](int t){return double(t*t);});



    std::cout << y.eval().sum() << std::endl;

/*#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size());
    matplotlib.RunGUILoop();
#endif*/

    return 0;
}
