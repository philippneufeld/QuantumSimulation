// Philipp Neufeld, 2021-2022

#include <iostream>
#include <random>
#include <Eigen/Dense>

#include <QSim/Math/LevenbergMarquardt.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

double test1(double x, const Vector4d& beta)
{
    return (((beta[0]/1e8)*x + beta[1])*x + beta[2])*x + beta[3];
}

double randf()
{
    return (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * 2;
}

int main()
{
    Vector4d beta = Vector4d{0.5*1e8, 3, -2, 1};

    constexpr int n = 100;
    double yerrMag = 5;
    VectorXd x = VectorXd::LinSpaced(n, -5, 5);
    VectorXd y = x.unaryExpr([&](double x){ return test1(x, beta) + yerrMag*randf(); });
    VectorXd yerrs = VectorXd::Ones(y.size())*yerrMag/1.5;
    Vector4d beta0 = Vector4d{0, 0, 0, 0};

    LevenbergMarquardt fit;
    auto [params, cov, good] = fit.CurveFitV(test1, x, y, yerrs, beta0);
    Vector4d errs = cov.diagonal().cwiseSqrt();

    VectorXd xfit = VectorXd::LinSpaced(500, -5, 5);
    VectorXd yfit = xfit.unaryExpr([p=params](double x){ return test1(x, p); });
    VectorXd yfit2 = xfit.unaryExpr([p=params+errs](double x){ return test1(x, p); });
    VectorXd yfit3 = xfit.unaryExpr([p=params-errs](double x){ return test1(x, p); });
    
    std::cout << params << std::endl << std::endl;
    std::cout << errs.cwiseQuotient(params) * 100 << std::endl << std::endl;
    std::cout << good << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    ax.Plot(x.eval().data(), y.data(), x.size(), "Data", "x");
    ax.Plot(xfit.data(), yfit.data(), xfit.size(), "Fit");
    ax.Plot(xfit.data(), yfit2.data(), xfit.size(), "Fit upper");
    ax.Plot(xfit.data(), yfit3.data(), xfit.size(), "Fit lower");
    fig.Legend();
    matplotlib.RunGUILoop();
#endif

    return 0;
}
