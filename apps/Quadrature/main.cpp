// Philipp Neufeld, 2021-2022

#ifdef QSIM_PYTHON3
#include <functional>
#include <QSim/Constants.h>
#include <QSim/Math/Quadrature.h>
#include <QSim/Python/Plotting.h>

using namespace QSim;

template<typename Quad>
auto createErrorFunction(std::function<double(double)> func, double a, double b, double exact) 
{
    return [=](double n)
    { 
        double res = Quad{}.Integrate(func, a, b, (std::size_t)n);
        return std::abs(exact - res); 
    };
}

void testIntegrators(const std::string& title, std::function<double(double)> func, 
    double a, double b, double exact)
{
    std::vector<double> ns(5000);
    for (std::size_t i = 0; i < ns.size(); i++)
        ns[i] = i + 1; 
    std::vector<double> errors(ns.size());
    
    QSim::PythonMatplotlib matplotlib;
    
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot(); 
    ax.SetTitle(title + " [" + std::to_string((float)a) + "; " + std::to_string((float)b) + "]");
    ax.SetYLog();
    ax.SetXLabel("Function evaluations");
    ax.SetYLabel("Absolute error");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<TQuadMidpoint<double>>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "midpoint");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<TQuadTrapezoidal<double>>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "trapezoid");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<TQuadSimpson<double>>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "simpson");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<TQuadSimpson38<double>>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "simpson38");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<TQuadBoole<double>>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "boole");

    // adaptive integrator
    auto ns_ad = std::vector<double>(2500);
    auto errs_ad = std::vector<double>(ns_ad.size());
    for (std::size_t j = 0; j < ns_ad.size(); j++)
    {
        double expmin = 2;
        double expmax = 14;
        double rtol = std::pow(10.0, -(j*(expmax - expmin)/ns_ad.size() + expmin));
        auto res = QSim::TQuadAdaptive<double>{}.IntegrateFevs(func, a, b, 250, rtol, std::min(1e-12, rtol), 5);
        errs_ad[j] = std::abs(exact - res.first);
        ns_ad[j] = res.second;
    }
    ax.Plot(ns_ad.data(), errs_ad.data(), ns_ad.size(), "adaptive");
    
    ax.Legend();
}


int main(int argc, const char* argv[])
{
    QSim::PythonMatplotlib matplotlib;

    testIntegrators("sin(x)*(1-x+x^2)", 
        [](double x){ return std::sin(x) * (1 - x + x*x); }, 
        0.0, 5.0, -15.01989999576955);

    testIntegrators("1/(pi*(x^2+1))", 
            [](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, 
            -2, 2, 0.7048327646991335);

    testIntegrators("1/(pi*(x^2+1))", 
            [](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, 
            -250.0, 250.0, 0.9974535344916211);

    testIntegrators("sin(1/x)", 
            [](double x){ return std::sin(1 / x); }, 
            0.025, 1.0, 0.5044592407911533);
    
    matplotlib.RunGUILoop();
}

#else
#include <iostream>
int main(int argc, const char* argv[])
{
    std::cout << "QSIM_PYTHON3 macro is not set. Aborting..." << std::endl;
}
#endif
