// Philipp Neufeld, 2021-2022

#include <iostream>

#ifdef QSIM_PYTHON3
#include <functional>
#include <QSim/Constants.h>
#include <QSim/Math/Quad.h>
#include <QSim/Math/QuadAdapt.h>
#include <QSim/Python/Plotting.h>

using namespace QSim;

template<typename Quad>
auto createErrorFunction(std::function<double(double)> func, double a, double b, double exact) 
{
    return [=](double n)
    { 
        double res = TQuadrature<Quad>::Integrate(func, a, b, (std::size_t)n);
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
    
    PythonMatplotlib matplotlib;
    
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot(); 
    ax.SetTitle(title + " [" + std::to_string((float)a) + "; " + std::to_string((float)b) + "]");
    ax.SetYLog();
    ax.SetXLabel("Function evaluations");
    ax.SetYLabel("Absolute error");


    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<QuadBoolePolicy>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "boole");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<QuadSimpsonPolicy>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "simpson");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<QuadMidpointPolicy>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "midpoint");
    
    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<QuadTrapezoidalPolicy>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "trapezoid");

    for (std::size_t i = 0; i < ns.size(); i++)
        errors[i] = std::invoke(createErrorFunction<QuadFixedAdaptivePolicy>(func, a, b, exact), ns[i]);
    ax.Plot(ns.data(), errors.data(), ns.size(), "ad2");

    // adaptive integrator
    std::size_t cnt = 2500;
    auto ns_ad1 = std::vector<double>(cnt);
    auto errs_ad1 = std::vector<double>(ns_ad1.size());
    auto ns_ad2 = std::vector<double>(cnt);
    auto errs_ad2 = std::vector<double>(ns_ad1.size());
    for (std::size_t j = 0; j < cnt; j++)
    {
        double expmin = 2;
        double expmax = 14;
        double rtol = std::pow(10.0, -(j*(expmax - expmin)/cnt + expmin));
        
        auto [res1, fevs1] = TQuadrature<QuadAdaptivePolicy>::IntegrateFevs(func, a, b, 250, rtol, std::min(1e-12, rtol), 10);
        errs_ad1[j] = std::abs(exact - res1);
        ns_ad1[j] = fevs1;
    }
    ax.Plot(ns_ad1.data(), errs_ad1.data(), ns_ad1.size(), "adaptive");
    
    ax.Legend();
}


int main(int argc, const char* argv[])
{
    // testIntegrators("sin(x)*(1-x+x^2)", 
    //     [](double x){ return std::sin(x) * (1 - x + x*x); }, 
    //     0.0, 5.0, -15.01989999576955);
    // testIntegrators("1/(pi*(x^2+1))", 
    //         [](double x){ return 1 / (Pi_v * (x*x + 1)); }, 
    //         -2, 2, 0.7048327646991335);
    testIntegrators("1/(pi*(x^2+1))", 
            [](double x){ return 1 / (Pi_v * (x*x + 1)); }, 
            -250.0, 250.0, 0.9974535344916211);
    // testIntegrators("sin(1/x)", 
    //         [](double x){ return std::sin(1 / x); }, 
    //         0.025, 1.0, 0.5044592407911533);


    PythonMatplotlib matplotlib;

    std::vector<double> samples;
    constexpr double x0 = 0.0;
    auto func = [&](double x){ return 1 / (Pi_v * ((x-x0)*(x-x0) + 1)); };
    auto func_tracked = [&](double x){ samples.push_back(x); return func(x); };
    std::size_t n = 100;

    TQuadrature<QuadFixedAdaptivePolicy> quadrature;
    auto [I2, fevs2] = quadrature.IntegrateFevs(func_tracked, -50.0, 50.0, n);

    auto fig = matplotlib.CreateFigure();
    auto ax1 = fig.AddSubplot(2, 1, 1);
    auto ax2 = fig.AddSubplot(2, 1, 2);

    auto samplesy = samples;
    for (auto& s: samplesy) s = func(s);
    ax1.Plot(samples.data(), samplesy.data(), samples.size(), "", ".");
    ax2.Plot(samples.data(), std::vector<double>(samples.size(), 0.0).data(), samples.size(), "", ".");

    matplotlib.RunGUILoop();
}

#else
int main(int argc, const char* argv[])
{
    std::cout << "QSIM_PYTHON3 macro is not set. Aborting..." << std::endl;
}
#endif
