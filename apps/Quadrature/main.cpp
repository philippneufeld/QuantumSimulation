// Philipp Neufeld, 2021

#include <functional>

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>

class CQuadratureApp : public QSim::CalcApp
{
public:

    template<typename Quad, typename Func>
    double GetQuadError(Func f, double exact, double x0, double x1, std::size_t cnt) 
    {
        return std::abs(exact - Quad{}.Integrate(f, x0, x1, cnt));
    }

    template<typename Quad, typename VT>
    void StoreQuadErrors(std::string name, const QSim::TVector<VT>& ns, 
        std::tuple<std::function<double(double)>, double, double, double>& desc)
    {
        auto errors = QSim::CreateZerosLike(ns);
        for (std::size_t i = 0; i < errors.Size(); i++)
        {
            auto func = std::get<0>(desc);
            auto x0 = std::get<1>(desc);
            auto x1 = std::get<2>(desc);
            auto exact = std::get<3>(desc);
            auto cnt = (~ns)(i);
            errors(i) = GetQuadError<Quad>(func, exact, x0, x1, cnt);
        }
        StoreMatrix(name, errors);
    }

    virtual void DoCalculation() override
    {
        std::size_t minFev = 1;
        std::size_t maxFev = 5000;
        auto ns = QSim::CreateLinspaceRow(1.0, static_cast<double>(maxFev), maxFev);
        StoreMatrix("ns", ns);

        std::tuple<std::function<double(double)>, double, double, double> funcs[] = { 
            {[](double x){ return std::sin(x) * (1 - x + x*x); }, 0.0, 5.0, -15.01989999576955},
            {[](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, -250.0, 250.0, 0.9974535344916211},
            {[](double x){ return std::sin(1 / x); }, 0.025, 1.0, 0.5044592407911533},
            {[](double x){ return abs(x) < 1 ? 1.0 : 0.0; }, -10.0, 10.0, 2.0}
        };

        for (std::size_t i = 0; i < sizeof(funcs)/sizeof(*funcs); i++)
        {
            StoreQuadErrors<QSim::QuadMidpoint>("midpoint" + std::to_string(i), ns, funcs[i]);
            StoreQuadErrors<QSim::QuadTrapezoidal>("trapezoid" + std::to_string(i), ns, funcs[i]);
            StoreQuadErrors<QSim::QuadSimpson>("simpson" + std::to_string(i), ns, funcs[i]);
            StoreQuadErrors<QSim::QuadSimpsonAlt>("simpsonAlt" + std::to_string(i), ns, funcs[i]);
            StoreQuadErrors<QSim::QuadBoole>("boole" + std::to_string(i), ns, funcs[i]);        
            StoreQuadErrors<QSim::QuadWeddle>("weddle" + std::to_string(i), ns, funcs[i]);        
        }
    }

    virtual void Plot() override
    {
        std::cout << QSim::QuadMidpoint{}.Integrate([](double x){ return abs(x) < 1 ? 1.0 : 0.0; }, -1.0, 1.0, 100000) << std::endl;

        std::string titles[] = {
            "f(x) = sin(x)*(1-x+x^2); [0, 5]",
            "f(x) = 1/(pi*(x^2+1)); [-250, 250]",
            "f(x) = sin(1/x); [0.025, 1]",
            "f(x) = Rechteck; [-250, 250]"
        };

        std::string cases[] = {
            "midpoint", "trapezoid", "simpson", "simpsonAlt", "boole", "weddle"};

        QSim::PythonMatplotlib matplotlib;
        auto ns = LoadMatrix("ns");
        for (std::size_t i = 0; i < sizeof(titles)/sizeof(*titles); i++)
        {
            auto fig = matplotlib.CreateFigure();
            auto ax = fig.AddSubplot(); 
            ax.SetTitle(titles[i]);
            for (std::size_t j = 0; j < sizeof(cases)/sizeof(*cases); j++)
            {
                auto name = cases[j] + std::to_string(i);
                ax.Plot(ns.Data(), LoadMatrix(name).Data(), ns.Size(), cases[j]); 
            }
            ax.SetYLog();
            ax.Legend();
            ax.SetXLabel("Function evaluations");
            ax.SetYLabel("Absolute error");
        }
        matplotlib.RunGUILoop();
    }
};


int main(int argc, const char* argv[])
{
    CQuadratureApp app;
    return app.Run(argc, argv);
}
