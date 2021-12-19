// Philipp Neufeld, 2021

#include <functional>

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>
#include <QSim/Math/Quadrature.h>

#include <QSim/Util/ConstList.h>

class CQuadratureApp : public QSim::CalcApp
{
    using MyParent = QSim::CalcApp;
public:

    CQuadratureApp() : MyParent()
    {
        auto res1 = QSim::TQuadAdaptive<double>{}.IntegrateFevs(
            [](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, -250, 250, 500, 1e-8, 1e-10, 10);
        std::cout << std::abs(0.9974535344916211 - res1.first) << " (" << res1.second << " fevs)" << std::endl;

        m_funcs.push_back({"sin(x)*(1-x+x^2)", 
            [](double x){ return std::sin(x) * (1 - x + x*x); }, 
            0.0, 5.0, -15.01989999576955});
        m_funcs.push_back({"1/(pi*(x^2+1))", 
            [](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, 
            -2, 2, 0.7048327646991335});
        m_funcs.push_back({"1/(pi*(x^2+1))", 
            [](double x){ return 1 / (QSim::Pi_v * (x*x + 1)); }, 
            -250.0, 250.0, 0.9974535344916211});
        m_funcs.push_back({"sin(1/x)", 
            [](double x){ return std::sin(1 / x); }, 
            0.025, 1.0, 0.5044592407911533});
    }

    template<typename Quad, typename VT>
    void StoreQuadErrors(std::string name, const QSim::TVector<VT>& ns, 
        std::tuple<std::string, std::function<double(double)>, double, double, double>& desc)
    {
        auto errors = QSim::CreateZerosLike(ns);
        for (std::size_t i = 0; i < errors.Size(); i++)
        {
            auto func = std::get<1>(desc);
            auto x0 = std::get<2>(desc);
            auto x1 = std::get<3>(desc);
            auto exact = std::get<4>(desc);
            auto cnt = (~ns)(i);
            errors(i) = std::abs(exact - Quad{}.Integrate(func, x0, x1, cnt));
        }
        StoreMatrix(name + "_n", ns);
        StoreMatrix(name + "_err", errors);
    }

    virtual void DoCalculation() override
    {
        std::size_t minFev = 1;
        std::size_t maxFev = 5000;
        auto ns = QSim::CreateLinspaceRow(1.0, static_cast<double>(maxFev), maxFev);
        StoreMatrix("ns", ns);

        for (std::size_t i = 0; i < m_funcs.size(); i++)
        {
            StoreQuadErrors<QSim::TQuadMidpoint<double>>("midpoint" + std::to_string(i), ns, m_funcs[i]);
            StoreQuadErrors<QSim::TQuadTrapezoidal<double>>("trapezoid" + std::to_string(i), ns, m_funcs[i]);
            StoreQuadErrors<QSim::TQuadSimpson<double>>("simpson" + std::to_string(i), ns, m_funcs[i]);
            StoreQuadErrors<QSim::TQuadSimpson38<double>>("simpson38" + std::to_string(i), ns, m_funcs[i]);
            StoreQuadErrors<QSim::TQuadBoole<double>>("boole" + std::to_string(i), ns, m_funcs[i]);        
        }

        // adaptive integrators
        for (std::size_t i = 0; i < m_funcs.size(); i++)
        {
            auto ns_ad1 = QSim::CreateZeros(2500);
            auto errs_ad1 = QSim::CreateZerosLike(ns_ad1);
            for (std::size_t j = 0; j < ns_ad1.Size(); j++)
            {
                auto func = std::get<1>(m_funcs[i]);
                auto x0 = std::get<2>(m_funcs[i]);
                auto x1 = std::get<3>(m_funcs[i]);
                auto exact = std::get<4>(m_funcs[i]);
                double expmin = 2;
                double expmax = 14;
                double rtol = std::pow(10.0, -(j*(expmax - expmin)/ns_ad1.Size() + expmin));
                auto res = QSim::TQuadAdaptive<double>{}.IntegrateFevs(func, x0, x1, 250, rtol, std::min(1e-12, rtol), 5);
                errs_ad1[j] = std::abs(exact - res.first);
                ns_ad1[j] = res.second;
            }
            StoreMatrix("adaptive" + std::to_string(i) + "_n", ns_ad1);
            StoreMatrix("adaptive" + std::to_string(i) + "_err", errs_ad1);
        }
    }

    virtual void Plot() override
    {
        std::string cases[] = { 
            "midpoint", "trapezoid", "simpson", 
            "simpson38", "boole", "adaptive"};

        QSim::PythonMatplotlib matplotlib;
        
        for (std::size_t i = 0; i < m_funcs.size(); i++)
        {
            auto fig = matplotlib.CreateFigure();
            auto ax = fig.AddSubplot(); 
            ax.SetTitle(
                std::get<0>(m_funcs[i]) + " [" + std::to_string(std::get<2>(m_funcs[i])) 
                + "; " + std::to_string(std::get<3>(m_funcs[i])) + "]");
            for (std::size_t j = 0; j < sizeof(cases)/sizeof(*cases); j++)
            {
                auto name = cases[j] + std::to_string(i);
                auto n = LoadMatrix(name + "_n");
                auto err = LoadMatrix(name + "_err");
                ax.Plot(n.Data(), err.Data(), n.Size(), cases[j]); 
            }

            ax.SetYLog();
            ax.Legend();
            ax.SetXLabel("Function evaluations");
            ax.SetYLabel("Absolute error");
        }
        matplotlib.RunGUILoop();
    }

private:
    std::vector<std::tuple<std::string, std::function<double(double)>, double, double, double>> m_funcs;
};


int main(int argc, const char* argv[])
{
    CQuadratureApp app;
    return app.Run(argc, argv);
}
