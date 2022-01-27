// Philipp Neufeld, 2021-2022

#include <chrono>

#include <QSim/Math/CurveFit.h>
#include <QSim/Executor/ThreadPool.h>
#include <QSim/Util/CLIProgressBar.h>

#include <QSim/NLevel/NLevelSystemQM.h>
#include <QSim/NLevel/Doppler.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace std::chrono_literals;

double func(double x, double a, double g, double x0)
{
    auto y = x - x0;
    return a / (y*y + g*g/4);
}

float randf()
{
    return 2 * (static_cast<double>(rand()) / RAND_MAX) - 1;
}




class Fitter
{
public:
    Fitter() 
    {
        m_doppler.SetMass(1.44316060e-25);
    }

    static TNLevelSystemQM<3> CreateSystem()
    {
        // calculate parameters
        constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
        constexpr double intProbe = GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = GetIntensityFromRabiFrequency(dip, 10.0e6);

        // Create system
        TNLevelSystemQM<3> system;
        system.SetLevel(0, -4.271e9);
        system.SetLevel(1, 2.563e9);
        system.SetLevel(2, SpeedOfLight_v / 780.241e-9);
        system.SetDipoleElement(0, 2, dip);
        system.SetDipoleElement(1, 2, dip);
        system.AddLaser(0, 2, intProbe, false);
        system.AddLaser(1, 2, intPump, false);
        system.SetDecay(2, 0, 3.0/8.0 * 6.065e6);
        system.SetDecay(2, 1, 5.0/8.0 * 6.065e6);

        return system;
    }

    static std::pair<Eigen::VectorXd, Eigen::VectorXd> GenerateTestData(ThreadPoolExecutor& pool, std::size_t n, double r)
    {
        constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
        constexpr double intProbe = GetIntensityFromRabiFrequency(dip, 3.5e6);
        
        // Create system
        TNLevelSystemQM<3> system = CreateSystem();
        
        DopplerIntegrator doppler;
        doppler.SetMass(1.44316060e-25);

        auto func = [&](double det)
        {
            return doppler.Integrate([&](double vel)
            { 
                Eigen::VectorXd dets(2);
                dets << det, 0;
                auto rho = system.GetDensityMatrixSS(dets, vel);
                return std::imag(rho(0, 2));
            }) + randf() * r;
        };

        Eigen::RowVectorXd detunings = Eigen::RowVectorXd::LinSpaced(n, -2e8, 2e8);
        Eigen::VectorXd absCoeffs(detunings.cols());
        auto genDetuning = [&detunings](){static int i=0; return detunings[i++]; };
        pool.MapG(func, absCoeffs, genDetuning, detunings.cols());
        
        return std::make_pair(detunings, absCoeffs);
    }

    double operator()(double x, double x0, double gamma1, double gamma2)
    {
        Eigen::Matrix<double, 1, 1> dets;
        dets << (x - x0);
        system.SetDecay(1, 0, gamma1);
        system.SetDecay(2, 0, gamma2);
        s_cnt++;
        return m_doppler.Integrate([&](double vel)
        { 
            auto rho = system.GetDensityMatrixSS(dets, vel);
            return std::imag(rho(0, 2));
        });
    }

    static std::atomic<std::size_t> s_cnt;

private:
    DopplerIntegrator m_doppler;
    thread_local static TNLevelSystemQM<3> system;
};

std::atomic<std::size_t> Fitter::s_cnt = 0;
thread_local TNLevelSystemQM<3> Fitter::system = Fitter::CreateSystem();


int main(int argc, const char* argv[])
{
    ThreadPoolExecutor pool;
    auto [x, y] = Fitter::GenerateTestData(pool, 501, 1e-3);

    double x0 = 5e6;
    double gamma1 = 1e6;
    double gamma2 = 1e6;
    Fitter fitter;
    CurveFit(pool, fitter, x, y, x0, gamma1, gamma2);

    std::cout << "gamma1: " << gamma1 / 1e6 << " MHz" << std::endl;
    std::cout << "gamma2: " << gamma1 / 1e6 << " MHz" << std::endl;
    std::cout << "x0: " << x0 / 1e6 << " MHz" << std::endl;

    std::cout << "fevs: " << fitter.s_cnt << std::endl;

    Eigen::VectorXd yfit = Eigen::VectorXd::Zero(x.size());
    pool.Map([&](double x) { return fitter(x, x0, gamma1, gamma2); }, 
        yfit, x.data(), x.data() + x.size());

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto figure = matplotlib.CreateFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(x.data(), y.data(), x.size(), "Data", "-");
    ax.Plot(x.data(), yfit.data(), x.size(), "Fit", "-");
    matplotlib.RunGUILoop();
#else
    std::cout << "Plotting disabled" << std::endl;
#endif

    return 0;
}
