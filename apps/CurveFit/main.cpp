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

std::pair<Eigen::VectorXd, Eigen::VectorXd> GenerateTestData(ThreadPoolExecutor& pool, std::size_t n, double r)
{
    constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
    constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
    
    // Create system
    TNLevelSystemQM<2> system;
    system.SetLevel(0, 0.0);
    system.SetLevel(1, QSim::SpeedOfLight_v / 780.241e-9);
    system.SetDipoleElement(0, 1, dip);
    system.AddLaser(0, 1, intProbe, false);
    system.SetDecay(1, 0, 6.065e6);

    DopplerIntegrator doppler;
    doppler.SetMass(1.44316060e-25);

    auto func = [&](auto dets)
    {
        return doppler.Integrate([&](double vel)
        { 
            auto rho = system.GetDensityMatrixSS(dets, vel);
            return std::imag(rho(0, 1));
        }) + randf() * r;
    };

    Eigen::RowVectorXd detunings = Eigen::RowVectorXd::LinSpaced(n, -2e9, 2e9);
    Eigen::VectorXd absCoeffs(detunings.cols());
    auto genDetuning = [&detunings](){static int i=0; return detunings.col(i++).eval(); };
    pool.MapG(func, absCoeffs, genDetuning, detunings.cols());
    
    return std::make_pair(detunings, absCoeffs);
}


class Fitter
{
public:
    Fitter() 
    {
        m_doppler.SetMass(1.44316060e-25);
    }

    static TNLevelSystemQM<2> CreateSystem()
    {
        double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        
        TNLevelSystemQM<2> system;
        system.SetLevel(0, 0.0);
        system.SetLevel(1, QSim::SpeedOfLight_v / 780.241e-9);
        system.SetDipoleElement(0, 1, dip);
        system.AddLaser(0, 1, intProbe, false);
        system.SetDecay(1, 0, 6.065e6);

        return system;
    }

    double operator()(double x, double x0, double gamma)
    {
        Eigen::Matrix<double, 1, 1> dets;
        dets << (x - x0);
        system.SetDecay(1, 0, gamma);
        s_cnt++;
        return m_doppler.Integrate([&](double vel)
        { 
            auto rho = system.GetDensityMatrixSS(dets, vel);
            return std::imag(rho(0, 1));
        });
    }

    static std::atomic<std::size_t> s_cnt;

private:
    DopplerIntegrator m_doppler;
    thread_local static TNLevelSystemQM<2> system;
};

std::atomic<std::size_t> Fitter::s_cnt = 0;
thread_local TNLevelSystemQM<2> Fitter::system = Fitter::CreateSystem();


int main(int argc, const char* argv[])
{
    ThreadPoolExecutor pool;
    auto [x, y] = GenerateTestData(pool, 501, 1e-3);

    double x0 = 10e7;
    double gamma = 1e6;
    // double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
    Fitter fitter;
    CurveFit(pool, fitter, x, y, x0, gamma);

    std::cout << "gamma: " << gamma / 1e6 << " MHz" << std::endl;
    // std::cout << "dip: " << dip / (QSim::ElementaryCharge_v * QSim::BohrRadius_v) << " e*a0" << std::endl;
    std::cout << "x0: " << x0 / 1e6 << " MHz" << std::endl;

    std::cout << "fevs: " << fitter.s_cnt << std::endl;

    Eigen::VectorXd yfit = Eigen::VectorXd::Zero(x.size());
    pool.Map([&](double x) { return fitter(x, x0, gamma); }, 
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
