// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/HydrogenicSystem.h>
#include <QSim/Rydberg/StarkMap.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Math/Wigner.h>
#include <chrono>

using namespace QSim;
using namespace Eigen;

double error(double exact, double val)
{
    return std::abs((exact-val)/exact);
}

int main(int argc, const char* argv[])
{
    HydrogenicSystem hyd;

    // std::cout << hyd.GetDipRadialME(23, 0, 0.5, 23, 1, 0.5) << std::endl;
    // std::cout << hyd.GetDipAngularME2(0, 0.5, 0.5, 1, 0.5, 0.5) << std::endl;
    // std::cout << hyd.GetDipoleME2(23, 0, 0.5, 0.5, 23, 1, 0.5, 0.5) << std::endl;
    // return 0;


    auto [r1, p1] = hyd.GetRadialWF(28, 8, 2.5, 0.1*BohrRadius_v, 28*(28+15)*4*BohrRadius_v, 2500);
    auto [r2, p2] = hyd.GetRadialWF(28, 3, 2.5, 0.1*BohrRadius_v, 28*(28+15)*4*BohrRadius_v, 2500);


    std::cout << hyd.GetDipRadialME(28, 2, 2.5, 28, 3, 2.5) << std::endl;
    std::cout << hyd.GetDipAngularME2(2, 2.5, 0.5, 3, 2.5, 0.5) << std::endl;
    std::cout << hyd.GetDipoleME2(28, 2, 2.5, 0.5, 28, 3, 2.5, 0.5) << std::endl;

    /*PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    r1 = r1 / BohrRadius_v;
    r2 = r2 / BohrRadius_v;*/

    /*ax.Plot(r1.data(), (r1.array().square()*p1.array().square()).eval().data(), r1.size());
    ax.Plot(r2.data(), (r2.array().square()*p2.array().square()).eval().data(), r2.size());*/

    //ax.Plot(r1.data(), (r1.array().pow(1.5)*p1.array()).square().eval().data(), r1.size());

    //matplotlib.RunGUILoop();


    StarkMap starkMap(hyd, 28, 0, 0.5, 0.5, 23, 31, 20);
    // StarkMap starkMap(hyd, 28, 0, 0.5, 0.5, 28, 28, 3);

    // for (auto state: starkMap.GetBasis())
    // {
    //     auto [n, l, j, mj] = state;
    //     std::cout << "n=" << n << " l=" << l << " j=" << j << " mj=" << mj << std::endl;
    // }
    // std::cout << starkMap.GetBasis().size() << " states in total" << std::endl;

    // std::cout << starkMap.GetDipoleOperator() << std::endl;

    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    double time = 0;

    int N = 180;
    for (int i=0; i < N; i++)
    {
        double el = i * 60000 / N;
        auto ts = std::chrono::high_resolution_clock::now();
        VectorXd energies = starkMap.GetEnergies(el) / (PlanckConstant_v*SpeedOfLight_v) / 100;
        time += (std::chrono::high_resolution_clock::now() - ts).count() / 1e6;
        VectorXd eField = VectorXd::Ones(cnt) * el / 100;
        ax.Plot(eField.data(), energies.data(), cnt, "", ".C0");
    }
    
    std::cout << time / N << "ms" << std::endl;

    matplotlib.RunGUILoop();
#endif

    return 0;
}
