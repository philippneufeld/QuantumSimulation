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

    // std::cout << hyd.GetDipoleMERad(28, 2, 2.5, 28, 3, 2.5) << std::endl;
    // std::cout << hyd.GetDipoleMEAng2(2, 2.5, 0.5, 3, 2.5, 0.5) << std::endl;
    // std::cout << hyd.GetDipoleME2(28, 2, 2.5, 0.5, 28, 3, 2.5, 0.5) << std::endl;

    StarkMap starkMap(hyd, 28, 0, 0.5, 0.5, 23, 31, 20);

    // for (auto state: starkMap.GetBasis())
    // {
    //     auto [n, l, j, mj] = state;
    //     std::cout << "n=" << n << " l=" << l << " j=" << j << " mj=" << mj << std::endl;
    // }

    // std::cout << starkMap.GetDipoleOperator() << std::endl;

    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    double time = 0;
    int N = 600;
    for (int i=0; i < N; i++)
    {
        double el = i * 60000 / N;
        auto ts = std::chrono::high_resolution_clock::now();
        VectorXd energies = starkMap.GetEnergies(el) / (PlanckConstant_v*SpeedOfLight_v) / 100;
        time += (std::chrono::high_resolution_clock::now() - ts).count() / 1e6;
        VectorXd eField = VectorXd::Ones(cnt) * el / 100;
        
        ax.Plot(eField.data(), energies.data(), cnt, "", ".C0");
        ax.SetYLimits(-150.0, -130.0);
        ax.SetXLimits(0, 600);
        ax.SetXLabel("Electric Field (V/cm)");
        ax.SetYLabel("Energy $\\frac{E}{hc}$ (cm$^{-1}$)");
    }
    
    std::cout << time / N << "ms per diagonalization step" << std::endl;

    matplotlib.RunGUILoop();
#endif

    return 0;
}
