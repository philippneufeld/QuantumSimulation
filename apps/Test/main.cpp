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
    StarkMap starkMap(hyd, 28, 0, 0, 23, 32, 20);
    // StarkMap starkMap(hyd, 2, 0, 0, 1, 3, 3);

    // for (auto state: starkMap.GetBasis())
    // {
    //     auto [n, l, m] = state;
    //     std::cout << "n=" << n << " l=" << l << " m=" << m << std::endl;
    // }
    // std::cout << starkMap.GetBasis().size() << " states in total" << std::endl;

    // std::cout << starkMap.GetDipoleOperator() << std::endl;

    int cnt = starkMap.GetBasis().size();

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
        VectorXd energies = starkMap.GetEnergies(el);
        time += (std::chrono::high_resolution_clock::now() - ts).count() / 1e6;
        VectorXd eField = VectorXd::Ones(cnt) * el;
        ax.Plot(eField.data(), energies.data(), cnt, "", ".k");
    }
    
    std::cout << time / N << "ms" << std::endl;

    matplotlib.RunGUILoop();
#endif

    return 0;
}
