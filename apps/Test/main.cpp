// Philipp Neufeld, 2021

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <QSim/Python/Plotting.h>
#include <QSim/Util/Argparse.h>
#include <QSim/StaticQSys.h>

int main(int argc, const char* argv[])
{

    std::map<std::string, double> levels;
    levels["S1_2"] = 0;
    levels["P3_2"] = 25;
    double mass = 1;
    
    QSim::TStaticQSys<2> system(levels, mass);
    system.SetDipoleMatrixElement("S1_2", "P3_2", 1.0e-37);
    // system.AddDecay("P3_2", "S1_2", 0.1);
    
    QSim::TStaticColVector<double, 1> intensity({10});
    QSim::TStaticColVector<double, 1> frequency({25.0});

    auto start_ts = std::chrono::high_resolution_clock::now();
    
    double dt = 0.005;
    auto rho0 = system.MakeGroundState();
    auto traj = system.GetTrajectoryNatural(frequency, intensity, rho0, dt, 100);

    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.MakeFigure();
    auto ax = figure.AddSubplot();

    auto x_axis = traj.first;
    auto y_axis = QSim::CreateZerosLike(x_axis);

    for (std::size_t i = 0; i < x_axis.Size(); i++)
    {
        // y_axis(i) = std::imag(traj[i](0, 1) * std::exp(std::complex<double>(1.0i*QSim::TwoPi_v*x_axis(i))));
        y_axis(i) = traj.second[i].GetPopulation("P3_2");
    }
    
    ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());

    matplotlib.RunGUILoop();

    return 0;
}
