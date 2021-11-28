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

    QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, 25});
    system.SetDipoleElementByName("S1_2", "P3_2", 1.0e-36);
    system.AddLaserByName("laser", "S1_2", "P3_2", 10.0, false);
    system.SetDecayByName("P3_2", "S1_2", 0.1);
    
    QSim::TStaticColVector<double, 1> detunings({0.0});

    auto start_ts = std::chrono::high_resolution_clock::now();
    
    double dt = 0.001;
    auto rho0 = system.CreateGroundState();
    auto traj = system.GetTrajectoryNatural(detunings, rho0, dt, 100);

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
        // y_axis(i) = std::real(traj.second[i].GetPopulation("P3_2") / std::exp(1.0i * std::complex<double>(QSim::TwoPi_v*25.0*i*dt)));
    }
    
    ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());

    matplotlib.RunGUILoop();

    return 0;
}
