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
    QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, 80});
    system.SetDipoleElementByName("S1_2", "P3_2", 1.0e-36);
    system.AddLaserByName("laser", "S1_2", "P3_2", 10.0, false);
    system.SetDecayByName("P3_2", "S1_2", 0.1);
    
    double dt = 0.05;
    auto rho0 = system.CreateGroundState();

    auto detunings = QSim::CreateLinspaceRow(-0.7, 0.7, 501);
    // QSim::TStaticColVector<double, 1> detunings({0.0});

    auto start_ts = std::chrono::high_resolution_clock::now();
    
    // auto traj = system.GetTrajectoryNatural(detunings, rho0, 150.0, dt);

    auto x_axis = detunings;
    auto y_axis = QSim::CreateZerosLike(x_axis);

    for (std::size_t i=0; i<detunings.Size(); i++)
    {
        auto dets = *(QSim::GetColIteratorBegin(detunings) + i);
        auto rho = system.EvolveNaturalDensityMatrix(
            dets, rho0, 0.0, 0.0, dt, std::size_t(250.0 / dt), 0.16 / dt);
        y_axis(i) = rho.GetAbsCoeff("S1_2", "P3_2");
    }


    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    QSim::PythonMatplotlib matplotlib;
    auto figure = matplotlib.MakeFigure();
    auto ax = figure.AddSubplot();

    // auto x_axis = traj.first;
    // auto y_axis = QSim::CreateZerosLike(x_axis);
    // for (std::size_t i = 0; i < x_axis.Size(); i++)
    // {
    //     // y_axis(i) = std::imag(traj.second[i](0, 1));// * std::exp(std::complex<double>(1.0i*QSim::TwoPi_v*x_axis(i))));
    //     y_axis(i) = traj.second[i].GetPopulation("P3_2");
    //     // y_axis(i) = std::real(traj.second[i].GetPopulation("P3_2") / std::exp(1.0i * std::complex<double>(QSim::TwoPi_v*25.0*i*dt)));
    // }
    
    ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());

    matplotlib.RunGUILoop();

    return 0;
}
