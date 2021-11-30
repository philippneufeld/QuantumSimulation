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
    QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, 100.0});
    system.SetDipoleElementByName("S1_2", "P3_2", 1.0e-37);
    system.AddLaserByName("Probe", "S1_2", "P3_2", 1.0, false);
    system.AddLaserByName("Pump", "S1_2", "P3_2", 1.0, true);
    system.SetDecayByName("P3_2", "S1_2", 0.2);
    
    double dt = 0.1;
    auto rho0 = system.CreateGroundState();

    double pdet = 0.3;
    auto laserDetunings = QSim::CreateLinspaceRow(-0.6, 0.6, 201);
    QSim::TDynamicMatrix<double> detunings(2, laserDetunings.Size());
    QSim::SetRow(detunings, laserDetunings, 0);
    QSim::SetRow(detunings, laserDetunings, 1);

    QSim::PythonMatplotlib matplotlib;
    auto figure1 = matplotlib.MakeFigure();
    auto figure2 = matplotlib.MakeFigure();
    auto ax1 = figure1.AddSubplot();
    auto ax2 = figure2.AddSubplot();

    auto start_ts = std::chrono::high_resolution_clock::now();
    double tmax = 100.0;

    auto traj1 = system.GetTrajectoryNatural(QSim::TStaticColVector<double, 2>({0, 0}), rho0, 0.0, tmax/10.0, dt/100.0);

    auto x_axis1 = laserDetunings;
    auto y_axis1 = QSim::CreateZerosLike(x_axis1);
    for (std::size_t i=0; i< detunings.Cols(); i++)
    {
        auto dets = *(QSim::GetColIteratorBegin(detunings) + i);
        auto rho = system.EvolveNaturalDensityMatrix(
            dets, rho0, 0.0, 0.0, dt, std::size_t(tmax / dt), 0.16 / dt);
        auto traj = system.GetTrajectoryNatural(dets, rho, tmax, tmax * 1.2, dt);
        for (const auto& r: traj.second)
            y_axis1(i) += r.GetPopulation("S1_2");
        y_axis1(i) /= traj.second.size();
        // y_axis1(i) += rho.GetPopulation("S1_2");
    }

    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    auto x_axis2 = traj1.first;
    auto y_axis2 = QSim::CreateZerosLike(x_axis2);
    for (std::size_t i = 0; i < x_axis1.Size(); i++)
    {
        // y_axis2(i) = std::imag(traj1.second[i](0, 1));// * std::exp(std::complex<double>(1.0i*QSim::TwoPi_v*x_axis(i))));
        y_axis2(i) = traj1.second[i].GetPopulation("P3_2");
        // y_axis2(i) = std::real(traj1.second[i].GetPopulation("P3_2") / std::exp(1.0i * std::complex<double>(QSim::TwoPi_v*25.0*i*dt)));
    }
    
    ax1.Plot(x_axis1.Data(), y_axis1.Data(), x_axis1.Size());
    ax2.Plot(x_axis2.Data(), y_axis2.Data(), x_axis2.Size());

    matplotlib.RunGUILoop();

    return 0;
}
