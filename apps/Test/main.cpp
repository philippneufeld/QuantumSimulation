// Philipp Neufeld, 2021

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <QSim/Python/Plotting.h>
#include <QSim/Util/Argparse.h>
#include <QSim/StaticQSys.h>
#include <QSim/Util/ThreadPool.h>

int main(int argc, const char* argv[])
{
    QSim::ThreadPool pool;

    QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, 1.0e8});
    system.SetDipoleElementByName("S1_2", "P3_2", 1.0e-35);
    system.AddLaserByName("Probe", "S1_2", "P3_2", 20, false);
    system.AddLaserByName("Pump", "S1_2", "P3_2", 2, true);
    system.SetDecayByName("P3_2", "S1_2", 1);
    system.SetMass(1.0e-23);

    constexpr auto rabi = QSim::GetRabiFrequencyFromIntensity(1e-35, 20);
    
    double dt = 3e-3;
    double tmax = 10;
    auto rho0 = system.CreateGroundState();

    auto laserDetunings = QSim::CreateLinspaceRow(-25.0, 25.0, 201);
    QSim::TDynamicMatrix<double> detunings(2, laserDetunings.Size());
    QSim::SetRow(detunings, laserDetunings, 0);
    QSim::SetRow(detunings, laserDetunings, 1);

    auto start_ts = std::chrono::high_resolution_clock::now();

    // auto trace = system.GetNaturalTrajectory(QSim::TStaticColVector<double, 2>({10.0, 60.0}), rho0, 0.0, 0, tmax, dt);
    // auto y_axis = QSim::CreateZerosLike(trace.first);
    // for (std::size_t i = 0; i < y_axis.Size(); i++)
    //     y_axis[i] = trace.second[i].GetPopulation("P3_2");

    auto func = [&](auto dets)
    { 
        auto rho = system.GetDensityMatrixAv(
            dets, rho0, 0.0, tmax, 0.25*tmax, dt);
        return rho.GetPopulation("P3_2");
    }; 
    auto absCoeffs = pool.Map(func, 
        QSim::GetColIteratorBegin(detunings), 
        QSim::GetColIteratorEnd(detunings));

    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    QSim::PythonMatplotlib matplotlib;

    auto figure = matplotlib.MakeFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(laserDetunings.Data(), absCoeffs.data(), laserDetunings.Size());
    // ax.Plot(trace.first.Data(), y_axis.Data(), y_axis.Size());

    matplotlib.RunGUILoop();

    return 0;
}
