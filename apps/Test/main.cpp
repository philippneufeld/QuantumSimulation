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

    QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, 1000.0});
    system.SetDipoleElementByName("S1_2", "P3_2", 1.0e-36);
    system.AddLaserByName("Probe", "S1_2", "P3_2", 1.0, false);
    system.AddLaserByName("Pump", "S1_2", "P3_2", 1.0, true);
    system.SetDecayByName("P3_2", "S1_2", 0.1);
    system.SetMass(1.0e-31);
    
    double dt = 1e-2;
    double tmax = 20;
    auto rho0 = system.CreateGroundState();

    auto laserDetunings = QSim::CreateLinspaceRow(-3.0, 3.0, 201);
    QSim::TDynamicMatrix<double> detunings(2, laserDetunings.Size());
    QSim::SetRow(detunings, laserDetunings, 0);
    QSim::SetRow(detunings, laserDetunings, 1);

    auto start_ts = std::chrono::high_resolution_clock::now();

    auto func = [&](auto dets)
    { 
        auto rho = system.GetDensityMatrixAv(
            dets, rho0, 0.0, tmax, 0.1*tmax, dt);
        return rho.GetAbsCoeff("S1_2", "P3_2");
    }; 
    auto absCoeffs = pool.Map(func, 
        QSim::GetColIteratorBegin(detunings), 
        QSim::GetColIteratorEnd(detunings));

    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    QSim::PythonMatplotlib matplotlib;

    auto figure = matplotlib.MakeFigure();
    auto ax = figure.AddSubplot();
    ax.Plot(laserDetunings.Data(), absCoeffs.data(), laserDetunings.Size());

    matplotlib.RunGUILoop();

    return 0;
}
