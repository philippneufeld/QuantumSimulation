// Philipp Neufeld, 2021

#include <iostream>
#include <chrono>
#include <fstream>

#include <QSim/Matrix.h>
#include <QSim/NLevelSystemStatic.h>
#include <QSim/Doppler.h>
#include <QSim/ThreadPool.h>

int main()
{
    QSim::ThreadPool pool;

    // Generate detuning axis
    auto detunings = QSim::CreateLinspaceRow(-1.0e9, 1.0e9, 501);

    // setup Rb87 parameters
    std::map<std::string, double> levels;
    levels["S1_2"] = 0;
    levels["P3_2"] = QSim::SpeedOfLight_v / 780.241e-9;
    double mass = 1.44316060e-25;
    double temperature = 300.0;

    // create system object
    QSim::TStaticNLevelSystem<2> system(levels, mass);
    system.AddTransition("S1_2", "P3_2", 3.5e6);
    system.AddDecay("P3_2", "S1_2", 6.065e6);
    system.SetTemperature(temperature);

    auto start_ts = std::chrono::high_resolution_clock::now();

    // auto absCoeffs = QSim::CreateZerosLike(detunings);
    // for (size_t i = 0; i < absCoeffs.Size(); i++)
    // {
    //     absCoeffs[i] = system.GetSteadyState(QSim::GetCol(detunings, i)).GetAbsCoeff("S1_2", "P3_2");
    // }
    

    auto absCoeffs = pool.Map([&](auto dets){ return system.GetSteadyState(dets).GetAbsCoeff("S1_2", "P3_2"); }, 
        QSim::GetColIteratorBegin(detunings), QSim::GetColIteratorEnd(detunings));
    
    pool.WaitUntilFinnished();
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;
    
    // Write to file
    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < detunings.Size(); i++)
        file << detunings[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
