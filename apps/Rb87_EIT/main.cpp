// Philipp Neufeld, 2021

#include <iostream>
#include <chrono>
#include <fstream>

#include <QSim/Math/Matrix.h>
#include <QSim/StaticQSysSS.h>
#include <QSim/Doppler.h>
#include <QSim/Util/ThreadPool.h>

int main()
{
    QSim::ThreadPool pool;

    // Generate detuning axis
    constexpr static std::size_t cnt = 501;
    QSim::TStaticMatrix<double, 2, cnt> detunings;
    auto probeDetunings = QSim::CreateLinspaceRow(-100.0e6, 100.0e6, cnt);
    QSim::SetRow(detunings, probeDetunings, 0);

    // setup Rb87 parameters
    std::map<std::string, double> levels;
    levels["S1_2_F1"] = -4.271e9;
    levels["S1_2_F2"] = 2.563e9;
    levels["P3_2"] = QSim::SpeedOfLight_v / 780.241e-9;
    double mass = 1.44316060e-25;
    double temperature = 300.0;

    // create system object
    QSim::TStaticQSysSS<3> system(levels, mass);
    system.AddTransition("S1_2_F1", "P3_2", 3.5e6);
    system.AddTransition("S1_2_F2", "P3_2", 10.0e6);
    system.AddDecay("P3_2", "S1_2_F1", 3.0/8.0 * 6.065e6);
    system.AddDecay("P3_2", "S1_2_F2", 5.0/8.0 * 6.065e6);
    system.SetTemperature(temperature);

    auto start_ts = std::chrono::high_resolution_clock::now();

    auto absCoeffs = pool.Map([&](auto dets){ return system.GetSteadyState(dets).GetAbsCoeff("S1_2_F1", "P3_2"); }, 
        QSim::GetColIteratorBegin(detunings), QSim::GetColIteratorEnd(detunings));

    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;
    
    // Write to file
    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < probeDetunings.Size(); i++)
        file << probeDetunings[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}

