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
    auto probeDetunings = QSim::CreateLinspace(-100.0e6, 100.0e6, 5001);

    std::map<std::string, double> levels;
    levels["g1"] = -4.271e9;
    levels["g2"] = 2.563e9;
    levels["e"] = QSim::SpeedOfLight_v / 780.241e-9;
    QSim::TStaticNLevelSystem<3> sys(levels);
    sys.AddTransition("g1", "e", 3.5e6);
    sys.AddTransition("g2", "e", 10.0e6);
    sys.AddDecay("e", "g1", 3.0/8.0 * 6.065e6);
    sys.AddDecay("e", "g2", 5.0/8.0 * 6.065e6);

    QSim::TDopplerIntegrator<double> doppler(1.44316060e-25, 300.0);

    auto start_ts = std::chrono::high_resolution_clock::now();

    std::vector<double> absCoeffs;
    absCoeffs.resize(probeDetunings.Size());
    for (std::size_t i = 0; i < absCoeffs.size(); i++)
    {
        auto task = [&, i]()
        { 
            QSim::TStaticColVector<double, 2> detunings({ probeDetunings[i], 0.0 });
            absCoeffs[i] = doppler.IntegrateAbsorptionCoefficient(sys, detunings, "g1", "e"); 
        };
        pool.AddTask(task);
    }
    
    pool.WaitUntilFinnished();
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;
    
    // Write to file
    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < probeDetunings.Size(); i++)
        file << probeDetunings[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
