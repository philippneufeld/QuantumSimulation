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
    double det_start = -100e6;
    double det_stop = 100e6;
    std::size_t det_steps = 501;
    std::vector<double> detuning;
    detuning.reserve(det_steps);
    for (std::size_t i = 0; i < det_steps; i++)
        detuning.push_back(det_start + i * (det_stop - det_start) / (det_steps - 1));

    std::map<std::string, double> levels;
    levels["g1"] = -4.271e9;
    levels["g2"] = 2.563e9;
    levels["e"] = QSim::SpeedOfLight_v / 780.241e-9;
    QSim::TStaticNLevelSystem<3> ssys(levels);
    ssys.AddTransition("g1", "e", 3.5e6);
    ssys.AddTransition("g2", "e", 10.0e6);
    ssys.AddDecay("e", "g1", 3.0/8.0 * 6.065e6);
    ssys.AddDecay("e", "g2", 5.0/8.0 * 6.065e6);

    QSim::TDopplerIntegrator<double> doppler(1.44316060e-25, 300.0);

    auto start_ts = std::chrono::high_resolution_clock::now();

    std::vector<double> absCoeffs;
    absCoeffs.resize(detuning.size());
    for (std::size_t i = 0; i < absCoeffs.size(); i++)
    {
        auto task = [&, i]()
        { 
            absCoeffs[i] = doppler.IntegrateAbsorptionCoefficient(
                ssys, QSim::TStaticMatrix<double, 2, 1>({ detuning[i], 0 }), "g1", "e"); 
        };
        pool.AddTask(task);
    }
    
    pool.WaitUntilFinnished();
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;
    
    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < det_steps; i++)
        file << detuning[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
