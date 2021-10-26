// Philipp Neufeld, 2021

#include <iostream>
#include <QSim/NLevelSystem.h>
#include <QSim/Matrix.h>
#include <QSim/TransitionTree.h>
#include <QSim/Doppler.h>
#include <QSim/ThreadPool.h>
#include <chrono>

#include <fstream>


std::mutex mutex;

void task()
{
    {
        std::unique_lock<std::mutex> lock(mutex);
        std::cout << "Starting task..." << std::endl;
    }

    std::this_thread::sleep_for(std::chrono::seconds(2));

    {
        std::unique_lock<std::mutex> lock(mutex);
        std::cout << "Finnishing task..." << std::endl;
    }
}

int main()
{
    
    QSim::ThreadPool pool;

    pool.AddTask(task);
    pool.AddTask(task);
    pool.AddTask(task);
    pool.WaitUntilFinnished();

    std::cout << "All tasks finnished." << std::endl;

    return 0;

    using Ty = double;

    Ty levels[] = { -4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9 };
    QSim::TNLevelSystem<3, Ty> system(levels);
    system.AddTransition(QSim::TTransition<Ty>{0, 2, 3.5e6});
    system.AddTransition(QSim::TTransition<Ty>{1, 2, 10e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 0, 3.0/8.0 * 6.065e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 1, 5.0/8.0 * 6.065e6});

    Ty mass = 1.44e-25;
    QSim::TDopplerIntegrator<Ty> dopplerIntegrator(mass, 300.0);
    
    Ty det_start = -100e6;
    Ty det_stop = 100e6;
    std::size_t det_steps = 501;
    std::vector<Ty> detuning;
    detuning.reserve(det_steps);
    for (std::size_t i = 0; i < det_steps; i++)
        detuning.push_back(det_start + i * (det_stop - det_start) / (det_steps - 1));

    auto start_ts = std::chrono::high_resolution_clock::now();

    std::vector<Ty> absCoeffs;
    absCoeffs.reserve(detuning.size());
    for (Ty det: detuning)
        absCoeffs.push_back(dopplerIntegrator.IntegrateAbsorptionCoefficient(system, QSim::TStaticMatrix<Ty, 2, 1>({ det, 0 }), 0, 2));
    
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < det_steps; i++)
        file << detuning[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
