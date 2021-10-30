// Philipp Neufeld, 2021

#include <iostream>
#include <chrono>
#include <fstream>

#include <QSim/NLevelSystem.h>
#include <QSim/Matrix.h>
#include <QSim/TransitionTree.h>
#include <QSim/Doppler.h>
#include <QSim/ThreadPool.h>
#include <QSim/NLevelSystem2.h>

#define ALT_IMPLD

int main()
{
    QSim::ThreadPool pool;

    double det_start = -100e6;
    double det_stop = 100e6;
    std::size_t det_steps = 501*1000;
    std::vector<double> detuning;
    detuning.reserve(det_steps);
    for (std::size_t i = 0; i < det_steps; i++)
        detuning.push_back(det_start + i * (det_stop - det_start) / (det_steps - 1));

#ifdef ALT_IMPLD
    QSim::NLevelSystem sys;
    sys.SetLevel("g1", -4.271e9);
    sys.SetLevel("g2", 2.563e9);
    sys.SetLevel("e", QSim::SpeedOfLight_v / 780.241e-9);

    std::cout << sys.GetLevel("e") << std::endl;

    bool success = true;
    success = sys.AddTransition("g1", "e", 1.0);
    success = sys.AddTransition("g2", "e", 1.0);
    success = sys.AddDecay("e", "g1", 1.0);
    success = sys.AddDecay("e", "g2", 1.0);
    std::cout << success << std::endl;

    auto start_ts = std::chrono::high_resolution_clock::now();
    std::vector<double> absCoeffs;
    absCoeffs.resize(detuning.size());
    for (std::size_t i = 0; i < absCoeffs.size(); i++)
    {
        auto task = [&, i]()
        { 
            absCoeffs[i] = sys.GetAbsorptionCoeff(
                QSim::TDynamicMatrix<double>(2, 1, { detuning[i], 0 }), 0.0, "g1", "e"); 
        };
        pool.AddTask(task);
    }
#else
    using Ty = double;
    Ty levels[] = { -4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9 };
    QSim::TNLevelSystem<3, Ty> system(levels);
    system.AddTransition(QSim::TTransition<Ty>{0, 2, 3.5e6});
    system.AddTransition(QSim::TTransition<Ty>{1, 2, 10e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 0, 3.0/8.0 * 6.065e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 1, 5.0/8.0 * 6.065e6});
    
    Ty mass = 1.44e-25;
    QSim::TDopplerIntegrator<Ty> dopplerIntegrator(mass, 300.0);
    
    auto start_ts = std::chrono::high_resolution_clock::now();

    std::vector<Ty> absCoeffs;
    absCoeffs.resize(detuning.size());
    for (std::size_t i = 0; i < absCoeffs.size(); i++)
    {
        auto task = [&, i]()
        { 
            absCoeffs[i] = dopplerIntegrator.IntegrateAbsorptionCoefficient(
                system, QSim::TStaticMatrix<Ty, 2, 1>({ detuning[i], 0 }), 0, 2); 
        };
        pool.AddTask(task);
    }
#endif
    
    pool.WaitUntilFinnished();
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < det_steps; i++)
        file << detuning[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
