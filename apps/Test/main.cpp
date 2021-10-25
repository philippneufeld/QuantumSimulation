// Philipp Neufeld, 2021

#include <iostream>
#include <QSim/NLevelSystem.h>
#include <QSim/Matrix.h>
#include <QSim/TransitionTree.h>
#include <chrono>

#include <fstream>

int main()
{
    
    using Ty = double;

    // Ty mass = 1.44e-25;
    Ty levels[] = { -4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9 };
    QSim::TNLevelSystem<3, Ty> system(levels);
    system.AddTransition(QSim::TTransition<Ty>{0, 2, 3.5e6});
    system.AddTransition(QSim::TTransition<Ty>{1, 2, 10e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 0, 3.0/8.0 * 6.065e6});
    system.AddDecay(QSim::TDecay<Ty>{2, 1, 5.0/8.0 * 6.065e6});
    
    Ty det_start = -0.025e9;
    Ty det_stop = 0.025e9;
    std::size_t det_steps = 501;
    std::vector<Ty> detuning;
    detuning.reserve(det_steps);
    for (std::size_t i = 0; i < det_steps; i++)
        detuning.push_back(det_start + i * (det_stop - det_start) / (det_steps - 1));

    auto start_ts = std::chrono::high_resolution_clock::now();

    std::vector<Ty> absCoeffs;
    absCoeffs.reserve(detuning.size());
    for (Ty det: detuning)
        absCoeffs.push_back(system.GetAbsorptionCoeff(QSim::TStaticMatrix<Ty, 2, 1>({ det, 0 }), 0, 2));
    
    std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

    std::ofstream file;
    file.open("data.txt", std::ios::out);
    for (std::size_t i = 0; i < det_steps; i++)
        file << detuning[i] << " " << absCoeffs[i] << std::endl;
    file.close();
    
    return 0;
}
