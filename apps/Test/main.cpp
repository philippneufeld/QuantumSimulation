// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Rydberg/RydbergDiatomic.h>

using namespace QSim;
using namespace Eigen;

int main()
{
    NitricOxide molecule;
    int n = 40;
    int l = 3;
    auto state = RydbergDiatomicState_t(n, l, 5, 5, 0);

    double hcR = PlanckConstant_v*SpeedOfLight_v*molecule.GetScaledRydbergConstant();

    double mu = molecule.GetQuantumDefect(state);
    double corr1 = -2*hcR*mu / std::pow(n, 3);
    double energy0 = -hcR / std::pow(n, 2);
    double energy = -hcR / std::pow(n-mu, 2);
    std::cout << "Hund's case (d): " 
        << mu << ", "
        << corr1 / EnergyGHz_v << " GHz, "
        << ((energy0+corr1)-energy) / EnergyGHz_v << " GHz, "
        << std::endl;

    return 0;
}
