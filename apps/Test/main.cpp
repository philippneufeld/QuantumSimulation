// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Rydberg/RydbergDiatomic.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char *argv[])
{
    NitricOxide molecule;

    // H state -> 3d state
    // C. Jungen, Rydberg Series in the NO Spectrum: An Interpretation of Quantum Defects and Intensities in the sand d Series
    double hstate = molecule.GetEnergy(std::make_tuple(3, 2, 5, 0, 0));

    double redLaserWL = 835.0;

    double lambda = std::numeric_limits<double>::infinity();
    int n = 4;
    for (int i = n; i<50; i++)
    {
        double rydberg = molecule.GetEnergy(std::make_tuple(i, 3, 5, 0, 0));
        double lambda0 = SpeedOfLight_v * PlanckConstant_v / (rydberg - hstate) * 1e9;
        double err = std::abs(lambda-redLaserWL);
        
        if (std::abs(lambda0-redLaserWL) < std::abs(lambda-redLaserWL))
        {
            lambda = lambda0;
            n = i;
        }
    }

    std::cout << "Principal quantum number: n=" << n << " (" << lambda << "nm)" << std::endl;

    return 0;
}
