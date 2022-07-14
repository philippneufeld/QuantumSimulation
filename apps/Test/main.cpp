// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/NLevel/Laser.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char *argv[])
{
    // add coupling lasers
    double uvInt = PowerToIntensity(0.05, 1e-3); // 50mW
    double greenInt = PowerToIntensity(1.0, 1e-3); // 1W
    double redInt = PowerToIntensity(0.5, 1e-3); // 500mW
    
    // print rabis
    std::cout << "UV Rabi:" << IntensityToRabi(0.1595 * Debye_v, uvInt) / 1e6 << "MHz" << std::endl;
    std::cout << "Green Rabi:" << IntensityToRabi(1e-2 * Debye_v, greenInt) / 1e6 << "MHz" << std::endl;
    std::cout << "Red Rabi:" << IntensityToRabi(1e-1 * Debye_v, redInt) / 1e6 << "MHz" << std::endl;

    return 0;
}
