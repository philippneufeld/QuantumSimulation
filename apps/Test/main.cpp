// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Math/Wigner.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char *argv[])
{
    std::cout << Wigner6j(1, 1.5, 0.5, 1, 0.5, 1.5) / 0.263523 << std::endl;
    std::cout << Wigner6j(1, 1.5, 0.5, 1, 1.5, 0.5) / -0.0833333 << std::endl;
    std::cout << Wigner6j(1, 0, 1, 0, 1, 0) * sqrt(3) << std::endl;
    std::cout << Wigner6j(2, 1, 1, 2, 1, 1) * 30 << std::endl;
    std::cout << Wigner6j(3, 2, 1, 1, 2, 3) * 5*sqrt(7.0/2.0) << std::endl;

    std::cout << (Wigner6j(1, 2, 1, 0, 0, 0) == 0) << std::endl;
    std::cout << (Wigner6j(1, 1, 1, 0, 0, 0) == 0) << std::endl;

    std::cout << Wigner6j(9, 15, 21, 21, 21, 22) / 0.00674869 << std::endl;
    std::cout << (Wigner6j(9, 15, 21, 21, 21, 21) == 0) << std::endl;

    return 0;
}
