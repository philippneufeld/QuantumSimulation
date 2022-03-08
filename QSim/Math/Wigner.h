// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Wigner_H
#define QSim_Math_Wigner_H

#include <cmath>
#include "Gamma.h"

namespace QSim
{
    double Wigner3j(
        double j1, double j2, double j3, 
        double m1, double m2, double m3);

    double ClebshGordan(
        double j1, double j2, double j3, 
        double m1, double m2, double m3);
}

#endif 
