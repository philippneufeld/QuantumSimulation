// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Wigner_H
#define QSim_Math_Wigner_H

#include <cmath>
#include "Gamma.h"

namespace QSim
{
    // Wigner angular momenta symbols
    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);
    double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6);

    // clebsh-gordan coefficients
    double ClebshGordan(double j1, double j2, double j3, 
        double m1, double m2, double m3);
}

#endif 
