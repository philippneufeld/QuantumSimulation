// Philipp Neufeld, 2021-2022

#ifndef QSim_QSim_Constants_H_
#define QSim_QSim_Constants_H_

#include "Platform.h"

namespace QSim
{
    // General constants
    constexpr static double Ln2_v = 0.69314718056;
    constexpr static double Pi_v = 3.14159265358979323846;
    constexpr static double TwoPi_v = 2 * Pi_v;

    // Physical constants
    constexpr static double SpeedOfLight_v = 2.99792458e8;
    constexpr static double BoltzmannConstant_v = 1.38064852e-23;
    constexpr static double VacuumPermittivity_v = 8.8541878128e-12;
    constexpr static double PlanckConstant_v = 6.62607004e-34;
    constexpr static double ReducedPlanckConstant_v = 1.054571817e-34;
    constexpr static double AtomicMassUnit_v = 1.66053906660e-27;
    constexpr static double ElementaryCharge_v = 1.602176462e-19;
    constexpr static double BohrRadius_v = 0.5291772083e-10;
    constexpr static double Debye_v = 0.39343 * ElementaryCharge_v * BohrRadius_v;

    // utility functions for constants
    constexpr double ConstexprSqrt(double a)
    {
        double x = a;
        for (std::size_t i = 0; i < 10; i++)
            x -= 0.5 * (x - a/x);
        return x;
    }
}

#endif
