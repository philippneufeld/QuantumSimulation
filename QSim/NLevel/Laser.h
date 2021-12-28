// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevel_Laser_H_
#define QSim_NLevel_Laser_H_

#include <cstdint>

#include "../Constants.h"

namespace QSim
{
    // General prupose conversion functions
    constexpr double GetIntensityFromPower(double power, double waistRadius);
    constexpr double GetElectricFieldFromIntensity(double intensity);
    constexpr double GetIntensityFromElectricField(double electricField);
    constexpr double GetRabiFrequencyFromIntensity(double dipole, double intensity);
    constexpr double GetIntensityFromRabiFrequency(double dipole, double rabi);

    namespace Internal
    {
        constexpr double ConstexprSqrt(double a)
        {
            double x = a;
            for (std::size_t i = 0; i < 10; i++)
                x -= 0.5 * (x - a/x);
            return x;
        }
    }

    constexpr double GetIntensityFromPower(double power, double waistRadius)
    {
        // I(r, 0) = I0 * exp(-2*r^2 / w0^2)
        // I0 = 2*P/(pi*w0^2)
        return 2 * power / (Pi_v * waistRadius * waistRadius);
    }

    constexpr double GetElectricFieldFromIntensity(double intensity)
    {
        constexpr double conv = 2 / (SpeedOfLight_v * VacuumPermittivity_v); 
        return Internal::ConstexprSqrt(conv * intensity);
    }

    constexpr double GetIntensityFromElectricField(double electricField)
    {
        constexpr double conv = 0.5 * SpeedOfLight_v * VacuumPermittivity_v; 
        return conv * electricField * electricField;
    }

    constexpr double GetRabiFrequencyFromIntensity(double dipole, double intensity)
    {
        return dipole * GetElectricFieldFromIntensity(intensity) / PlanckConstant_v;
    }

    constexpr double GetIntensityFromRabiFrequency(double dipole, double rabi)
    {
        return GetIntensityFromElectricField(PlanckConstant_v * rabi / dipole);
    }
}

#endif
