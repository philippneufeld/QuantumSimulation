// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevel_Laser_H_
#define QSim_NLevel_Laser_H_

#include <cstdint>
#include <cmath>
#include <functional>

#include "../Constants.h"

namespace QSim
{

    double IntensityToEField(double intensity)
    {
        constexpr double conv = 2 / (SpeedOfLight_v * VacuumPermittivity_v); 
        return std::sqrt(conv * intensity);
    }
    
    double EFieldToIntensity(double eField)
    {
        constexpr double conv = 0.5 * SpeedOfLight_v * VacuumPermittivity_v; 
        return conv * eField * eField;
    }

    double PowerToIntensity(double power, double waistRadius)
    {
        // I(r, 0) = I0 * exp(-2*r^2 / w0^2)
        // I0 = 2*P/(pi*w0^2)
        return 2 * power / (Pi_v * waistRadius * waistRadius);
    }

    double RabiToIntensity(double dipole, double rabi)
    {
        return EFieldToIntensity(PlanckConstant_v * rabi / dipole);
    }
    
    double IntensityToRabi(double dipole, double intensity)
    {
        return dipole * IntensityToEField(intensity) / PlanckConstant_v;
    }

}

#endif
