// Philipp Neufeld, 2021

#ifndef QSim_Doppler_H_
#define QSim_Doppler_H_

#include <cstdint>
#include <functional>

#include "../Math/Quadrature.h"
#include "../Constants.h"

namespace QSim
{

    template<typename Ty>
    class TDopplerIntegrator
    {
    public:
        // constructors
        TDopplerIntegrator() : TDopplerIntegrator(1.674e-27, 300) { }
        TDopplerIntegrator(Ty mass, Ty temperature) : TDopplerIntegrator(mass, temperature, 1000) { }
        TDopplerIntegrator(Ty mass, Ty temperature, std::size_t steps)
            : m_mass(mass), m_temperature(temperature), m_steps(steps) { }

        // copy operators
        TDopplerIntegrator(const TDopplerIntegrator&) = default;
        TDopplerIntegrator& operator=(const TDopplerIntegrator&) = default;

        // Parameter getter/setter
        void SetMass(double mass) { m_mass = mass; }
        double GetMass() const { return m_mass; }
        void SetTemperature(double temp) { m_temperature = temp; }
        double GetTemperature() const { return m_temperature; }
        void SetIntegrationSteps(std::size_t steps) { m_steps = steps; }
        std::size_t GetIntegrationSteps() const { return m_steps; }

        double GetDopplerWidth(double frequency) const;

        template<typename Lambda, typename Ret = decltype(std::declval<Lambda>()(std::declval<Ty>()))>
        Ret Integrate(Lambda func) const;
    
    private:
        Ty m_mass;
        Ty m_temperature;
        std::size_t m_steps;
    };

    template<typename Ty>
    template<typename Lambda, typename Ret>
    Ret TDopplerIntegrator<Ty>::Integrate(Lambda func) const
    {
        const static Ty pi = std::acos(-1.0);
        if (m_temperature > 0 && m_mass > 0 && m_steps > 1)
        {
            Ty sigma = std::sqrt(BoltzmannConstant_v * m_temperature / m_mass);
            Ty sigma2SqRec = 1 / (2 * sigma * sigma);
            Ty norm = 1 / (std::sqrt(2*pi)*sigma);
            
            QuadSimpsonAlt integrator;
            auto convFunc = [=](double v){ return norm * std::exp(-v*v*sigma2SqRec) * func(v); };
            return integrator.Integrate(convFunc, -3.5 * sigma, 3.5 * sigma, m_steps);
        }
        else
            return func(0);
    }
    
    template<typename Ty>
    double TDopplerIntegrator<Ty>::GetDopplerWidth(double frequency) const
    {
        constexpr double constants = 8*BoltzmannConstant_v*Ln2_v / (SpeedOfLight_v * SpeedOfLight_v); 
        return std::sqrt(constants * m_temperature / m_mass) * frequency;
    }

}

#endif
