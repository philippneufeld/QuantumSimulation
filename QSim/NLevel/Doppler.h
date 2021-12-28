// Philipp Neufeld, 2021-2022

#ifndef QSim_Doppler_H_
#define QSim_Doppler_H_

#include <cstdint>
#include <functional>

#include "../Math/Quadrature.h"
#include "../Constants.h"

namespace QSim
{

    class DopplerIntegrator
    {
    public:
        // constructors
        DopplerIntegrator() : DopplerIntegrator(1.674e-27, 300) { }
        DopplerIntegrator(double mass, double temperature) : DopplerIntegrator(mass, temperature, 250) { }
        DopplerIntegrator(double mass, double temperature, std::size_t steps)
            : m_mass(mass), m_temperature(temperature), m_steps(steps), 
            m_rtol(1e-7), m_atol(1e-9), m_maxDepth(10), m_sigmas(3.5) { }

        // copy operators
        DopplerIntegrator(const DopplerIntegrator&) = default;
        DopplerIntegrator& operator=(const DopplerIntegrator&) = default;

        // Thermal parameters
        void SetMass(double mass) { m_mass = mass; }
        double GetMass() const { return m_mass; } 
        void SetTemperature(double temp) { m_temperature = temp; }
        double GetTemperature() const { return m_temperature; }
        
        // Integration parameters
        void SetIntegrationSteps(std::size_t steps) { m_steps = steps; }
        std::size_t GetIntegrationSteps() const { return m_steps; }
        void SetRTol(double rtol) { m_rtol = rtol; }
        double GetRTol() const { return m_rtol; }
        void SetATol(double atol) { m_atol = atol; }
        double GetATol() const { return m_atol; }
        void SetMaxRecursionDepth(double depth) { m_maxDepth = depth; }
        double GetMaxRecursionDepth() const { return m_maxDepth; }
        void SetIntegrationWidth(double sigmas) { m_sigmas = sigmas; }
        double GetIntegrationWidth(double sigmas) const { return m_sigmas; }

        double GetDopplerWidth(double frequency) const;

        template<typename Lambda, typename Ret = decltype(std::declval<Lambda>()(std::declval<double>()))>
        Ret Integrate(Lambda func) const;
    
    private:
        double m_mass;
        double m_temperature;
        double m_rtol;
        double m_atol;
        double m_maxDepth;
        double m_sigmas;
        std::size_t m_steps;
    };

    template<typename Lambda, typename Ret>
    Ret DopplerIntegrator::Integrate(Lambda func) const
    {
        const static double pi = std::acos(-1.0);
        if (m_temperature > 0 && m_mass > 0 && m_steps > 1)
        {
            double sigma = std::sqrt(BoltzmannConstant_v * m_temperature / m_mass);
            double sigma2SqRec = 1 / (2 * sigma * sigma);
            double norm = 1 / (std::sqrt(2*pi)*sigma);

            double vmin = -m_sigmas * sigma;
            double vmax = m_sigmas * sigma;

            TQuadAdaptive<double> integrator;
            auto convFunc = [=](double v){ return norm * std::exp(-v*v*sigma2SqRec) * func(v); };
            return integrator.Integrate(convFunc, vmin, vmax, m_steps, m_rtol, m_atol, m_maxDepth);
        }
        else
            return func(0);
    }
    
    double DopplerIntegrator::GetDopplerWidth(double frequency) const
    {
        constexpr double constants = 8*BoltzmannConstant_v*Ln2_v / (SpeedOfLight_v * SpeedOfLight_v); 
        return std::sqrt(constants * m_temperature / m_mass) * frequency;
    }

}

#endif
