// Philipp Neufeld, 2021-2022

#ifndef QSim_Doppler_H_
#define QSim_Doppler_H_

#include <cstdint>
#include <functional>

#include "../Math/MatrixTraits.h"
#include "../Math/Quadrature.h"
#include "../Constants.h"

namespace QSim
{

    namespace Internal
    {
        template<typename QP>
        class TDopplerIntegratorHelper
        {
        protected:
            ~TDopplerIntegratorHelper() { }
            
        public:
            TDopplerIntegratorHelper() { }
        };

        template<>
        class TDopplerIntegratorHelper<QuadAdaptivePolicy> : public QuadAdaptivePolicy
        {
        protected:
            ~TDopplerIntegratorHelper() { }
            
        public:
            TDopplerIntegratorHelper() 
            {
                // maximal machine precision is not required here
                this->SetIntegrationRTol(1e-6);
                this->SetIntegrationDepth(10);
            }
        };

    }

    template<typename QuadPolicy=QuadAdaptivePolicy>
    class TDopplerIntegrator : public Internal::TDopplerIntegratorHelper<QuadAdaptivePolicy>
    {
    public:
        // constructors
        TDopplerIntegrator() : TDopplerIntegrator(1.674e-27, 300) { }
        TDopplerIntegrator(double mass, double temperature) : TDopplerIntegrator(mass, temperature, 250) { }
        TDopplerIntegrator(double mass, double temperature, std::size_t steps)
            : m_mass(mass), m_temperature(temperature), m_steps(steps), m_sigmas(3.5) { }

        // copy operators
        TDopplerIntegrator(const TDopplerIntegrator&) = default;
        TDopplerIntegrator& operator=(const TDopplerIntegrator&) = default;

        // Thermal parameters
        void SetMass(double mass) { m_mass = mass; }
        double GetMass() const { return m_mass; } 
        void SetTemperature(double temp) { m_temperature = temp; }
        double GetTemperature() const { return m_temperature; }
        
        // Integration parameters
        void SetIntegrationSteps(std::size_t steps) { m_steps = steps; }
        std::size_t GetIntegrationSteps() const { return m_steps; }
        void SetIntegrationWidth(double sigmas) { m_sigmas = sigmas; }
        double GetIntegrationWidth(double sigmas) const { return m_sigmas; }

        double GetDopplerWidth(double frequency) const;

        template<typename Lambda, typename Ret = decltype(std::declval<Lambda>()(std::declval<double>()))>
        Ret Integrate(Lambda func) const;
    
    private:
        double m_mass;
        double m_temperature;
        double m_sigmas;
        std::size_t m_steps;
    };

    template<typename QuadPolicy>
    template<typename Lambda, typename Ret>
    Ret TDopplerIntegrator<QuadPolicy>::Integrate(Lambda func) const
    {
        const static double pi = std::acos(-1.0);
        if (m_temperature > 0 && m_mass > 0 && m_steps > 1)
        {
            double sigma = std::sqrt(BoltzmannConstant_v * m_temperature / m_mass);
            double sigma2SqRec = 1 / (2 * sigma * sigma);
            double norm = 1 / (std::sqrt(2*pi)*sigma);

            double vmin = -m_sigmas * sigma;
            double vmax = m_sigmas * sigma;

            auto convFunc = [=](double v) 
            { 
                auto res =  norm * std::exp(-v*v*sigma2SqRec) * func(v);
                return static_cast<TMatrixEvalType_t<decltype(res)>>(res);
            };
            return QuadPolicy::Integrate(convFunc, vmin, vmax, m_steps);
        }
        else
            return func(0);
    }
    
    template<typename QuadPolicy>
    double TDopplerIntegrator<QuadPolicy>::GetDopplerWidth(double frequency) const
    {
        constexpr double constants = 8*BoltzmannConstant_v*Ln2_v / (SpeedOfLight_v * SpeedOfLight_v); 
        return std::sqrt(constants * m_temperature / m_mass) * frequency;
    }

}

#endif
