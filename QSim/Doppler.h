// Philipp Neufeld, 2021

#ifndef QSim_Doppler_H_
#define QSim_Doppler_H_

#include <cstdint>
#include <functional>

namespace QSim
{

    constexpr static double BoltzmannConstant_v = 1.38064852e-23;

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
            int steps = static_cast<int>(m_steps / 2);
            Ty sigma = std::sqrt(BoltzmannConstant_v * m_temperature / m_mass);
            Ty sigma2SqRec = 1 / (2 * sigma * sigma);
            Ty norm = 1 / (std::sqrt(2*pi)*sigma);
            Ty vel_step = 3.5 * sigma / steps;
  
            Ret integrated;
            for (auto i = -steps; i <= steps; i++)
            {
                Ty velocity = i * vel_step;
                Ty thermal = norm * std::exp(-velocity*velocity * sigma2SqRec);
                Ret absCoeff = func(velocity);

                if (i == -steps) 
                    integrated = thermal * absCoeff;
                else
                    integrated += thermal * absCoeff;
            }
            integrated *= vel_step;
            return integrated;
        }
        else
            return func(0);
    }
}

#endif
