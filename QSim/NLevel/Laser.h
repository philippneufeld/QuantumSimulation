// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevel_Laser_H_
#define QSim_NLevel_Laser_H_

#include <cstdint>
#include <functional>

#include "../Constants.h"

namespace QSim
{
    
    // Helper that handles the laser modulation
    namespace Internal
    {
        template<bool AM>
        class TNLevelLaserModulator
        {
        protected:
            ~TNLevelLaserModulator() = default;
        public:
            TNLevelLaserModulator() 
                : m_modulation([](double){ return 1.0; }) {}

            template<typename Func, typename=std::enable_if_t<std::is_invocable_r_v<double, Func, double>>>
            void SetModulationFunc(Func&& func) { m_modulation = std::forward<Func>(func); };
            double GetModulation(double t) const { return std::invoke(m_modulation, t); }
        private:
            std::function<double(double)> m_modulation;
        };

        template<>
        class TNLevelLaserModulator<false>
        {
        protected:
            ~TNLevelLaserModulator() = default;
        public:
            TNLevelLaserModulator() {}
            double GetModulation(double t) const { return 1.0; }
        };
    }

    template<bool AM=false> // allow amplitude modulation
    class TNLevelLaser : public Internal::TNLevelLaserModulator<AM>
    {
    public:
        TNLevelLaser(std::pair<unsigned int, unsigned int> lvls, 
            double intensity=0.0, double prop=1.0);

        // property setter
        void SetElAmplitude(double el) { m_elAmplitude = el; }
        void SetIntensity(double intensity);
        void SetPropDirection(double dir) { m_propDir = std::max(-1.0, std::min(1.0, dir)); }

        // property getter
        auto GetLevels() const { return m_lvls; }
        double GetElAmplitude() const { return m_elAmplitude; }
        double GetIntensity() const;
        double GetPropDirection() const { return m_propDir; }

        double GetModElAmplitude(double t) const { return m_elAmplitude * this->GetModulation(t); }
        double GetModIntensity(double t) const;

        // static conversion functions
        static double IntensityToEField(double intensity);
        static double EFieldToIntensity(double eField);
        static double PowerToIntensity(double power, double waistRadius);
        static double RabiToIntensity(double dipole, double rabi);
        static double IntensityToRabi(double dipole, double intensity);
        
    private:
        const std::pair<unsigned int, unsigned int> m_lvls;
        double m_elAmplitude;
        double m_propDir;
    };

    template<bool AM>
    TNLevelLaser<AM>::TNLevelLaser(std::pair<unsigned int, unsigned int> lvls, 
        double intensity, double prop)
        : m_lvls(lvls), m_elAmplitude(0), m_propDir(prop) 
    {
        SetIntensity(intensity);
    }

    template<bool AM>
    void TNLevelLaser<AM>::SetIntensity(double intensity)
    {
        m_elAmplitude = IntensityToEField(intensity);
    }

    template<bool AM>
    double TNLevelLaser<AM>::GetIntensity() const 
    {
        return EFieldToIntensity(m_elAmplitude);
    }

    template<bool AM>
    double TNLevelLaser<AM>::GetModIntensity(double t) const
    {
        double el = GetModElAmplitude(t);
        return EFieldToIntensity(el);
    }

    template<bool AM>
    double TNLevelLaser<AM>::IntensityToEField(double intensity)
    {
        constexpr double conv = 2 / (SpeedOfLight_v * VacuumPermittivity_v); 
        return std::sqrt(conv * intensity);
    }
    
    template<bool AM>
    double TNLevelLaser<AM>::EFieldToIntensity(double eField)
    {
        constexpr double conv = 0.5 * SpeedOfLight_v * VacuumPermittivity_v; 
        return conv * eField * eField;
    }

    template<bool AM>
    double TNLevelLaser<AM>::PowerToIntensity(double power, double waistRadius)
    {
        // I(r, 0) = I0 * exp(-2*r^2 / w0^2)
        // I0 = 2*P/(pi*w0^2)
        return 2 * power / (Pi_v * waistRadius * waistRadius);
    }

    template<bool AM>
    double TNLevelLaser<AM>::RabiToIntensity(double dipole, double rabi)
    {
        return EFieldToIntensity(PlanckConstant_v * rabi / dipole);
    }
    
    template<bool AM>
    double TNLevelLaser<AM>::IntensityToRabi(double dipole, double intensity)
    {
        return dipole * IntensityToEField(intensity) / PlanckConstant_v;
    }


    // typedefs
    using ModulatedNLevelLaser = TNLevelLaser<true>;
    using NLevelLaser = TNLevelLaser<false>;
    
}

#endif
