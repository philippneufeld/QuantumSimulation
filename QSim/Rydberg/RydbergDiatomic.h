// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_RydbergDiatomic_H_
#define QSim_Rydberg_RydbergDiatomic_H_

// Includes
#include "RydbergSystem.h"

namespace QSim
{

    // n, l, ml, N+, mN+
    using RydbergDiatomicState_t = std::tuple<int, int, int, int, int>;
    
    class RydbergDiatomic : public TRydbergSystem<RydbergDiatomicState_t>
    {
    public:
        RydbergDiatomic(double mass);

        virtual double GetQuantumDefect(const RydbergDiatomicState_t& state) const override = 0;
        virtual double GetRotationalConstant() const = 0;

        virtual double GetEnergy(const RydbergDiatomicState_t& state) const override;
        virtual double GetPotential(double r, const RydbergDiatomicState_t& state) const override;
        virtual double GetDipoleME(const RydbergDiatomicState_t& state1, const RydbergDiatomicState_t& state2) const override;
    };

    class NitricOxide : public RydbergDiatomic
    {
    public:
        NitricOxide();
        
        virtual double GetScaledRydbergConstant() const override;
        virtual double GetQuantumDefect(const RydbergDiatomicState_t& state) const override;
        virtual double GetRotationalConstant() const override;
    };

}

#endif
