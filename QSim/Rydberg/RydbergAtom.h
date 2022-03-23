// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_RydbergAtom_H_
#define QSim_Rydberg_RydbergAtom_H_

#include "RydbergSystem.h"

namespace QSim
{

    // n, l, j, mj
    using RydbergAtomState_t = std::tuple<int, int, double, double>;
    
    class RydbergAtom : public TRydbergSystem<RydbergAtomState_t>
    {
    public:
        RydbergAtom(double mass);

        virtual double GetQuantumDefect(const RydbergAtomState_t& state) const override = 0;

        virtual double GetEnergy(const RydbergAtomState_t& state) const override;
        virtual double GetPotential(double r, const RydbergAtomState_t& state) const override;
        virtual double GetDipoleME(const RydbergAtomState_t& state1, const RydbergAtomState_t& state2) const override;
    };

    class Hydrogen : public RydbergAtom
    {
    public:
        Hydrogen();
        virtual double GetQuantumDefect(const RydbergAtomState_t& state) const override;
    };

    class Rubidium : public RydbergAtom
    {
    public:
        Rubidium();
        virtual double GetQuantumDefect(const RydbergAtomState_t& state) const override;
    };


}

#endif
