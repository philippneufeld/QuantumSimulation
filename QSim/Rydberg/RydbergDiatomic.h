// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_RydbergDiatomic_H_
#define QSim_Rydberg_RydbergDiatomic_H_

// Includes
#include "RydbergSystem.h"

namespace QSim
{

    // Hund's case (d) basis
    // n, l, R, N, mN
    // n: Principle quantum number of the Rydberg electron
    // l: Orbital momentum quantum number of the Rydberg electron
    // R: Rotational core angular momentum
    // N: Total angular momentum (excluding spin)
    // mN: Projection of N onto the z-axis
    using RydbergDiatomicState_t = std::tuple<int, int, int, int, int>;
    
    class RydbergDiatomic : public TRydbergSystem<RydbergDiatomicState_t>
    {
    public:
        RydbergDiatomic(double mass);

        // Retrieve molecule specific constants
        virtual double GetRotationalConstant() const = 0;
        virtual double GetCentrifugalDistConstant() const = 0;
        virtual double GetCoreDipoleMoment() const = 0;

        virtual double GetEnergy(const RydbergDiatomicState_t& state) const override;
        virtual double GetPotential(double r, const RydbergDiatomicState_t& state) const override;
        
        virtual double GetDipoleME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const override;

        double GetSelfDipoleME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const;
    };

    class NitricOxide : public RydbergDiatomic
    {
    public:
        NitricOxide();
        
        virtual double GetQuantumDefect(const RydbergDiatomicState_t& state) const override;
        virtual double GetRotationalConstant() const override;
        virtual double GetCentrifugalDistConstant() const override;
        virtual double GetCoreDipoleMoment() const override;
    };

}

#endif
