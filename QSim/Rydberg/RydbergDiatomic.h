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

        // coefficients needed for the multielectron correction
        virtual double GetCoreInteractionME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const = 0;

        virtual double GetEnergy(const RydbergDiatomicState_t& state) const override;
        virtual double GetPotential(double r, const RydbergDiatomicState_t& state) const override;
        
        virtual double GetDipoleME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const override;

        double GetSelfDipoleME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const;
        
    protected:
        double GetHcbToHcdCoeff(int N, int l, int R, int lambda) const;
    };

    class NitricOxide : public RydbergDiatomic
    {
    private:
        // DOI: 10.1039/d1cp01930a
        constexpr static std::array<double, 1> s_l0quantumDefects = { 0.210 };
        constexpr static std::array<double, 2> s_l1quantumDefects = { 0.7038, 0.7410 };
        constexpr static std::array<double, 3> s_l2quantumDefects = { 0.050, -0.053, 0.089 };
        constexpr static std::array<double, 4> s_l3quantumDefects = { 0.0182, 0.0172, 0.00128, 0.0057 };
        constexpr static int s_quantumDefectsLmax = 3;
    public:
        NitricOxide();
        
        virtual double GetQuantumDefect(const RydbergDiatomicState_t& state) const override;
        virtual double GetRotationalConstant() const override;
        virtual double GetCentrifugalDistConstant() const override;
        virtual double GetCoreDipoleMoment() const override;

        virtual double GetCoreInteractionME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const override;
        double GetCoreInteractionHcb(int n1, int l1, int R1, int n2, int l2, int R2, int lambda, int N) const;
        
    };

}

#endif
