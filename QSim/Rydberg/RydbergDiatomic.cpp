// Philipp Neufeld, 2021-2022

#include "../Constants.h"
#include "RydbergDiatomic.h"

namespace QSim
{

    RydbergDiatomic::RydbergDiatomic(double mass) 
        : TRydbergSystem<RydbergDiatomicState_t>(mass) {}
    
    double RydbergDiatomic::GetEnergy(const RydbergDiatomicState_t& state) const
    {
        double energy = 0.0;
        const auto [n, l, R, N, mN] = state;
        
        // rotational term
        constexpr double hc = PlanckConstant_v * SpeedOfLight_v;
        double rotTerm = R*(R+1);
        energy += hc * this->GetRotationalConstant() * rotTerm;
        energy += hc * this->GetCentrifugalDistConstant() * rotTerm*rotTerm;

        // rydberg term
        energy += this->GetRydbergEnergy(n, state);

        return energy;
    }

    double RydbergDiatomic::GetPotential(double r, const RydbergDiatomicState_t& state) const
    {
        double potential = 0.0;
        const auto [n, l, R, N, mN] = state;
        
        // rotational term
        constexpr double hc = PlanckConstant_v * SpeedOfLight_v;
        double rotTerm = R*(R+1);
        potential += hc * this->GetRotationalConstant() * rotTerm;
        potential += hc * this->GetCentrifugalDistConstant() * rotTerm*rotTerm;

        // atomic potential term
        potential += this->GetAtomicPotential(r, n, l);

        return potential;
    }

    double RydbergDiatomic::GetDipoleME(const RydbergDiatomicState_t& state1, const RydbergDiatomicState_t& state2) const
    {
        const auto [n1, l1, ml1, N1, mN1] = state1;
        const auto [n2, l2, ml2, N2, mN2] = state2;
        
        // rotational wavefunction integral
        if (N1 != N2 || mN1 != mN2)
            return 0.0;

        // TODO: angular part
        double dip = 0.0;
        
        // radial part
        if (dip != 0)
        {
            double rmax1 = 3*(n1+15)*n1*BohrRadius_v;
            double rmax2 = 3*(n2+15)*n2*BohrRadius_v;
            dip *= this->GetDipMeRadHelper(state1, state2, rmax1, rmax2, 50);
        }

        return dip * ElementaryCharge_v;
    }

    //
    // NitricOxide
    //

    NitricOxide::NitricOxide() 
        : RydbergDiatomic(30.006 * AtomicMassUnit_v) {}

    double NitricOxide::GetQuantumDefect(const RydbergDiatomicState_t& state) const
    {
        return 0.0;
    }

    double NitricOxide::GetRotationalConstant() const
    {
        // Phys. Chem. Chem. Phys., 2021, 23, 18806; doi: 10.1039/d1cp01930a
        return 198.7825; // m^-1
    }

    double NitricOxide::GetCentrifugalDistConstant() const
    {
        // Phys. Chem. Chem. Phys., 2021, 23, 18806; doi: 10.1039/d1cp01930a
        return 5.64e-4; // m^-1
    }



}
