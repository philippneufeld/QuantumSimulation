// Philipp Neufeld, 2021-2022

#include "RydbergDiatomic.h"

#include "../Constants.h"
#include "QuantumDefects.h"


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
        const auto [n1, l1, R1, N1, mN1] = state1;
        const auto [n2, l2, R2, N2, mN2] = state2;
        
        // selection rules
        if (std::abs(l1-l2) != 1 || R1 != R2)
            return 0.0;

        // angular part
        double dip = std::sqrt((2*N1+1)*(2*N2+1));
        dip *= Wigner3j(N1, 1, N2, -mN1, 0, mN2);
        dip *= Wigner6j(l1, N1, R1, N2, l2, 1);
        dip *= (N1-mN1+N2+R1+l1+1) % 2 == 0 ? 1.0 : -1.0;

        // check if angular momentum algebra already results 
        // in a vanishing matrix element
        if (dip == 0) 
            return 0.0;

        // radial rydberg matrix element
        double rmax1 = 3*(n1+15)*n1*BohrRadius_v;
        double rmax2 = 3*(n2+15)*n2*BohrRadius_v;
        dip *= this->GetDipMeRadHelper(state1, state2, rmax1, rmax2, 50);
        dip *= std::sqrt(std::max(l1, l2)) * (l2 > l1 ? -1 : 1);

        return dip * ElementaryCharge_v;
    }

    double RydbergDiatomic::GetSelfDipoleME(const RydbergDiatomicState_t& state1, 
        const RydbergDiatomicState_t& state2) const
    {
        const auto [n1, l1, R1, N1, mN1] = state1;
        const auto [n2, l2, R2, N2, mN2] = state2;

        // selection rules
        if (N1 != N2 || mN1 != mN2 || std::abs(l2-l1) != 1 || std::abs(R2-R1) != 1)
            return 0.0;

        double f = ((l1 + l2 + N2) % 2 == 0 ? 1.0 : -1.0);
        f *= Wigner6j(N2, R1, l1, 1, l2, R2);
        f *= std::sqrt((2*R1+1)*(2*R2+1)) * Wigner3j(R1, 1, R2, 0, 0, 0);
        f *= std::sqrt((2*l1+1)*(2*l2+1)) * Wigner3j(l1, 1, l2, 0, 0, 0);

        constexpr double k = -ElementaryCharge_v / (4*Pi_v*VacuumPermittivity_v);
        double mu = this->GetCoreDipoleMoment();
        double a0 = BohrRadius_v * ElectronMass_v / this->GetReducedMass();
            
        // adjusted quantum numbers l and n
        double mu1 = this->GetQuantumDefect(state1);
        double mu2 = this->GetQuantumDefect(state1);
        double nu1 = n1 - mu1;
        double nu2 = n1 - mu2;
        double lambda1 = l1 - mu1;
        double lambda2 = l1 - mu2;

        // approximation: see Phys. Chem. Chem. Phys.,2021, 23, 18806
        double rad = 2 / (a0*a0 * std::pow(nu1*nu2, 1.5) * (lambda1 + lambda2 + 1));
        rad *= std::sin(Pi_v * (lambda1 - lambda2)) / (Pi_v*(lambda1 - lambda2));

        return k * mu * rad * f;
    }

    //
    // NitricOxide
    //

    NitricOxide::NitricOxide() 
        : RydbergDiatomic(30.006 * AtomicMassUnit_v) {}

    double NitricOxide::GetQuantumDefect(const RydbergDiatomicState_t& state) const
    {
        const auto [n, l, R, N, mN] = state;

        // DOI: 10.1039/d1cp01930a
        std::array<double, 1> l0quantumDefects = { 0.210 };
        std::array<double, 2> l1quantumDefects = { 0.7038, 0.7410 };
        std::array<double, 3> l2quantumDefects = { 0.050, -0.053, 0.089 };
        std::array<double, 4> l3quantumDefects = { 0.0182, 0.0172, 0.00128, 0.0057 };

        std::array<double*, 4> quantumDefects = {
            l0quantumDefects.data(), l1quantumDefects.data(),
            l2quantumDefects.data(), l3quantumDefects.data()
        };

        constexpr int lmax = quantumDefects.size() - 1;
        int l0 = std::min(l, lmax);

        double defect = 0;
        for (int Lambda=0; Lambda <= l0; Lambda++)
        {
            // projection coefficient
            double kronL0 = (Lambda == 0 ? 1 : 0);
            double aCoeff = std::sqrt((2*R + 1) * 2.0 / (1 + kronL0));
            aCoeff *= Wigner3j(l0, N, R, -Lambda, Lambda, 0);
            aCoeff *= (l0+Lambda-R) % 2 == 0 ? 1 : -1;

            defect += aCoeff*aCoeff*quantumDefects[l0][Lambda];
        }

        if (l0 < l)
            defect = 0.0; // ExtrapolateQuantumDefect(defect, l0, l);

        return defect;
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

    double NitricOxide::GetCoreDipoleMoment() const
    {
        // Phys. Chem. Chem. Phys., 2021, 23, 18806; doi: 10.1039/d1cp01930a
        return 0.4 * Debye_v;
    }



}
