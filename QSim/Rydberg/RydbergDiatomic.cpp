// Philipp Neufeld, 2021-2022

#include "RydbergDiatomic.h"

#include <iostream>
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
        energy -= hc * this->GetCentrifugalDistConstant() * rotTerm*rotTerm;

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
        potential -= hc * this->GetCentrifugalDistConstant() * rotTerm*rotTerm;

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
        double mu2 = this->GetQuantumDefect(state2);
        double nu1 = n1 - mu1;
        double nu2 = n2 - mu2;
        double lambda1 = l1 - mu1;
        double lambda2 = l2 - mu2;

        // approximation: see Phys. Chem. Chem. Phys.,2021, 23, 18806
        double rad = 2.0 / (a0*a0 * std::pow(nu1*nu2, 1.5) * (lambda1 + lambda2 + 1));
        rad *= std::sin(Pi_v * (lambda1 - lambda2)) / (Pi_v*(lambda1 - lambda2));

        return k * mu * rad * f;
    }

    double RydbergDiatomic::GetSelfMultiElectronME(const RydbergDiatomicState_t& state1, 
        const RydbergDiatomicState_t& state2) const
    {
        const auto [n1, l1, R1, N1, mN1] = state1;
        const auto [n2, l2, R2, N2, mN2] = state2;

        // selection rules
        if (N1 != N2 || mN1 != mN2 || std::abs(R1-R2) != 2) // TODO: Delta R selection rule correct here?
            return 0.0;

        double hcoeff = this->GetConfigurationMixingCoeff(l1, R1, l2, R2, N1);
        if (hcoeff == 0.0)
            return 0.0;
            
        // adjusted quantum number n
        double nu1 = n1 - this->GetQuantumDefect(state1);
        double nu2 = n2 - this->GetQuantumDefect(state2);

        constexpr double k = -2 * PlanckConstant_v * SpeedOfLight_v;
        return k * this->GetScaledRydbergConstant() / std::pow(nu1*nu2, 1.5) * hcoeff;
    }

    double RydbergDiatomic::GetHcbToHcdCoeff(int N, int l, int R, int lambda) const
    {
        double kronL0 = (lambda == 0 ? 1 : 0);
        double coeff = std::sqrt((2*R + 1) * 2.0 / (1 + kronL0));
        coeff *= Wigner3j(l, N, R, -lambda, lambda, 0);
        coeff *= (l+lambda-R) % 2 == 0 ? 1 : -1;
        return coeff;
    }

    //
    // NitricOxide
    //

    NitricOxide::NitricOxide() 
        : RydbergDiatomic(30.006 * AtomicMassUnit_v) {}

    double NitricOxide::GetQuantumDefect(const RydbergDiatomicState_t& state) const
    {
        const auto [n, l, R, N, mN] = state;

        const static std::array<const double*, 4> quantumDefects = {
            s_l0quantumDefects.data(), s_l1quantumDefects.data(),
            s_l2quantumDefects.data(), s_l3quantumDefects.data()
        };

        constexpr int lmax = quantumDefects.size() - 1;
        int l0 = std::min(l, lmax);

        double defect = 0;
        for (int Lambda=0; Lambda <= l0; Lambda++)
        {
            double aCoeff = GetHcbToHcdCoeff(N, l0, R, Lambda);
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

    double NitricOxide::GetConfigurationMixingCoeff(int l1, int R1, int l2, int R2, int N) const
    {
        const static std::array<const double*, 4> quantumDefects = {
            s_l0quantumDefects.data(), s_l1quantumDefects.data(),
            s_l2quantumDefects.data(), s_l3quantumDefects.data()
        };

        constexpr double mixingAngle = -38.7 * Pi_v / 180.0;
        double c = std::cos(mixingAngle);
        double s = std::sin(mixingAngle);
        double s2 = 0.5 * std::sin(2*mixingAngle);

        double result = 0.0;
        for (int lambda = 0; lambda < quantumDefects.size(); lambda++)
        {
            double hinner = 0.0;
            if (lambda == 0 && ((l1 == 0 && l2 == 2) || (l1 == 2 && l2 == 0)))
                    hinner = s2 * (quantumDefects[0][0] + quantumDefects[2][0]);
            else if (lambda == 0 && l1 == 0 && l2 == 0)
                    hinner = c*c*quantumDefects[0][0] + s*s*quantumDefects[2][0];
            else if (lambda == 0 && l1 == 2 && l2 == 2)
                    hinner = s*s*quantumDefects[0][0] + c*c*quantumDefects[2][0];
            else if (l1 == l2 && l1 <= lambda)
                hinner = quantumDefects[l1][lambda];
            else
                continue;
            
            result += hinner * 
                this->GetHcbToHcdCoeff(N, l1, R1, lambda) * 
                this->GetHcbToHcdCoeff(N, l2, R2, lambda);
        }

        return result;
    }


}
