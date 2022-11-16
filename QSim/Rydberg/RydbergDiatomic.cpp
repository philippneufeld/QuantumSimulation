// Philipp Neufeld, 2021-2022

#include "RydbergDiatomic.h"

#include <iostream>
#include "../Constants.h"
#include "QuantumDefects.h"

#define DEBUG_PRINTING
#ifdef DEBUG_PRINTING
#include <iostream>
#endif

#define QSIM_STARK_VRAKKING
// #define QSIM_STARK_HOGAN
// #define QSIM_STARK_QNOSE

namespace QSim
{

    RydbergDiatomic::RydbergDiatomic(double mass) 
        : TRydbergSystem<RydbergDiatomicState_t>(mass) {}
    
    int RydbergDiatomic::GetPrincipalQN(const RydbergDiatomicState_t& state) const
    {
        const auto [n, l, R, N, mN] = state;
        return n;
    }

    double RydbergDiatomic::GetQuantumDefectHcbMixing(int l1, int l2, int Lambda) const
    {
        // s-d mixing in the Lambda=0 channel
        if ((l1==0||l1==2) && (l2==0||l2==2) && Lambda == 0)
        {
            double sdAngle = this->GetSDMixingAngle();
            double c2 = [](double x){ return x*x; }(std::cos(sdAngle));
            double s2 = 1 - c2;

            double mus = this->GetQuantumDefectHcb(0, 0);
            double mud = this->GetQuantumDefectHcb(2, 0);

            if (l1 == 0 && l2 == 0)
                return (c2*mus + s2*mud);
            else if (l1 == 2 && l2 == 2)
                return (c2*mud + s2*mus);
            else
                return 0.5*std::sin(2*sdAngle) * (mus - mud);

        }
        else if (l1 == l2)
        {
            return this->GetQuantumDefectHcb(l1, Lambda);
        }

#ifdef DEBUG_PRINTING
        std::cout << "Something bad happened" << l1 << " " << l2 << std::endl;
        abort();
#endif
    }

    double RydbergDiatomic::GetQuantumDefect(const RydbergDiatomicState_t& state) const
    {
        // Calculates the Hund's case (d) quantum defect
        const auto [n, l, R, N, mN] = state;

        if (l > this->GetQuantumDefectHcbMaxLambda())
            return 0.0;

        double defect = 0;
        for (int Lambda=0; Lambda <= l; Lambda++)
        {
            double aCoeff = GetHcbToHcdCoeff(N, l, R, Lambda);
            defect += aCoeff*aCoeff*this->GetQuantumDefectHcbMixing(l, l, Lambda);
        }

        return defect;
    }

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
        double dlambda = lambda1 - lambda2;

        // approximation: see Phys. Chem. Chem. Phys.,2021, 23, 18806
        double rad = 2.0 / (a0*a0 * std::pow(nu1*nu2, 1.5) * (lambda1 + lambda2 + 1));
        rad *= std::sin(Pi_v*dlambda) / (Pi_v*dlambda);

        return k * mu * rad * f;
    }

    double RydbergDiatomic::GetCoreInteractionHcb(int n1, int l1, int R1, 
        int n2, int l2, int R2, int lambda, int N) const
    {
        // Hogen paper: H_ll^Lambda
        return this->GetQuantumDefectHcbMixing(l1, l2, lambda);
    }

    double RydbergDiatomic::GetCoreInteractionME(const RydbergDiatomicState_t& state1, 
            const RydbergDiatomicState_t& state2) const
    {
        // Hogen paper: H_multi Hc (d)
        
        const auto [n1, l1, R1, N1, mN1] = state1;
        const auto [n2, l2, R2, N2, mN2] = state2;

        // diagonal in N and mN
        if (N1 != N2 || mN1 != mN2)
            return 0.0;


        // if (std::abs(R1-R2) != 2)
        //     return 0.0;

        // convert from Hund's case (b) to (d)
        double result = 0;
        for (int lambda = 0; lambda <= this->GetQuantumDefectHcbMaxLambda(); lambda++)
        {
            double hcb = this->GetCoreInteractionHcb(n1, l1, R1, n2, l2, R2, lambda, N1);
            double a1 = this->GetHcbToHcdCoeff(N1, l1, R1, lambda);
            double a2 = this->GetHcbToHcdCoeff(N2, l2, R2, lambda);
            result += hcb * a1 * a2;
        }

        
        constexpr double k = -2 * PlanckConstant_v * SpeedOfLight_v;
        double RNO = this->GetScaledRydbergConstant();
        double nu1 = n1 - this->GetQuantumDefect(state1);
        double nu2 = n2 - this->GetQuantumDefect(state2);
        result *= k * RNO / std::pow(nu1*nu2, 1.5);

        return result;
    }

    double RydbergDiatomic::GetDipoleME(const RydbergDiatomicState_t& state1, const RydbergDiatomicState_t& state2) const
    {
        const auto [n1, l1, R1, N1, mN1] = state1;
        const auto [n2, l2, R2, N2, mN2] = state2;
        
        // selection rules
        if (std::abs(l1-l2) != 1 || R1 != R2 || mN1 != mN2)
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
        double rmax1 = this->GetIntegrationRange(n1);
        double rmax2 = this->GetIntegrationRange(n2);
        dip *= this->GetDipMeRadHelper(state1, state2, 1.0, rmax1, rmax2, s_defaultIntStepsPerOsc);
        dip *= std::sqrt(std::max(l1, l2)) * (l2 > l1 ? -1 : 1);

        return dip * ElementaryCharge_v;
    }

    double RydbergDiatomic::GetHcbToHcdCoeff(int N, int l, int R, int lambda)
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

    double NitricOxide::GetSDMixingAngle() const
    {
        // Phys. Chem. Chem. Phys., 2021, 23, 18806; doi: 10.1039/d1cp01930a
        return -38.7 * Pi_v / 180.0;
    }

    double NitricOxide::GetQuantumDefectHcb(int l, int Lambda) const
    {
        // DOI: 10.1039/d1cp01930a
        constexpr static std::array<double, 1> s_l0quantumDefects = { 0.210 };
        constexpr static std::array<double, 2> s_l1quantumDefects = { 0.7038, 0.7410 };
        constexpr static std::array<double, 3> s_l2quantumDefects = { 0.050, -0.053, 0.089 };
        constexpr static std::array<double, 4> s_l3quantumDefects = { 0.0182, 0.0172, 0.00128, 0.0057 };

        assert(Lambda <= l && l <= this->GetQuantumDefectHcbMaxLambda());

        const static std::array<const double*, 4> quantumDefects = {
            s_l0quantumDefects.data(), s_l1quantumDefects.data(),
            s_l2quantumDefects.data(), s_l3quantumDefects.data()
        };
        return quantumDefects[l][Lambda];
    }

    int NitricOxide::GetQuantumDefectHcbMaxLambda() const
    {
        return 3;
    }

}
