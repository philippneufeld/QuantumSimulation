// Philipp Neufeld, 2021-2022

#include "RydbergAtom.h"
#include "QuantumDefects.h"

namespace QSim
{

    //
    // RydbergAtom
    //

    RydbergAtom::RydbergAtom(double mass) 
        : TRydbergSystem<RydbergAtomState_t>(mass) {}

    double RydbergAtom::GetEnergy(const RydbergAtomState_t& state) const
    {
        const auto [n, l, j, mj] = state;
        return this->GetRydbergEnergy(n, state);
    }

    double RydbergAtom::GetPotential(double r, const RydbergAtomState_t& state) const
    {
        const auto [n, l, j, mj] = state;
        return this->GetAtomicPotentialFS(r, n, l, j);
    }

    double RydbergAtom::GetDipoleME(const RydbergAtomState_t& state1, const RydbergAtomState_t& state2) const
    {
        const auto [n1, l1, j1, mj1] = state1;
        const auto [n2, l2, j2, mj2] = state2;
        
        double dip = this->GetDipMEAngular(l1, j1, mj1, l2, j2, mj2);
        
        if (dip != 0)
        {
            // radial dipole matrix element
            double rmax1 = 3*(n1+15)*n1*BohrRadius_v;
            double rmax2 = 3*(n2+15)*n2*BohrRadius_v;
            dip *= this->GetDipMeRadHelper(state1, state2, rmax1, rmax2, 50);
        }

        return dip * ElementaryCharge_v;
    }

    double RydbergAtom::GetDipMEAngular(
        int l1, double j1, double mj1, int l2, double j2, double mj2) const
    {
        // PRA 20.6 (1979)
        // <l,m| cos \Theta |l-1, m> = sqrt((l^2-m^2)/((2*l+1)*(2*l-1)))
        // Then summation over the clebsch gordan coefficients

        if (std::round(2*mj1) != std::round(2*mj2))
            return 0.0;

        double result = 0.0;
        double s = 0.5;
        for(int i = 0; i<static_cast<int>(std::round(2*s+1)); i++)
        {
            double ml = mj1 - s + i;
            if (std::abs(ml)-0.1 < l1 && std::abs(ml)-0.1 < l2)
            {
                int twoL1 = static_cast<int>(std::round(2*l1));
                int twoL2 = static_cast<int>(std::round(2*l2));
                if (std::abs(twoL1 - twoL2) == 2)
                {
                    double lmax = std::max(l1, l2);
                    double ang = std::sqrt((lmax*lmax - ml*ml) / ((2*lmax+1)*(2*lmax-1)));
                    double cg1 = ClebshGordan(l1, s, j1, ml, mj1 - ml, mj1);
                    double cg2 = ClebshGordan(l2, s, j2, ml, mj2 - ml, mj2);
                    result += cg1 * cg2 * ang;
                }
            }
        }

        return result;
    }

    //
    // Hydrogen
    //

    Hydrogen::Hydrogen() 
        : RydbergAtom(AtomicMassUnit_v) {}

    double Hydrogen::GetQuantumDefect(const RydbergAtomState_t& state) const
    {
        return 0.0;
    }

    //
    // Rubidium
    //

    Rubidium::Rubidium() 
        : RydbergAtom(84.9117897379*AtomicMassUnit_v) {}

    double Rubidium::GetQuantumDefect(const RydbergAtomState_t& state) const
    {
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.67.052502 (l=0...2)
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.74.054502 (l=3)
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.102.062817 (l=4)

        const auto [n, l, j, mj] = state;

        using DefectCoeffs = std::array<double, 2>;
        constexpr DefectCoeffs coeffs1[] = { // coeffs for j = l - s
            {3.1311804, 0.1784}, // s
            {2.6548849, 0.2900}, // p
            {1.34809171, -0.60286}, // d
            {0.0165192, -0.085}, // f
            {0.0039990, -0.0202}, // g
        };
        constexpr DefectCoeffs coeffs2[] = { // coeffs for j = l + s
            {3.1311804, 0.1784}, // s
            {2.6416737, 0.2950}, // p
            {1.34646572, -0.59600}, // d
            {0.0165437, -0.086}, // f
            {0.0039990, -0.0202}, // g
        };

        // select right table
        double s = 0.5;
        int idx = static_cast<int>(std::round(((j - l) / s + 1) * 0.5));

        const DefectCoeffs* pDefectTable = nullptr;
        if (idx == 0) pDefectTable = coeffs1;
        else if(idx == 1) pDefectTable = coeffs2;
        else return 0.0;

        int ldef = std::min(l, 4);
        double defect = GetQuantumDefectFromCoeffs(n, pDefectTable[ldef]);

        if (l > 4)
            defect = ExtrapolateQuantumDefect(defect, ldef, l);
        
        return defect;
    }

}
