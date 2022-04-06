// Philipp Neufeld, 2021-2022

#include "Wigner.h"

#include <array>
#include <algorithm>

namespace QSim
{
    // References:
    
    //
    // Wigner 3j symbol
    //
    // References:
    // https://doi.org/10.1002/qua.560520303
    //

    double WignerLnDeltaHelper(int twoA, int twoB, int twoC)
    {
        double result = 0.0;
        result += GammaFunction::GammaLn(1 + (twoA + twoB - twoC)/2.0);
        result += GammaFunction::GammaLn(1 + (twoA - twoB + twoC)/2.0);
        result += GammaFunction::GammaLn(1 + (-twoA + twoB + twoC)/2.0);
        result -= GammaFunction::GammaLn(2 + (twoA + twoB + twoC)/2.0);
        return result;
    }

    double Wigner3jInt(int twoJ1, int twoJ2, int twoJ3, int twoM1, int twoM2, int twoM3)
    {
        
    
        // m1+m2+m3 must equal to 0
        if ((twoM1 + twoM2 + twoM3) != 0.0)
            return 0.0;

        // (j1+j2+j3) must be an integer
        int twoJsum = twoJ1 + twoJ2 + twoJ3;
        if (twoJsum % 2 != 0)
            return 0.0;

        // (j1+j2+j3) must be an even integer if (m1=m2=m3=0)
        if (twoM1 == twoM2 && twoM1 == twoM3 && twoJsum % 4 != 0)
            return 0.0;

        // js must be non-negative
        if (twoJ1 < 0 || twoJ2 < 0 || twoJ3 < 0)
            return 0.0;

        // triangle condition satisfied
        if (twoJ3 > (twoJ1 + twoJ2) || twoJ3 < std::abs(twoJ1 - twoJ2))
            return 0.0;

        // m_i in { -j_i, -j_i+1, ..., j_i-1, j_i }
        if (std::abs(twoM1) > twoJ1 || std::abs(twoM2) > twoJ2 || std::abs(twoM3) > twoJ3)
            return 0.0;

        // calculate prefactor
        double lnPre = WignerLnDeltaHelper(twoJ1, twoJ2, twoJ3);
        lnPre += GammaFunction::GammaLn(1 + (twoJ3 - twoM3)/2.0);
        lnPre += GammaFunction::GammaLn(1 + (twoJ3 + twoM3)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ2 - twoM2)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ2 + twoM2)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ1 - twoM2 - twoM3)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ1 + twoM2 + twoM3)/2.0);
        double prefactor = std::exp(0.5*lnPre);

        // calculate sum part
        double sumPart = 0.0;
        double u = (twoJ2 - twoJ1 + twoJ3) / 2.0;
        int twoZmax = std::min(twoJ3 - twoM3, twoJ2 - twoJ1 + twoJ3);
        int zmax = static_cast<int>(std::round(twoZmax / 2.0));
        for (int z=0; z<=zmax; z++)
        {
            double zterm = 0.0;
            zterm += GammaFunction::GammaLn(1 - z + (twoJ2 + twoJ3 - twoM2 - twoM3)/2.0);
            zterm += GammaFunction::GammaLn(1 + z + (twoJ1 + twoM2 + twoM3)/2.0);
            zterm -= GammaFunction::GammaLn(1 + z);
            zterm -= GammaFunction::GammaLn(1 - z + (twoJ3 - twoM3)/2.0);
            zterm -= GammaFunction::GammaLn(1 - z + (twoJ2 - twoJ1 + twoJ3)/2.0);
            zterm -= GammaFunction::GammaLn(1 + z + (twoJ1 - twoJ2 + twoM3)/2.0);
            sumPart += (z % 2 == 0 ? 1.0 : -1.0) * std::exp(zterm);
        }

        double phase = (twoJ2 - (twoJ1+twoM1)/2) % 2 == 0 ? 1.0 : -1.0;
        return phase * prefactor * sumPart;
    }

    //
    // Wigner 6j symbol
    //
    // References:
    // https://doi.org/10.1016/S0010-4655(98)00065-4
    // https://doi.org/10.1137/15M1021908
    //

    bool Wigner6jCheckTriple(int twoJ1, int twoJ2, int twoJ3)
    {
        // sum must be integer
        if ((twoJ1 + twoJ2 + twoJ3) % 2 == 1)
            return false;

        if (twoJ1 < std::abs(twoJ2-twoJ3) || twoJ1 > (twoJ2 + twoJ3))
            return false;

        return true;
    }

    double Wigner6jInt(int twoJ1, int twoJ2, int twoJ3, int twoJ4, int twoJ5, int twoJ6)
    { 
        
    
        // triangle conditions
        if (!Wigner6jCheckTriple(twoJ1, twoJ2, twoJ3) || 
            !Wigner6jCheckTriple(twoJ4, twoJ5, twoJ3) ||
            !Wigner6jCheckTriple(twoJ1, twoJ5, twoJ6) || 
            !Wigner6jCheckTriple(twoJ4, twoJ2, twoJ6)) 
            return 0.0;
        
        std::array<int, 4> alpha = {
            twoJ1 + twoJ2 + twoJ3,
            twoJ4 + twoJ5 + twoJ3,
            twoJ1 + twoJ5 + twoJ6,
            twoJ4 + twoJ2 + twoJ6
        };
        
        std::array<int, 3> beta = {
            twoJ1 + twoJ2 + twoJ4 + twoJ5,
            twoJ1 + twoJ3 + twoJ4 + twoJ6,
            twoJ2 + twoJ3 + twoJ5 + twoJ6
        };

        // alphas and betas must be integer
        if (alpha[0]%2 || alpha[1]%2 || alpha[2]%2 || alpha[3]%2 ||
            beta[0]%2 || beta[1]%2 || beta[2]%2)
            return 0.0;

        for (auto& v: alpha) v /= 2;
        for (auto& v: beta) v /= 2;
        
        double lnPre = 0.0;
        lnPre += WignerLnDeltaHelper(twoJ1, twoJ2, twoJ3);
        lnPre += WignerLnDeltaHelper(twoJ4, twoJ5, twoJ3);
        lnPre += WignerLnDeltaHelper(twoJ1, twoJ5, twoJ6);
        lnPre += WignerLnDeltaHelper(twoJ4, twoJ2, twoJ6);
        double pre = std::exp(0.5 * lnPre);

        double sumPart = 0.0;
        int kmin = *std::max_element(alpha.begin(), alpha.end());
        int kmax = *std::min_element(beta.begin(), beta.end());

        for (int k=kmin; k<=kmax; k++)
        {
            double kterm = GammaFunction::GammaLn(2 + k);
            for (auto& v: alpha) kterm -= GammaFunction::GammaLn(1 + k - v);
            for (auto& v: beta) kterm -= GammaFunction::GammaLn(1 - k + v);
            sumPart += (k % 2 == 1 ? -1 : 1) * std::exp(kterm);
        }

        return pre * sumPart;
    }

    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
    {
        // round because of inprecise floating point arithmetic
        // angular momenta can only be an integer or half-integer
        int twoJ1 = static_cast<int>(std::round(2*j1));
        int twoJ2 = static_cast<int>(std::round(2*j2));
        int twoJ3 = static_cast<int>(std::round(2*j3));
        int twoM1 = static_cast<int>(std::round(2*m1));
        int twoM2 = static_cast<int>(std::round(2*m2));
        int twoM3 = static_cast<int>(std::round(2*m3));

        return Wigner3jInt(twoJ1, twoJ2, twoJ3, twoM1, twoM2, twoM3);
    }

    double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6)
    {
        // round because of inprecise floating point arithmetic
        // angular momenta can only be an integer or half-integer
        int twoJ1 = static_cast<int>(std::round(2*j1));
        int twoJ2 = static_cast<int>(std::round(2*j2));
        int twoJ3 = static_cast<int>(std::round(2*j3));
        int twoJ4 = static_cast<int>(std::round(2*j4));
        int twoJ5 = static_cast<int>(std::round(2*j5));
        int twoJ6 = static_cast<int>(std::round(2*j6));

        return Wigner6jInt(twoJ1, twoJ2, twoJ3, twoJ4, twoJ5, twoJ6);
    }

    double ClebshGordan(double j1, double j2, double j3, double m1, double m2, double m3)
    {
        double phase = static_cast<int>(j2-j1-m3) % 2 == 0 ? 1.0 : -1.0;
        return phase * std::sqrt(2*j3+1) * Wigner3j(j1, j2, j3, m1, m2, -m3);
    }

}
