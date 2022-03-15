// Philipp Neufeld, 2021-2022

#include "Wigner.h"

namespace QSim
{

    double Wigner3jInt(int twoJ1, int twoJ2, int twoJ3, int twoM1, int twoM2, int twoM3)
    {
        // Reference:
        // Shan-Tao Lai, 
        // Computation of Algebraic Formulas for Wigner 
        // 3-j, 6-j, and 9-j Symbols by Maple

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

        double lnDelta = 0.0;
        lnDelta += GammaFunction::GammaLn(1 + (twoJ1 + twoJ2 - twoJ3)/2.0);
        lnDelta += GammaFunction::GammaLn(1 + (twoJ1 - twoJ2 + twoJ3)/2.0);
        lnDelta += GammaFunction::GammaLn(1 + (-twoJ1 + twoJ2 + twoJ3)/2.0);
        lnDelta -= GammaFunction::GammaLn(2 + (twoJ1 + twoJ2 + twoJ3)/2.0);
        
        double lnPre = 0.0;
        lnPre += GammaFunction::GammaLn(1 + (twoJ3 - twoM3)/2.0);
        lnPre += GammaFunction::GammaLn(1 + (twoJ3 + twoM3)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ2 - twoM2)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ2 + twoM2)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ1 - twoM2 - twoM3)/2.0);
        lnPre -= GammaFunction::GammaLn(1 + (twoJ1 + twoM2 + twoM3)/2.0);

        double prefactor = std::exp(0.5*(lnDelta + lnPre));

        double sumPart = 0.0;
        double u = (twoJ2 - twoJ1 + twoJ3) / 2.0;
        double zmax = std::min(twoJ3 - twoM3, twoJ2 - twoJ1 + twoJ3) / 2.0;
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

    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
    {
        return Wigner3jInt(
            static_cast<int>(2*j1), static_cast<int>(2*j2), static_cast<int>(2*j3),
            static_cast<int>(2*m1), static_cast<int>(2*m2), static_cast<int>(2*m3));
    }

    double ClebshGordan(double j1, double j2, double j3, double m1, double m2, double m3)
    {
        double phase = static_cast<int>(j2-j1-m3) % 2 == 0 ? 1.0 : -1.0;
        return phase * std::sqrt(2*j3+1) * Wigner3j(j1, j2, j3, m1, m2, -m3);
    }

}
