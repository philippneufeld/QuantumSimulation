// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_QuantumDefect_H_
#define QSim_Rydberg_QuantumDefect_H_

// Includes
#include <cmath>
#include <array>
#include <type_traits>

namespace QSim
{
    template<std::size_t N>
    double GetQuantumDefectFromCoeffs(int n, const std::array<double, N>& coeffs)
    {
        // Quantum defect expansion in n
        // d = d_0 + d_2 / (n - d_0)^2 + d_4 / (n - d_0)^4 + ...

        double defect = 0.0;
        if (N > 0)
        {
            double d0 = coeffs[0];
            defect += d0;
            if (N > 1)
            {
                double factor = 1.0 / ((n - d0)*(n - d0));
                double weight = factor;

                for (int i=1; i < N; i++, weight*=factor)
                    defect += weight * coeffs[i];
            }
        }
        return defect;
    }

    inline double ExtrapolateQuantumDefect(double defect, int l0, int l)
    {
        // d_l = d_l0 * (l0 / l)^5
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.74.062712
        return defect * std::pow(static_cast<double>(l0) / l, 5);
    }
}

#endif
