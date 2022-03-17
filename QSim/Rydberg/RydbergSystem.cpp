// Philipp Neufeld, 2021-2022

#include <utility>
#include <cmath>
#include <complex>
#include <algorithm>

#include "RydbergSystem.h"

namespace QSim
{
    //
    // Class that detects if the Numerov algorithm diverges
    //
    class HydSysNumerovMonitor
    {
    public:
        HydSysNumerovMonitor(int distThreshold) 
            : m_distThreshold(distThreshold), m_maxIdx(0), m_maxVal(0.0) {};

        template<typename Xs, typename Ys>
        bool operator()(int i, Xs& xs, Ys& ys);

    private:
        const int m_distThreshold;
        int m_maxIdx;
        double m_maxVal;
    };


    //
    // RydbergSystem methods
    //

    RydbergSystem::RydbergSystem(double mass)
    {
        m_reducedMass = (mass * ElectronMass_v) / (mass + ElectronMass_v);

    }

    double RydbergSystem::GetQuantumDefectFromCoeffs(int n, const double* coeffs, int coeffCnt)
    {
        double defect = 0.0;
        if (coeffCnt > 0)
        {
            double d0 = coeffs[0];
            defect += d0;
            if (coeffCnt > 1)
            {
                double factor = 1.0 / ((n - d0)*(n - d0));
                double weight = factor;

                for (int i=1; i<coeffCnt; i++, weight*=factor)
                    defect += weight * coeffs[i];
            }
        }
        return defect;
    }

    double RydbergSystem::GetExtrapolatedQuantumDefect(double defect, int l, int lnew)
    {
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.74.062712
        return defect * std::pow(static_cast<double>(l) / lnew, 5);
    }

    double RydbergSystem::GetEnergy(int n, int l, double j) const
    {
        double nAdj = n - GetQuantumDefect(n, l, j);
        return -RydbergEnergy_v / (nAdj*nAdj);
    }

    double RydbergSystem::CorePotential(double r, int n, int l) const
    {
        // k1 = e^2/(4*pi*eps0); k2 = hbar^2/(2*mu)
        constexpr double k1 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);
        
        constexpr double k2Helper = ConstexprPow(ReducedPlanckConstant_v, 2) / 2;
        double k2 = k2Helper / m_reducedMass;

        return -k1/r + k2 * l*(l+1) / (r*r);
    }

    double RydbergSystem::FSPotential(double r, int n, int l, double j) const
    {
        // fs = alpha * hbar^3 / (4*me^2*c) = e^2 / (4 pi eps0) * (gs-1) / (4*me^2*c^2)
        constexpr double ls = FineStructureConstant_v * ConstexprPow(ReducedPlanckConstant_v, 3) / 
            (4*ConstexprPow(ElectronMass_v, 2)*SpeedOfLight_v);

        double r3 = r*r*r;
        return CorePotential(r, n, l) + ls / r3 * (j*(j+1) - l*(l+1) - 0.75);
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd> RydbergSystem::GetRadialWF(int n, int l, double j,
        double rInner, double rOuter, std::size_t steps) const
    {
        // Variable transformation:
        // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
        auto [rs, rads] = GetRadialWFTransformed(n, l, j, std::sqrt(rInner), std::sqrt(rOuter), steps);

        // transform back
        rads = rads.cwiseQuotient(rs.cwiseProduct(rs.cwiseSqrt()));
        rs = rs.cwiseProduct(rs);

        return std::make_pair(rs, rads);
    }

    double RydbergSystem::GetDipoleMERad(int n1, int l1, double j1, int n2, int l2, double j2) const
    {
        if (n1 > n2) return GetDipoleMERad(n2, l2, j2, n1, l1, j1);
        
        double dx = std::sqrt(BohrRadius_v / 1000);
        double xmax1 = std::sqrt(3*(n1+15)*n1*BohrRadius_v);
        double xmax2 = std::sqrt(3*(n2+15)*n2*BohrRadius_v);
        int cnt1 = static_cast<int>(std::ceil(xmax1 / dx));
        int cnt2 = static_cast<int>(std::ceil(xmax2 / dx));

        auto [x1, f1] = GetRadialWFTransformed(n1, l1, j1, dx, cnt1*dx, cnt1);
        auto [x2, f2] = GetRadialWFTransformed(n2, l2, j2, dx, cnt2*dx, cnt2);
        
        auto x1Quad = x1.array().square().square().matrix();
        auto overlap = f1.cwiseProduct(f2.tail(cnt1));
        auto integrand = overlap.cwiseProduct(x1Quad);

        return 2*QuadSimpsonPolicy::Integrate(integrand, dx);
    }

    double RydbergSystem::GetDipoleMEAng(int l1, double j1, double mj1, int l2, double j2, double mj2) const
    {
        // PRA 20.6 (1979)
        // <l,m| cos \Theta |l-1, m> = sqrt((l^2-m^2)/((2*l+1)*(2*l-1)))

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

    double RydbergSystem::GetDipoleME(int n1, int l1, double j1, double m1, int n2, int l2, double j2, double m2) const
    {
        double res = GetDipoleMEAng(l1, j1, m1, l2, j2, m2);
        if (res != 0)
            res *= GetDipoleMERad(n1, l1, j1, n2, l2, j2);
        return res * ElementaryCharge_v;
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd> RydbergSystem::GetRadialWFTransformed(int n, int l, double j,
        double xInner, double xOuter, std::size_t steps) const
    {
        // Variable transformation:
        // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
        constexpr double k1 = -2*ElectronMass_v / ConstexprPow(ReducedPlanckConstant_v, 2);
        auto kfunc = [&](double x){ 
            double r = x*x;
            return -(0.75/r + 4*r*k1*(GetEnergy(n, l, j) - FSPotential(r, n, l, j))); 
        };
        
        HydSysNumerovMonitor monitor(steps / (4*n));
        auto [xs, fs] = Numerov::Integrate(kfunc, xOuter, xInner, steps, 0.01, 0, monitor);

        // Normalization
        // int_0^\infty r^2 (R(r))^2 dr = 2 * int_0^infty x^2 (f(x))^2 = 1
        double dx = (xOuter - xInner) / (xs.size() - 1);
        double normSq = QuadSimpsonPolicy::Integrate((xs.cwiseProduct(fs)).array().square().matrix(), 2*dx);
        fs /= std::sqrt(normSq);

        return std::make_pair(xs, fs);
    }


    template<typename Xs, typename Ys>
    bool HydSysNumerovMonitor::operator()(int i, Xs& xs, Ys& ys)
    {
        double absYi = std::abs(ys[i]);
        if (i - m_maxIdx < m_distThreshold)
        {
            // get height of first peak
            if (absYi > m_maxVal)
            {
                m_maxVal = absYi;
                m_maxIdx = i;
            }
        }
        else
        {
            // subsequent peaks may not be higher than the first peak
            if (absYi > m_maxVal)
            {
                // clean up (set to zero until node where divergence started)
                for(;i > 1 && std::abs(ys[i]) > std::abs(ys[i-1]); i--);
                while(i < ys.size()) ys[i++] = 0;
                return false;
            }
        }

        return true;
    }

}
