// Philipp Neufeld, 2021-2022

#include "HydrogenicSystem.h"

namespace QSim
{
    
    double HydrogenicSystem::GetEnergy(int n) const
    {
        return -RydbergEnergy_v / (n*n);
    }

    double HydrogenicSystem::CorePotential(double r, int n, int l) const
    {
        // k1 = hbar^2/(2*mu); k2 = e^2/(4*pi*eps0)
        constexpr double k1 = ConstexprPow(ReducedPlanckConstant_v, 2) / (2*ElectronMass_v);
        constexpr double k2 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);
        return -k2/r + k1 * l*(l+1) / (r*r);
    }

    double HydrogenicSystem::FSPotential(double r, int n, int l, double j) const
    {
        return CorePotential(r, n, l);

        // TODO
        constexpr double k1 = 1.0 / (2*ConstexprPow(ElectronMass_v, 2)*ConstexprPow(SpeedOfLight_v, 2));
        constexpr double k2 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);
        constexpr double k = k1*k2 * ConstexprPow(ReducedPlanckConstant_v, 2) / 2;

        double r3 = r*r*r;
        return CorePotential(r, n, l) + k / r3 * (j*(j+1) - l*(l+1) - 0.75);
    }

    /*std::pair<Eigen::VectorXd, Eigen::VectorXd> HydrogenicSystem::GetRadialWFLinear(int n, int l, 
        double rInner, double rOuter, std::size_t steps) const
    {
        // P(r) = r*R(r)
        // Solve P''(r) = (E-V(r))*P(r)
        constexpr double k1 = -2*ElectronMass_v / ConstexprPow(ReducedPlanckConstant_v, 2);
        auto kfunc = [&](double r){ return -k1*(GetEnergy(n) - CorePotential(r, n, l)); };
        auto [rs, rads] = Numerov::Integrate(kfunc, rOuter, rInner, steps, 0.01, 0);

        // Normalization
        // int_0^\infty r^2 |R(r)|^2 dr = int_0^\infty |P(r)|^2 dr
        double dr = (rOuter - rInner) / (rs.size() - 1);
        double tmp = rads.norm() * std::sqrt(dr);

        // P(r) => R(r)
        rads = rads.cwiseQuotient(tmp * rs); 

        return std::make_pair(rs, rads);
    }*/

    std::pair<Eigen::VectorXd, Eigen::VectorXd> HydrogenicSystem::GetRadialWF(int n, int l, double j,
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

    double HydrogenicSystem::GetDipoleMERad(int n1, int l1, double j1, int n2, int l2, double j2) const
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

    double HydrogenicSystem::GetDipoleMEAng(int l1, int m1, int l2, int m2) const
    {
        // sqrt(4pi/3) * int dOm Y_{l1,m1}(Om) * Y_{l2,m2}(Om) * Y_{1,0}(Om)

        // selection rules
        if (std::abs(l2-l1) != 1 || std::abs(m2-m1) > 1)
            return 0.0;

        double res = std::sqrt((2*l1+1)*(2*l2+1));
        res *= Wigner3j(l1, l2, 1, m1, m2, 0);

        if (res != 0.0)
            res *= Wigner3j(l1, l2, 1, 0, 0, 0);
            
        return res;
    }
    
    double HydrogenicSystem::GetDipoleMEAng2(int l1, double j1, double mj1, int l2, double j2, double mj2) const
    {
        double result = 0.0;
        int multiplicity = 2; // 2*s+1
        for(int i = 0; i<multiplicity; i++)
        {
            double ml = mj1 - 0.5 + i;
            if (std::abs(ml)-0.1 < l1 && std::abs(ml)-0.1 < l2)
            {
                double ang = 0.0;
                if (std::abs(l1 - l2 - 1) < 0.1)
                {
                    ang = std::sqrt((l1*l1 - ml*ml) / ((2*l1+1)*(2*l1-1)));
                }
                else if (std::abs(l1 - l2 + 1) < 0.1)
                {
                    ang = std::sqrt((l2*l2 - ml*ml) / ((2*l2+1)*(2*l2-1)));
                }
                else
                {
                    continue;
                }
                
                double cg1 = ClebshGordan(l1, 0.5, j1, ml, mj1 - ml, mj1);
                double cg2 = ClebshGordan(l2, 0.5, j2, ml, mj2 - ml, mj2);
                result += cg1 * cg2 * ang;
            }
        }

        return result;
    }

    /*double HydrogenicSystem::GetDipoleME(int n1, int l1, int m1, int n2, int l2, int m2) const
    {
        double res = GetDipAngularME(l1, m1, l2, m2);
        if (res != 0)
            res *= GetDipRadialME(n1, l1, n2, l2);
        return res * ElementaryCharge_v;
    }*/

    double HydrogenicSystem::GetDipoleME2(int n1, int l1, double j1, double m1, int n2, int l2, double j2, double m2) const
    {
        double res = GetDipoleMEAng2(l1, j1, m1, l2, j2, m2);
        if (res != 0)
            res *= GetDipoleMERad(n1, l1, j1, n2, l2, j2);
        return res * ElementaryCharge_v;
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd> HydrogenicSystem::GetRadialWFTransformed(int n, int l, double j,
        double xInner, double xOuter, std::size_t steps) const
    {
        // Variable transformation:
        // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
        constexpr double k1 = -2*ElectronMass_v / ConstexprPow(ReducedPlanckConstant_v, 2);
        auto kfunc = [&](double x){ 
            double r = x*x;
            return -(0.75/r + 4*r*k1*(GetEnergy(n) - FSPotential(r, n, l, j))); 
        };
        
        Internal::HydSysNumerovMonitor monitor(steps / (4*n));
        auto [xs, fs] = Numerov::Integrate(kfunc, xOuter, xInner, steps, 0.01, 0, monitor);

        // Normalization
        // int_0^\infty r^2 (R(r))^2 dr = 2 * int_0^infty x^2 (f(x))^2 = 1
        double dx = (xOuter - xInner) / (xs.size() - 1);
        double normSq = QuadSimpsonPolicy::Integrate((xs.cwiseProduct(fs)).array().square().matrix(), 2*dx);
        fs /= std::sqrt(normSq);

        return std::make_pair(xs, fs);
    }

}
