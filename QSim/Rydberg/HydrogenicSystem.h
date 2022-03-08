// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_HydrogenicSystem_H_
#define QSim_Rydberg_HydrogenicSystem_H_

#include <cstdint>
#include <utility>
#include <cmath>
#include <complex>
#include <algorithm>

#include <Eigen/Dense>

#include <iostream>
#include "../Math/Quad.h"
#include "../Math/Numerov.h"

namespace QSim
{

    class HydrogenicSystem
    {
    public:

        double GetEnergy(unsigned int n)
        {
            return -RydbergEnergy_v / (n*n);
        }

        double CorePotential(double r, unsigned int n, unsigned int l)
        {
            // k1 = hbar^2/(2*mu); k2 = e^2/(4*pi*eps0)
            constexpr double k1 = ConstexprPow(ReducedPlanckConstant_v, 2) / (2*ElectronMass_v);
            constexpr double k2 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);
            return -k2/r - k1 * l*(l+1) / (r*r);
        }

        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFLinear(unsigned int n, unsigned int l, 
            double rInner, double rOuter, std::size_t steps)
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
        }

        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWF(unsigned int n, unsigned int l, 
            double rInner, double rOuter, std::size_t steps)
        {
            // Variable transformation:
            // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
            auto [rs, rads] = GetRadialWFTransformed(n, l, std::sqrt(rInner), std::sqrt(rOuter), steps);

            // transform back
            rads = rads.cwiseQuotient(rs.cwiseProduct(rs.cwiseSqrt()));
            rs = rs.cwiseProduct(rs);

            return std::make_pair(rs, rads);
        }


    public:
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFTransformed(unsigned int n, unsigned int l, 
            double xInner, double xOuter, std::size_t steps)
        {
            // Variable transformation:
            // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
            constexpr double k1 = -2*ElectronMass_v / ConstexprPow(ReducedPlanckConstant_v, 2);
            auto kfunc = [&](double x){ 
                double r = x*x;
                return -(0.75/r + 4*r*k1*(GetEnergy(n) - CorePotential(r, n, l))); 
            };
            
            auto [xs, fs] = Numerov::Integrate(kfunc, xOuter, xInner, steps, 0.01, 0);

            // Normalization
            // int_0^\infty r^2 (R(r))^2 dr = 2 * int_0^infty x^2 (f(x))^2 = 1
            double dx = (xOuter - xInner) / (xs.size() - 1);
            double normSq = QuadSimpsonPolicy::Integrate((xs.cwiseProduct(fs)).array().square().matrix(), 2*dx);
            fs /= std::sqrt(normSq);

            return std::make_pair(xs, fs);
        }

        /*template<typename T>
        double GetRadialMatrixElementLinear(unsigned int n1, unsigned int l1, unsigned int n2, unsigned int l2, T& ax)
        {
            if (n1 > n2) return GetRadialMatrixElementLinear(n2, l2, n1, l1, ax);
            
            double dr = BohrRadius_v / 100;
            int cnt1 = static_cast<int>(std::ceil(3.5*(n1+5)*n1*BohrRadius_v / dr));
            int cnt2 = static_cast<int>(std::ceil(3.5*(n2+5)*n2*BohrRadius_v / dr));

            auto [r1, psi1] = GetRadialWFLinear(n1, l1, dr, cnt1*dr, cnt1);
            auto [r2, psi2] = GetRadialWFLinear(n2, l2, dr, cnt2*dr, cnt2);
            
            auto r1Cb = (r1.cwiseProduct(r1)).cwiseProduct(r1);
            auto overlap = psi1.cwiseProduct(psi2.tail(cnt1));
            auto integral = r1Cb.cwiseProduct(overlap);

            ax.Plot(r1.data(), psi1.data(), r1.size());
            ax.Plot(r2.tail(cnt1).eval().data(), psi2.tail(cnt1).eval().data(), r2.tail(cnt1).size());

            return QuadSimpsonPolicy::Integrate(integral, dr);
        }*/

        double GetRadialMatrixElement(unsigned int n1, unsigned int l1, unsigned int n2, unsigned int l2)
        {
            if (n1 > n2) return GetRadialMatrixElement(n2, l2, n1, l1);
            
            double dx = std::sqrt(BohrRadius_v / 1000);
            double xmax1 = std::sqrt(3*(n1+15)*n1*BohrRadius_v);
            double xmax2 = std::sqrt(3*(n2+15)*n2*BohrRadius_v);
            int cnt1 = static_cast<int>(std::ceil(xmax1 / dx));
            int cnt2 = static_cast<int>(std::ceil(xmax2 / dx));

            auto [x1, f1] = GetRadialWFTransformed(n1, l1, dx, cnt1*dx, cnt1);
            auto [x2, f2] = GetRadialWFTransformed(n2, l2, dx, cnt2*dx, cnt2);
            
            auto x1Quad = x1.array().square().square().matrix();
            auto overlap = f1.cwiseProduct(f2.tail(cnt1));
            auto integral = overlap.cwiseProduct(x1Quad);

            return 2*QuadSimpsonPolicy::Integrate(integral, dx);
        }

    };
    
}

#endif
