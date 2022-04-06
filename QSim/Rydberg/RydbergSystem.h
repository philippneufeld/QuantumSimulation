// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_RydbergSystem_H_
#define QSim_Rydberg_RydbergSystem_H_

#include <cstdint>

#include <Eigen/Dense>

#include "../Constants.h"
#include "../Math/Quad.h"
#include "../Math/Numerov.h"
#include "../Math/Wigner.h"

namespace QSim
{

    template<typename State>
    class TRydbergSystem
    {
    public:
        TRydbergSystem(double mass);
        virtual ~TRydbergSystem();

        // Public interface -> functions must be overwritten by child class
        virtual double GetQuantumDefect(const State& state) const = 0;
        virtual double GetEnergy(const State& state) const = 0;
        virtual double GetPotential(double r, const State& state) const = 0;
        virtual double GetDipoleME(const State& state1, const State& state2) const = 0;

    protected:
        // Rydberg energy helpers
        double GetScaledRydbergConstant() const;
        double GetRydbergEnergy(int n, const State& state) const;

        // potential helpers
        double GetCoulombPotential(double r, int n, int l) const;
        double GetAtomicPotential(double r, int n, int l) const;
        double GetAtomicPotentialFS(double r, int n, int l, double j) const;

        // wavefunction in transformed variables (see function definition)
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFTransformed(
            const State& state, double xInner, double xOuter, 
            std::size_t steps, std::size_t peakStepThreshold) const;

        // Helper for the calculation of the radial dipole matrix elements
        double GetDipMeRadHelper(
            const State& state1, const State& state2,
            double rmax1, double rmax2, std::size_t stepsPerOscillation) const;

    private:
        double m_reducedMass;
    };

    namespace Internal
    {
        // Class that detects if the Numerov algorithm diverges
        class RydSysNumerovMonitor
        {
        public:
            RydSysNumerovMonitor(int distThreshold) 
                : m_distThreshold(distThreshold), m_maxIdx(0), m_maxVal(0.0) {};

            template<typename Xs, typename Ys>
            bool operator()(int i, Xs& xs, Ys& ys);

        private:
            const int m_distThreshold;
            int m_maxIdx;
            double m_maxVal;
        };
    }


    //
    // Function definitions
    //

    template<typename State>
    TRydbergSystem<State>::TRydbergSystem(double mass)
    {
        m_reducedMass = (mass * ElectronMass_v) / (mass + ElectronMass_v);
    }

    template<typename State>
    TRydbergSystem<State>::~TRydbergSystem() { }
    
    template<typename State>
    double TRydbergSystem<State>::GetScaledRydbergConstant() const
    {
        return RydbergConstant_v * (m_reducedMass / ElectronMass_v);
    }

    template<typename State>
    double TRydbergSystem<State>::GetRydbergEnergy(int n, const State& state) const
    {
        constexpr double hc = PlanckConstant_v * SpeedOfLight_v;
        double nAdj = n - this->GetQuantumDefect(state);
        return -hc * this->GetScaledRydbergConstant() / (nAdj*nAdj);
    }

    template<typename State>
    double TRydbergSystem<State>::GetAtomicPotential(double r, int n, int l) const
    {
        double potential = 0.0;
        double r2 = r * r;

        // Electrostatic potential
        // k1 = e^2/(4*pi*eps0)
        constexpr double k1 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);
        potential -= k1 / r;

        // Orbital potential term
        // k2 = hbar^2/(2*mu)
        constexpr double ck2 = ConstexprPow(ReducedPlanckConstant_v, 2) / 2;
        const double k2 = ck2 / m_reducedMass;
        potential += k2 * l*(l+1) / r2;

        return potential;
    }

    template<typename State>
    double TRydbergSystem<State>::GetAtomicPotentialFS(double r, int n, int l, double j) const
    {
        // Atomic potential (electrostatic and orbital term)
        double potential = GetAtomicPotential(r, n, l);
        double r3 = r * r * r;

        // Fine-structure (L-S coupling term)
        // k3 = alpha * hbar^3 / (4*me^2*c) = e^2 / (4 pi eps0) * (gs-1) / (4*me^2*c^2)
        constexpr double k3 = FineStructureConstant_v * ConstexprPow(ReducedPlanckConstant_v, 3) / 
            (4*ConstexprPow(ElectronMass_v, 2)*SpeedOfLight_v);
        potential += k3 / r3 * (j*(j+1) - l*(l+1) - 0.75);

        return potential;
    }

    template<typename State>
    std::pair<Eigen::VectorXd, Eigen::VectorXd> TRydbergSystem<State>::GetRadialWFTransformed(
        const State& state, double xInner, double xOuter, 
        std::size_t steps, std::size_t peakStepThreshold) const
    {
        // Variable transformation:
        // f(x) = r^(3/4)*R(r)  with  x = sqrt(r)
        
        // k1 = - 2*mu / hbar^2
        constexpr double ck1 = -2 / ConstexprPow(ReducedPlanckConstant_v, 2);
        const double k1 = m_reducedMass * ck1;

        auto kfunc = [&](double x){ 
            double r = x*x;
            return -(0.75/r + 4*r*k1*(this->GetEnergy(state) - this->GetPotential(r, state))); 
        };
        
        // Do Numerov integration (with divergence control)
        Internal::RydSysNumerovMonitor monitor(peakStepThreshold);
        auto [xs, fs] = Numerov::Integrate(kfunc, xOuter, xInner, steps, 0.01, 0, monitor);

        // Normalization
        // int_0^\infty r^2 (R(r))^2 dr = 2 * int_0^infty x^2 (f(x))^2 = 1
        double dx = (xOuter - xInner) / (xs.size() - 1);
        double normSq = QuadSimpsonPolicy::Integrate((xs.cwiseProduct(fs)).array().square().matrix(), 2*dx);
        fs /= std::sqrt(normSq);

        return std::make_pair(xs, fs);
    }

    template<typename State>
    double TRydbergSystem<State>::GetDipMeRadHelper(
        const State& state1, const State& state2,
        double rmax1, double rmax2, std::size_t stepsPerOscillation) const
    {
        if (rmax1 > rmax2)
            return this->GetDipMeRadHelper(state2, state1, rmax2, rmax1, stepsPerOscillation);

        double dx = std::sqrt(BohrRadius_v / (stepsPerOscillation * stepsPerOscillation));
        int cnt1 = static_cast<int>(std::ceil(std::sqrt(rmax1) / dx));
        int cnt2 = static_cast<int>(std::ceil(std::sqrt(rmax2) / dx));

        std::size_t peakStepThreshold = std::max<std::size_t>(stepsPerOscillation / 4, 3);

        auto [x1, f1] = GetRadialWFTransformed(state1, dx, cnt1*dx, cnt1, peakStepThreshold);
        auto [x2, f2] = GetRadialWFTransformed(state2, dx, cnt2*dx, cnt2, peakStepThreshold);
        
        auto x1Quad = x1.array().square().square().matrix();
        auto overlap = f1.cwiseProduct(f2.tail(cnt1));
        auto integrand = overlap.cwiseProduct(x1Quad);

        return 2 * QuadSimpsonPolicy::Integrate(integrand, dx);
    }

    namespace Internal
    {
        template<typename Xs, typename Ys>
        bool RydSysNumerovMonitor::operator()(int i, Xs& xs, Ys& ys)
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

}

#endif
