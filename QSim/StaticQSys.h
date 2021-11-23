// Philipp Neufeld, 2021

#ifndef QSim_StaticQSys_H_
#define QSim_StaticQSys_H_

#include <cstdint>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>

#include "Math/Matrix.h"
#include "Doppler.h"
#include "DensityMatrix.h"

namespace QSim
{
    constexpr static double Pi_v = 3.14159265358979323846;
    constexpr static double TwoPi_v = 2 * Pi_v;
    constexpr static double SpeedOfLight2_v = 2.99792458e8;
    constexpr static double VacuumPermittivity_v = 8.8541878128e-12;
    constexpr static double ReducedPlanckConstant_v = 1.054571817e-34;

    //
    // N-level quantum system solver
    //

    
    struct EulerIntegrator
    {
        template<typename XTy, typename YTy, typename Func>
        static void Step(YTy& y, XTy x, XTy dx, Func func)
        {
            y += func(x, y) * dx;
        }
    };

    struct RK4Integrator
    {
        template<typename XTy, typename YTy, typename Func>
        static void Step(YTy& y, XTy x, XTy dx, Func func)
        {
            double dx2 = dx/2;
            auto k1 = func(x, y);
            auto k2 = func(x + dx2, y + dx2*k1);
            auto k3 = func(x + dx2, y + dx2*k2);
            auto k4 = func(x + dx, y + dx*k3);
            y += dx/6 * (k1 + 2*k2 + 2*k3 + k4);
        }
    };





    namespace Internal
    {
        struct StaticQSysDecayDesc
        {
            std::size_t from;
            std::size_t to;
            double rate;
        };
    }

    template<std::size_t N>
    class TStaticQSys
    {
        template<typename InputIt>
        using EnableIfLvlIt_t = std::enable_if_t<
            std::is_same<std::string, std::decay_t<decltype(std::declval<InputIt>()->first)>>::value &&
            std::is_same<double, std::decay_t<decltype(std::declval<InputIt>()->second)>>::value>;
    public:
        // constructors
        TStaticQSys(const std::map<std::string, double>& levels, double mass)
            : TStaticQSys(levels.begin(), mass) { assert(levels.size() == N); }
        template<typename InputIt, typename=EnableIfLvlIt_t<InputIt>>
        TStaticQSys(InputIt levelIterator, double mass);

        // copy operations
        TStaticQSys(const TStaticQSys&) = default;
        TStaticQSys& operator=(const TStaticQSys&) = default;

        // level name control
        bool HasLevel(const std::string& name) const { return m_levelNames.find(name) != m_levelNames.end(); }
        std::size_t GetLevelIndexByName(const std::string& name) const { return m_levelNames.at(name); }

        // functions to change the system properties 
        bool SetDipoleMatrixElement(const std::string& lvl1, const std::string& lvl2, double dip);
        bool AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rate);

        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetHamiltonian(
            const TColVector<VT>& laserFreqs,
            const TColVector<VT>& laserIntensities, 
            double velocity, double t) const;

        // thermal environment
        void SetMass(double mass) { m_doppler.SetMass(mass); }
        void SetTemperature(double temp) { m_doppler.SetTemperature(temp); }

        double GetMass() const { return m_doppler.GetMass(); }
        double GetTemperature() const { return m_doppler.GetTemperature(); }

        TStaticDensityMatrix<N> MakeGroundState() const;

        template<typename VT>
        std::vector<TStaticDensityMatrix<N>> GetTrajectoryNatural(
            const TColVector<VT>& laserFreqs,
            const TColVector<VT>& laserIntensities, 
            const TStaticDensityMatrix<N>& initial, 
            double dt, std::size_t steps);

    private:
        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetDensityOpDerivative(
            const TStaticMatrix<std::complex<double>, N, N>& rho,
            const TColVector<VT>& laserFreqs,
            const TColVector<VT>& laserIntensities, 
            double velocity, double t) const;

    private:
        // Map from the level name to their index
        std::map<std::string, std::size_t> m_levelNames;

        // Properties of the system
        TStaticColVector<double, N> m_levels;
        std::vector<Internal::StaticQSysDecayDesc> m_decays;
        TStaticMatrix<std::complex<double>, N, N> m_dipoleOperator;

        // thermal environment
        TDopplerIntegrator<double> m_doppler;
    };

    template<std::size_t N>
    template<typename InputIt, typename>
    TStaticQSys<N>::TStaticQSys(InputIt levelIterator, double mass)
        : m_doppler(mass, 300.0) // use room temperature as default
    {
        // levelIterator is a pair containing (name, level)
        for (std::size_t i = 0; i < N; i++, levelIterator++)
        {
            m_levelNames[levelIterator->first] = i;
            m_levels[i] = levelIterator->second;
        }
    }

    template<std::size_t N>
    bool TStaticQSys<N>::SetDipoleMatrixElement(const std::string& lvl1, const std::string& lvl2, double dip)
    {
        if (!HasLevel(lvl1) || !HasLevel(lvl2) || lvl1 == lvl2)
            return false;

        std::size_t idx1 = GetLevelIndexByName(lvl1);
        std::size_t idx2 = GetLevelIndexByName(lvl2);
        m_dipoleOperator(idx1, idx2) = dip;
        m_dipoleOperator(idx2, idx1) = dip;
        return true;
    }

    template<std::size_t N>
    bool TStaticQSys<N>::AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rate)
    {
        if (!HasLevel(lvlFrom) || !HasLevel(lvlTo) || lvlFrom == lvlTo)
            return false;

        Internal::StaticQSysDecayDesc decay;
        decay.from = GetLevelIndexByName(lvlFrom);
        decay.to = GetLevelIndexByName(lvlTo);
        decay.rate = rate;

        m_decays.push_back(decay);
        return true;
    }

    template<std::size_t N>
    TStaticDensityMatrix<N> TStaticQSys<N>::MakeGroundState() const
    {
        TStaticDensityMatrix<N> rho(m_levelNames);
        auto gs = std::min_element(m_levels.Data(), m_levels.Data() + N) - m_levels.Data();
        rho(gs, gs) = 1.0;
        return rho;
    }

    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TStaticQSys<N>::GetHamiltonian(
        const TColVector<VT>& laserFreqs,
        const TColVector<VT>& laserIntensities, 
        double velocity, double t) const
    {
        assert((~laserFreqs).Rows() == (~laserIntensities).Rows());
        
        using HamiltonianType = TStaticMatrix<std::complex<double>, N, N>;
        HamiltonianType hamiltonian(N, N);

        // Atom hamiltonian
        for (std::size_t i = 0; i < N; i++)
            hamiltonian(i, i) = TwoPi_v * m_levels[i];  

        // Calculate doppler shifted laser frequencies
        VT laserFreqsDoppler = laserFreqs * (1 - velocity / SpeedOfLight2_v);
        
        // Rotating frame
        // hamiltonian(0, 0) -= TwoPi_v * (~laserFreqsDoppler)(0);

        // Calculate electric field
        constexpr double twoOverEps0c = 2 / (VacuumPermittivity_v * SpeedOfLight2_v);
        std::complex<double> electricField = 0;
        for (std::size_t i = 0; i < (~laserFreqsDoppler).Size(); i++)
        {
            double E0 = twoOverEps0c * std::sqrt((~laserIntensities)(i));
            electricField += E0 * std::cos(TwoPi_v * (~laserFreqsDoppler)(i) * t);
            // electricField += E0 * 0.5 * std::exp(std::complex<double>(1.0i * TwoPi_v * (~laserFreqsDoppler)(i) * t));
            // electricField += E0 * 0.5 * (1.0 + std::exp(std::complex<double>(-2.0i * TwoPi_v * (~laserFreqsDoppler)(i) * t)));
        }

        // System-Light interaction
        hamiltonian -= (electricField / ReducedPlanckConstant_v) * m_dipoleOperator;

        return hamiltonian; 
    }
    
    template<std::size_t N>
    template<typename VT>
    std::vector<TStaticDensityMatrix<N>> TStaticQSys<N>::GetTrajectoryNatural(
        const TColVector<VT>& laserFreqs,
        const TColVector<VT>& laserIntensities, 
        const TStaticDensityMatrix<N>& initial, 
        double dt, std::size_t steps)
    {
        std::vector<TStaticDensityMatrix<N>> trajectory;
        trajectory.reserve(steps + 1);
        trajectory.push_back(initial);
        
        TStaticMatrix<std::complex<double>, N, N> rho = initial;
        for (std::size_t i = 1; i <= steps; i++)
        {
            double t = i * dt;
            auto func = [&](auto x, const auto& y) { return GetDensityOpDerivative(y, laserFreqs, laserIntensities, 0.0, x); };

            
            // EulerIntegrator::Step(rho, t, dt, func);
            RK4Integrator::Step(rho, t, dt, func);
            
            trajectory.emplace_back(m_levelNames, rho);
        }
        
        return trajectory;
    }

    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TStaticQSys<N>::GetDensityOpDerivative(
        const TStaticMatrix<std::complex<double>, N, N>& rho,
        const TColVector<VT>& laserFreqs,
        const TColVector<VT>& laserIntensities, 
        double velocity, double t) const
    {
        // von Neumann term
        auto h = GetHamiltonian(laserFreqs, laserIntensities, 0, t);
        TStaticMatrix<std::complex<double>, N, N> rhoPrime = -1.0i * (h * rho - rho * h);

        // add lindblad dissipation term
        for (const auto& decay: m_decays)
        {
            std::complex<double> popDecayRate = decay.rate * rho(decay.from, decay.from);
            rhoPrime(decay.from, decay.from) -= popDecayRate;
            rhoPrime(decay.to, decay.to) += popDecayRate;
            rhoPrime(decay.from, decay.to) -= 0.5 * decay.rate * rho(decay.from, decay.to);
            rhoPrime(decay.to, decay.from) -= 0.5 * decay.rate * rho(decay.to, decay.from);
        }

        return rhoPrime;
    }
}

#endif
