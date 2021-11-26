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
#include "Math/Integrator.h"
#include "Doppler.h"
#include "DensityMatrix.h"

#include <iostream>

namespace QSim
{
    constexpr static double Pi_v = 3.14159265358979323846;
    constexpr static double TwoPi_v = 2 * Pi_v;
    constexpr static double SpeedOfLight2_v = 2.99792458e8;
    constexpr static double VacuumPermittivity_v = 8.8541878128e-12;
    constexpr static double PlanckConstant_v = 6.62607004e-34;
    constexpr static double ReducedPlanckConstant_v = 1.054571817e-34;
    constexpr static double AtomicMassUnit_v = 1.66053906660e-27;

    //
    // N-level quantum system solver
    //

    namespace Internal
    {
        struct StaticQSysDecayDesc
        {
            std::size_t from;
            std::size_t to;
            double rate;
        };
    }

    template<std::size_t N, typename MyT>
    class TStaticQLvlSys
    {
        using IndexPair = std::pair<std::size_t, std::size_t>;
        template<typename InputIt>
        using EnableIfLvlIt_t = std::enable_if_t<
            std::is_same<std::string, std::decay_t<decltype(std::declval<InputIt>()->first)>>::value &&
            std::is_same<double, std::decay_t<decltype(std::declval<InputIt>()->second)>>::value>;
    public:
        // constructors
        TStaticQLvlSys();
        TStaticQLvlSys(const std::array<double, N>& levels);
        TStaticQLvlSys(const std::array<std::string, N>& lvlNames);
        TStaticQLvlSys(const std::array<std::string, N>& lvlNames, const std::array<double, N>& levels);

        // copy operations
        TStaticQLvlSys(const TStaticQLvlSys&) = default;
        TStaticQLvlSys& operator=(const TStaticQLvlSys&) = default;

        // CRTP operators
        MyT& operator~() { return static_cast<MyT&>(*this); }
        const MyT& operator~() const { return static_cast<const MyT&>(*this); }

        // level names
        std::string GetLevelNameByIndex(std::size_t idx) const;
        std::size_t GetLevelIndexByName(const std::string& name) const;
        bool SetLevelName(std::size_t idx, const std::string& newName);

        // levels
        double GetLevelByName(const std::string& name) const;
        bool SetLevelByName(const std::string& name, double level);
        double GetLevel(std::size_t idx) const;
        bool SetLevel(std::size_t idx, double level);

        // decay rates due to spontaneous emission
        bool SetDecay(std::size_t from, std::size_t to, double rate);
        bool SetDecayByName(const std::string& from, std::string& to, double rate);
        double GetDecay(std::size_t from, std::size_t to) const;
        double GetDecayByName(const std::string& from, std::string& to) const;

        // Transition dipole operator
        bool SetTransitionDipole(std::size_t from, std::size_t to, double dip);
        bool SetTransitionDipoleByName(const std::string& from, std::string& to, double rate);
        double GetTransitionDipole(std::size_t from, std::size_t to) const;
        double GetTransitionDipoleByName(const std::string& from, std::string& to) const;

        // thermal environment and properties needed for the doppler integration
        void SetMass(double mass) { m_doppler.SetMass(mass); }
        void SetTemperature(double temp) { m_doppler.SetTemperature(temp); }
        double GetMass() const { return m_doppler.GetMass(); }
        double GetTemperature() const { return m_doppler.GetTemperature(); }

        // create stecific density matrices
        TStaticDensityMatrix<N> CreateGroundState() const;
        TStaticDensityMatrix<N> CreateThermalState() const;

    private:
        TStaticDensityMatrix<N> CreateThermalStateHelper(double temperature) const;

        std::array<double, N> GenerateDefaultLevels() const;
        std::array<std::string, N> GenerateDefaultLevelNames() const;

    protected:
        // Properties of the system
        std::array<std::string, N> m_levelNames;
        TStaticColVector<double, N> m_levels;
        std::map<IndexPair, double> m_decays;
        TStaticMatrix<std::complex<double>, N, N> m_dipoleOperator;

        // thermal environment
        TDopplerIntegrator<double> m_doppler;
    };

    template<std::size_t N, typename MyT>
    TStaticQLvlSys<N, MyT>::TStaticQLvlSys()
        : TStaticQLvlSys(GenerateDefaultLevelNames(), GenerateDefaultLevels()) { }

    template<std::size_t N, typename MyT>
    TStaticQLvlSys<N, MyT>::TStaticQLvlSys(const std::array<double, N>& levels)
        : TStaticQLvlSys(GenerateDefaultLevelNames(), levels) { }

    template<std::size_t N, typename MyT>
    TStaticQLvlSys<N, MyT>::TStaticQLvlSys(const std::array<std::string, N>& lvlNames)
        : TStaticQLvlSys(lvlNames, GenerateDefaultLevels()) { }

    template<std::size_t N, typename MyT>
    TStaticQLvlSys<N, MyT>::TStaticQLvlSys(const std::array<std::string, N>& lvlNames, 
        const std::array<double, N>& levels)
        : m_levelNames(lvlNames), m_doppler(AtomicMassUnit_v, 300.0)
    {
        for (size_t i = 0; i < N; i++)
            m_levels(i) = levels[i];
    }

    template<std::size_t N, typename MyT>
    std::string TStaticQLvlSys<N, MyT>::GetLevelNameByIndex(std::size_t idx) const
    {
        return idx < N ? m_levelNames[idx] : std::string();
    }

    template<std::size_t N, typename MyT>
    std::size_t TStaticQLvlSys<N, MyT>::GetLevelIndexByName(const std::string& name) const
    {
        auto it = std::find(m_levelNames.begin(), m_levelNames.end(), name);
        return it != m_levelNames.end() ? it - m_levelNames.begin() : -1;
    }
    
    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetLevelName(std::size_t idx, const std::string& newName)
    {
        if (GetLevelIndexByName(newName) < N)
            return false;  // level name already present
        m_levelNames[idx] = newName;
        return true;
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetLevelByName(const std::string& name) const
    {
        return GetLevel(GetLevelIndexByName(name));
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetLevelByName(const std::string& name, double level)
    {
        return SetLevel(GetLevelIndexByName(name), level);
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetLevel(std::size_t idx) const
    {
        return idx < N ? m_levels[idx] : 0.0;
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetLevel(std::size_t idx, double level)
    {
        if (idx >= N)
            return false;
        m_levels[idx] = level;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetDecay(std::size_t from, std::size_t to, double rate)
    {
        if (from >= N || to >= N)
            return false;  // index out of bound
        m_decays[std::make_pair(from, to)] = rate;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetDecayByName(const std::string& from, std::string& to, double rate)
    {
        return SetDecay(GetLevelIndexByName(from), GetLevelIndexByName(to), rate);
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetDecay(std::size_t from, std::size_t to) const
    {
        auto it = m_decays.find(std::make_pair(from, to));
        return it != m_decays.end() ? it->second : 0.0;
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetDecayByName(const std::string& from, std::string& to) const
    {
        return GetDecay(GetLevelIndexByName(from), GetLevelIndexByName(to));
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetTransitionDipole(std::size_t from, std::size_t to, double dip)
    {
        if (from >= N || to >= N)
            return false;  // index out of bound
        m_dipoleOperator(from, to) = dip;
        m_dipoleOperator(to, from) = dip;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TStaticQLvlSys<N, MyT>::SetTransitionDipoleByName(const std::string& from, std::string& to, double dip)
    {
        return SetTransitionDipole(GetLevelIndexByName(from), GetLevelIndexByName(to), dip);
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetTransitionDipole(std::size_t from, std::size_t to) const
    {
        return (from < N && to < N) ? m_dipoleOperator(from, to) : 0.0;
    }

    template<std::size_t N, typename MyT>
    double TStaticQLvlSys<N, MyT>::GetTransitionDipoleByName(const std::string& from, std::string& to) const
    {
        return GetTransitionDipole(GetLevelIndexByName(from), GetLevelIndexByName(to));
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TStaticQLvlSys<N, MyT>::CreateGroundState() const
    { 
        return CreateThermalStateHelper(0.0);
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TStaticQLvlSys<N, MyT>::CreateThermalState() const
    { 
        return CreateThermalStateHelper(GetTemperature()); 
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TStaticQLvlSys<N, MyT>::CreateThermalStateHelper(double temperature) const
    { 
        TStaticDensityMatrix<N> res;

        TStaticColVector<double, N> lvlOffset;
        std::size_t gsIdx = std::min_element(m_levels.begin(), m_levels.end()) - m_levels.begin();
        for (size_t i = 0; i < N; i++)
            lvlOffset(i) = m_levels(i) - m_levels(gsIdx);

        TStaticColVector<double, N> thermalWeights;
        if (temperature <= 0.0)
        {
            for (size_t i = 0; i < N; i++)
            {
                if (lvlOffset(i) == 0.0)
                    thermalWeights(i) = 1.0;
            }
        }
        else
        {
            double tmp = PlanckConstant_v / (BoltzmannConstant_v * temperature);
            for (size_t i = 0; i < N; i++)
                thermalWeights(i) = std::exp(-tmp * lvlOffset(i));
        }

        double norm = 0.0;
        for (size_t i = 0; i < N; i++)
            norm += thermalWeights(i);
        thermalWeights /= norm;

        for (size_t i = 0; i < N; i++)
            res(i, i) = thermalWeights(i);
        return res;
    }

    template<std::size_t N, typename MyT>
    std::array<double, N> TStaticQLvlSys<N, MyT>::GenerateDefaultLevels() const
    {
        std::array<double, N> lvls;
        std::fill(lvls.begin(), lvls.end(), 0.0);
        return lvls;
    }

    template<std::size_t N, typename MyT>
    std::array<std::string, N> TStaticQLvlSys<N, MyT>::GenerateDefaultLevelNames() const
    {
        std::array<std::string, N> names;
        for (size_t i = 0; i < N; i++)
            names[i] = std::to_string(i);
        return names;
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
        std::pair<TDynamicColVector<double>, std::vector<TStaticDensityMatrix<N>>> GetTrajectoryNatural(
            const TColVector<VT>& laserFreqs,
            const TColVector<VT>& laserIntensities, 
            const TStaticDensityMatrix<N>& initial, 
            double dt, double tmax);

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
        hamiltonian(0, 0) -= TwoPi_v * (~laserFreqsDoppler)(0);

        // Calculate electric field
        constexpr double twoOverEps0c = 2 / (VacuumPermittivity_v * SpeedOfLight2_v);
        std::complex<double> electricField = 0;
        for (std::size_t i = 0; i < (~laserFreqsDoppler).Size(); i++)
        {
            double E0 = twoOverEps0c * std::sqrt((~laserIntensities)(i));
            // electricField += E0 * std::cos(TwoPi_v * (~laserFreqsDoppler)(i) * t);
            // electricField += E0 * 0.5 * std::exp(std::complex<double>(1.0i * TwoPi_v * (~laserFreqsDoppler)(i) * t));
            electricField += E0 * 0.5 * (1.0 + std::exp(std::complex<double>(-2.0i * TwoPi_v * (~laserFreqsDoppler)(i) * t)));
        }

        // System-Light interaction
        hamiltonian -= (electricField / ReducedPlanckConstant_v) * m_dipoleOperator;

        return hamiltonian; 
    }
    
    template<std::size_t N>
    template<typename VT>
    std::pair<TDynamicColVector<double>, std::vector<TStaticDensityMatrix<N>>> TStaticQSys<N>::GetTrajectoryNatural(
        const TColVector<VT>& laserFreqs,
        const TColVector<VT>& laserIntensities, 
        const TStaticDensityMatrix<N>& initial, 
        double dt, double tmax)
    {
        std::size_t steps = static_cast<std::size_t>(std::ceil(tmax / dt));

        std::vector<TStaticDensityMatrix<N>> trajectory;
        trajectory.reserve(steps + 1);
        trajectory.push_back(initial);
        
        using YType = TStaticMatrix<std::complex<double>, N, N>;
        RK4Integrator<double, YType> integrator;
        YType rho = initial;

        for (std::size_t i = 1; i <= steps; i++)
        {
            double t = i * dt;
            auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(y, laserFreqs, laserIntensities, 0.0, x); };

            rho += integrator.Step(rho, t, dt, func);         
            trajectory.emplace_back(m_levelNames, rho);
        }
        
        return {QSim::CreateLinspaceCol(0.0, steps*dt, steps + 1), trajectory};
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
