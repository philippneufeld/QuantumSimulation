// Philipp Neufeld, 2021

#ifndef QSim_QSim_StaticQSysSS_H_
#define QSim_QSim_StaticQSysSS_H_

#include <cstdint>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <cassert>

#include "Math/Matrix.h"
#include "Doppler.h"

namespace QSim
{
    constexpr static double SpeedOfLight_v = 2.99792458e8;

    //
    // Density matrix class
    //

    template<std::size_t N>
    class TStaticDensityMatrix
    {
    public:
        TStaticDensityMatrix(std::map<std::string, std::size_t> levelNames, 
            const TStaticMatrix<std::complex<double>, N, N>& densityMatrix)
            : m_levelNames(levelNames), m_densityMatrix(densityMatrix) {}

        TStaticDensityMatrix(const TStaticDensityMatrix&) = default;
        TStaticDensityMatrix& operator=(const TStaticDensityMatrix&) = default;
        
        double GetPopulation(const std::string& lvl) const;
        double GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const;

        const TStaticMatrix<std::complex<double>, N, N>& GetMatrix() const { return m_densityMatrix; }

    private:
        std::map<std::string, std::size_t> m_levelNames;
        TStaticMatrix<std::complex<double>, N, N> m_densityMatrix;
    };

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetPopulation(const std::string& lvl) const 
    { 
        std::size_t idx = m_levelNames.at(lvl);
        return std::real(m_densityMatrix(idx, idx)); 
    }

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const 
    { 
        std::size_t idx1 = m_levelNames.at(lvl1);
        std::size_t idx2 = m_levelNames.at(lvl2);
        return std::imag(m_densityMatrix(idx1, idx2)); 
    }


    //
    // N-level quantum system steady state solver
    //

    template<std::size_t N>
    class TStaticQSysSS
    {
        template<typename InputIt>
        using EnableIfLvlIt_t = std::enable_if_t<
            std::is_same<std::string, std::decay_t<decltype(std::declval<InputIt>()->first)>>::value &&
            std::is_same<double, std::decay_t<decltype(std::declval<InputIt>()->second)>>::value>;

    public:
        // constructors
        TStaticQSysSS(const std::map<std::string, double>& levels, double mass)
            : TStaticQSysSS(levels.begin(), mass) { assert(levels.size() == N); }
        template<typename InputIt, typename=EnableIfLvlIt_t<InputIt>>
        TStaticQSysSS(InputIt levelIterator, double mass);

        // copy operations
        TStaticQSysSS(const TStaticQSysSS&) = default;
        TStaticQSysSS& operator=(const TStaticQSysSS&) = default;

        // get index by name
        std::size_t GetLevelIndexByName(const std::string& name) { return m_levelNames.at(name); }

        // functions to change the system properties 
        bool AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi);
        bool AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi);

        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetHamiltonian(const TColVector<VT>& detunings, double velocity) const;

        // thermal environment
        void SetMass(double mass) { m_doppler.SetMass(mass); }
        void SetTemperature(double temp) { m_doppler.SetTemperature(temp); }

        double GetMass() const { return m_doppler.GetMass(); }
        double GetTemperature() const { return m_doppler.GetTemperature(); }

        // Steady state
        template<typename VT>
        TStaticDensityMatrix<N> GetSteadyStateNatural(const TColVector<VT>& detunings, double velocity) const;
        template<typename VT>
        TStaticDensityMatrix<N> GetSteadyState(const TColVector<VT>& detunings) const;
        
    private:
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
            std::set<std::size_t>& visitedLevels, std::size_t transFrom);

    private:
        // Map from the level name to their index
        std::map<std::string, std::size_t> m_levelNames;

        // Properties of the system
        TStaticColVector<double, N> m_levels;
        std::vector<std::tuple<std::size_t, std::size_t, double>> m_transitions;
        std::vector<std::tuple<std::size_t, std::size_t, double>> m_decays;

        // thermal environment
        TDopplerIntegrator<double> m_doppler;

        // auxilliary variables
        TDynamicColVector<double> m_transResonances;
        TDynamicColVector<double> m_dopplerFactors;
        TDynamicMatrix<double> m_photonBasis;
        TStaticMatrix<std::complex<double>, N, N> m_hamiltonianNoLight;
    };

    template<std::size_t N>
    template<typename InputIt, typename>
    TStaticQSysSS<N>::TStaticQSysSS(InputIt levelIterator, double mass)
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
    bool TStaticQSysSS<N>::AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi)
    {
        // check that both levels exist in the system
        if (m_levelNames.find(lvl1) == m_levelNames.end() || m_levelNames.find(lvl2) == m_levelNames.end())
            return false;

        m_transitions.emplace_back(m_levelNames[lvl1], m_levelNames[lvl2], rabi);
        
        // try to preapre the calculation
        if (!PrepareCalculation())
        {
            m_transitions.pop_back();
            PrepareCalculation();
            return false;
        }

        return true;
    }
    
    template<std::size_t N>
    bool TStaticQSysSS<N>::AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi)
    {
        // check that both levels exist in the system
        if (m_levelNames.find(lvlFrom) == m_levelNames.end() || m_levelNames.find(lvlTo) == m_levelNames.end())
            return false;

        m_decays.emplace_back(m_levelNames[lvlFrom], m_levelNames[lvlTo], rabi);
        return true;
    }

    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TStaticQSysSS<N>::GetHamiltonian(const TColVector<VT>& detunings, double velocity) const
    {
        assert((~detunings).Rows() == m_transResonances.Rows());
        auto hamiltonian = m_hamiltonianNoLight;

        // Calculate doppler shifted laser frequencies
        auto laserFreqs = m_transResonances + detunings;
        for (std::size_t i = 0; i < (~laserFreqs).Rows(); i++)
            laserFreqs(i, 0) *= (1 - m_dopplerFactors(i, 0) * velocity);

        // Add light term to the hamiltonian
        TStaticColVector<double, N> photonTerms;
        MatrixMul(photonTerms, m_photonBasis, laserFreqs);
        for (std::size_t i = 0; i < (~hamiltonian).Rows(); i++)
            hamiltonian(i, i) += photonTerms(i, 0); 

        return hamiltonian; 
    }
    
    template<std::size_t N>
    template<typename VT>
    TStaticDensityMatrix<N> TStaticQSysSS<N>::GetSteadyStateNatural(const TColVector<VT>& detunings, double velocity) const
    {
        const TStaticMatrix<std::complex<double>, N, N>& h = GetHamiltonian(detunings, velocity);
        TStaticMatrix<std::complex<double>, N*N + 1, N*N> A;

        // von Neumann part of the evolution operator
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                for (std::size_t k = 0; k < N; k++)
                {
                    A(i*N+k, j*N+k) -= std::complex<double>(0, 1) * h(i, j);
                    A(k*N+i, k*N+j) += std::complex<double>(0, 1) * h(j, i);
                }
            }
        }
        
        // Lindblad part of the evolution operator
        for(const auto& decay: m_decays)
        {
            std::size_t i = std::get<0>(decay);
            std::size_t f = std::get<1>(decay);
            auto rate = std::get<2>(decay);
            A(f*N+f, i*N+i) += rate;
            for (std::size_t j = 0; j < N; j++)
            {
                A(j*N+i, j*N+i) -= 0.5*rate;
                A(i*N+j, i*N+j) -= 0.5*rate;
            }
        }

        // Normalization condition
        for (std::size_t i = 0; i < N; i++)
            A(N*N, i*N+i) += 1.0;
        
        TStaticMatrix<std::complex<double>, N*N + 1, 1> b;
        b(N*N, 0) = 1.0;
        
        auto A_adj = Adjoint(A);
        auto x = LinearSolve(A_adj*A, A_adj*b);

        TStaticMatrix<std::complex<double>, N, N> ss;
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
                ss(i, j) = x(i*N+j, 0);
        }
        
        return TStaticDensityMatrix<N>(m_levelNames, ss);
    }

    template<std::size_t N>
    template<typename VT>
    TStaticDensityMatrix<N> TStaticQSysSS<N>::GetSteadyState(const TColVector<VT>& detunings) const
    {
        auto ss = m_doppler.Integrate([&](double vel){ return this->GetSteadyStateNatural(detunings, vel).GetMatrix(); });
        return TStaticDensityMatrix<N>(m_levelNames, ss);
    }

    template<std::size_t N>
    bool TStaticQSysSS<N>::PrepareCalculation()
    {
        // calculate transition splittings
        std::size_t transCnt = m_transitions.size();
        m_transResonances.Resize(transCnt);
        for (std::size_t i = 0; i < transCnt; i++)
        {
            auto lvl1 = std::get<0>(m_transitions[i]);
            auto lvl2 = std::get<1>(m_transitions[i]);
            m_transResonances(i, 0) = std::abs(m_levels[lvl1] - m_levels[lvl2]);
        }

        // calulate doppler factors
        m_dopplerFactors.Resize(transCnt);
        for (std::size_t i = 0; i < transCnt; i++)
            m_dopplerFactors(i, 0) = 1.0 / SpeedOfLight_v;

        // create photon basis matrix
        m_photonBasis.Resize(N, transCnt);
        m_photonBasis.SetZero();
        std::vector<std::size_t> trans_path;
        std::set<std::size_t> visited_levels;
        trans_path.reserve(m_transitions.size());

        while (visited_levels.size() < N)
        {
            std::size_t head = 0;
            for(; visited_levels.count(head) != 0; head++);

            if (!PreparePhotonBasis(trans_path, visited_levels, head))
            {
                m_transResonances.Resize(0);
                m_dopplerFactors.Resize(0);
                m_photonBasis.Resize(0, 0);
                return false;
            }
        }

        // Atom hamiltonian
        m_hamiltonianNoLight.SetZero();
        for (std::size_t i = 0; i < N; i++)
            m_hamiltonianNoLight(i, i) += m_levels[i];
        
        // Interaction hamiltonian
        for (const auto& trans: m_transitions)
        {
            std::size_t l1Idx = std::get<0>(trans);
            std::size_t l2Idx = std::get<1>(trans);
            auto rabi = std::get<2>(trans);
            m_hamiltonianNoLight(l1Idx, l2Idx) += 0.5*rabi;
            m_hamiltonianNoLight(l2Idx, l1Idx) += 0.5*rabi;
        }

        return true;
    }

    template<std::size_t N>
    bool TStaticQSysSS<N>::PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
        std::set<std::size_t>& visitedLvls, std::size_t transFrom)
    {
        // cehck for circular transition path
        if (visitedLvls.count(transFrom) != 0)
            return false;
        else
            visitedLvls.insert(transFrom);

        // Iterate over all transitions
        for (std::size_t tIdx = 0; tIdx < m_transitions.size(); tIdx++)
        {
            // skip if the transition was already in the current transition path
            if (std::find(trans_path.begin(), trans_path.end(), tIdx) != trans_path.end())
                continue;

            // skip if currLevel is not involved in the current transition
            // otherwise check to which level the transition leads
            std::size_t transTo;
            if (std::get<0>(m_transitions[tIdx]) == transFrom)
                transTo = std::get<1>(m_transitions[tIdx]);
            else if (std::get<1>(m_transitions[tIdx]) == transFrom)
                transTo = std::get<0>(m_transitions[tIdx]);
            else
                continue;
            
            // check whether a photon is absorbed or emitted
            double relPhotonNumber = 1;
            if (m_levels[transTo] > m_levels[transFrom])
                relPhotonNumber = -1;

            // Add(Subtract) one photon to(from) the "from" state to obtain the "to" state
            for (std::size_t j=0; j<m_photonBasis.Cols(); j++)
                m_photonBasis(transTo, j) = m_photonBasis(transFrom, j);
            m_photonBasis(transTo, tIdx) += relPhotonNumber;

            // continue recursively
            trans_path.push_back(tIdx);
            auto success = PreparePhotonBasis(trans_path, visitedLvls, transTo);
            trans_path.pop_back();

            if (!success)
                return false;
        }

        return true;
    }
}

#endif
