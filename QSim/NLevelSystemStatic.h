// Philipp Neufeld, 2021

#ifndef QSIM_NLevelSystemStatic_H_
#define QSIM_NLevelSystemStatic_H_

#include <cstdint>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <cassert>

#include "Matrix.h"

namespace QSim
{
    constexpr static double SpeedOfLight_v = 2.99792458e8;

    template<std::size_t N>
    class TStaticNLevelSystem
    {
        template<typename InputIt>
        using EnableIfLvlIt_t = std::enable_if_t<
            std::is_same<std::string, std::decay_t<decltype(std::declval<InputIt>()->first)>>::value &&
            std::is_same<double, std::decay_t<decltype(std::declval<InputIt>()->second)>>::value>;

    public:
        // constructors
        TStaticNLevelSystem(const std::map<std::string, double>& levels)
            : TStaticNLevelSystem(levels.begin()) { assert(levels.size() == N); }
        template<typename InputIt, typename=EnableIfLvlIt_t<InputIt>>
        TStaticNLevelSystem(InputIt levelIterator);

        // copy operations
        TStaticNLevelSystem(const TStaticNLevelSystem&) = default;
        TStaticNLevelSystem& operator=(const TStaticNLevelSystem&) = default;

        // functions to change the system properties 
        bool AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi);
        bool AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi);

        template<typename VT>
        TStaticMatrix<double, N, N> GetHamiltonian(const TMatrix<VT>& detunings, double velocity) const;
        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetSteadyState(const TMatrix<VT>& detunings, double velocity) const;
        template<typename VT>
        double GetAbsorptionCoeff(const TMatrix<VT>& detunings, double velocity, const std::string& lvl1, const std::string lvl2) const;

    private:
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
            std::set<std::size_t>& visitedLevels, std::size_t transFrom);

    private:
        // Map from the level name to their index
        std::map<std::string, std::size_t> m_levelNames;

        // Properties of the system
        std::array<double, N> m_levels;
        std::vector<std::tuple<std::size_t, std::size_t, double>> m_transitions;
        std::vector<std::tuple<std::size_t, std::size_t, double>> m_decays;

        // auxilliary variables
        TDynamicMatrix<double> m_transResonances;
        TDynamicMatrix<double> m_dopplerFactors;
        TDynamicMatrix<double> m_photonBasis;
        TStaticMatrix<double, N, N> m_hamiltonianNoLight;
    };

    template<std::size_t N>
    template<typename InputIt, typename>
    TStaticNLevelSystem<N>::TStaticNLevelSystem(InputIt levelIterator)
    {
        // levelIterator is a pair containing (name, level)
        for (std::size_t i = 0; i < N; i++, levelIterator++)
        {
            m_levelNames[levelIterator->first] = i;
            m_levels[i] = levelIterator->second;
        }
    }

    template<std::size_t N>
    bool TStaticNLevelSystem<N>::AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi)
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
    bool TStaticNLevelSystem<N>::AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi)
    {
        // check that both levels exist in the system
        if (m_levelNames.find(lvlFrom) == m_levelNames.end() || m_levelNames.find(lvlTo) == m_levelNames.end())
            return false;

        m_decays.emplace_back(m_levelNames[lvlFrom], m_levelNames[lvlTo], rabi);
        return true;
    }

    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<double, N, N> TStaticNLevelSystem<N>::GetHamiltonian(const TMatrix<VT>& detunings, double velocity) const
    {
        assert((~detunings).Rows() == m_transResonances.Rows());
        TStaticMatrix<double, N, N> hamiltonian = m_hamiltonianNoLight;

        // Calculate doppler shifted laser frequencies
        auto laserFreqs = m_transResonances;
        MatrixAdd(laserFreqs, laserFreqs, detunings);
        for (std::size_t i = 0; i < (~laserFreqs).Rows(); i++)
            laserFreqs(i, 0) *= (1 - m_dopplerFactors(i, 0) * velocity);

        // Add light term to the hamiltonian
        TStaticMatrix<double, N, 1> photonTerms;
        MatrixMul(photonTerms, m_photonBasis, laserFreqs);
        for (std::size_t i = 0; i < (~hamiltonian).Rows(); i++)
            hamiltonian(i, i) += photonTerms(i, 0); 

        return hamiltonian; 
    }
    
    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TStaticNLevelSystem<N>::GetSteadyState(const TMatrix<VT>& detunings, double velocity) const
    {
        TStaticMatrix<double, N, N> h = GetHamiltonian(detunings, velocity);
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
        
        return ss;
    }

    template<std::size_t N>
    template<typename VT>
    double TStaticNLevelSystem<N>::GetAbsorptionCoeff(const TMatrix<VT>& detunings, 
        double velocity, const std::string& lvl1, const std::string lvl2) const
    {
        return std::imag(GetSteadyState(detunings, velocity)(m_levelNames.at(lvl1), m_levelNames.at(lvl2)));
    }

    template<std::size_t N>
    bool TStaticNLevelSystem<N>::PrepareCalculation()
    {
        // calculate transition splittings
        std::size_t transCnt = m_transitions.size();
        m_transResonances.Resize(transCnt, 1);
        for (std::size_t i = 0; i < transCnt; i++)
        {
            auto lvl1 = std::get<0>(m_transitions[i]);
            auto lvl2 = std::get<1>(m_transitions[i]);
            m_transResonances(i, 0) = std::abs(m_levels[lvl1] - m_levels[lvl2]);
        }

        // calulate doppler factors
        m_dopplerFactors.Resize(transCnt, 1);
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
                m_transResonances.Resize(0, 0);
                m_dopplerFactors.Resize(0, 0);
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
    bool TStaticNLevelSystem<N>::PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
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
