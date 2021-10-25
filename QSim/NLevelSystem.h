// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystem_H_
#define QSim_NLevelSystem_H_

#include <cstdint>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <complex>

#include "Matrix.h"
#include "Transition.h"
#include "TransitionTree.h"
#include "Decay.h"


namespace QSim
{

    constexpr static double SpeedOfLight_v = 2.99792458e8;

    template<std::size_t N, typename Ty = double>
    class TNLevelSystem
    {
    public:
        // constructors
        TNLevelSystem(const Ty (&levels)[N]) : TNLevelSystem(static_cast<const Ty*>(levels)) { } 
        TNLevelSystem(const std::initializer_list<Ty> levels) : TNLevelSystem(levels.begin()) { assert(levels.size() == N); }
        template<typename VT>
        TNLevelSystem(const TMatrix<VT>& levels) : m_levels(levels) {}
        TNLevelSystem(const Ty* levels); 
        ~TNLevelSystem();

        // copy operations
        TNLevelSystem(const TNLevelSystem& rhs) = default;
        TNLevelSystem& operator=(const TNLevelSystem& rhs) = default;

        // Transition management functions
        bool AddTransition(const TTransition<Ty>& trans);
        const TTransition<Ty>& GetTransition(std::size_t idx) const { return m_transitions[idx]; }
        std::size_t GetTransitionCount() const { return m_transitions.size(); }
        std::size_t GetTransitionIdx(const TTransition<Ty>& trans) const;

        // Decays
        void AddDecay(const TDecay<Ty>& decay) { m_decays.push_back(decay); }

        // Solver
        template<typename VT>
        TStaticMatrix<std::complex<Ty>, N, N> GetSteadyState(const TMatrix<VT>& detunings, Ty velocity);
        template<typename VT>
        Ty GetAbsorptionCoeff(const TMatrix<VT>& detunings, std::size_t lIdx1, std::size_t lIdx2);

        const TDynamicMatrix<Ty>& GetPhotonBasis() const { return m_photonBasis; }
        template<typename VT>
        TStaticMatrix<Ty, N, N> GetHamiltonian(const TMatrix<VT>& detunings, Ty velocity) const;

    private:
        TStaticMatrix<Ty, N, 1> m_levels;
        std::vector<TTransition<Ty>> m_transitions;
        std::vector<TDecay<Ty>> m_decays;
        
        TDynamicMatrix<Ty> m_photonBasis;
        TDynamicMatrix<Ty> m_transitionFreqs;
        TDynamicMatrix<Ty> m_dopplerFactors;
        TStaticMatrix<Ty, N, N> m_partialHamiltonian;
    };

    template<std::size_t N, typename Ty>
    TNLevelSystem<N, Ty>::TNLevelSystem(const Ty* levels)
    {
        for (std::size_t i = 0; i < N; i++)
            m_levels(i, 0) = levels[i];
    }

    template<std::size_t N, typename Ty>
    TNLevelSystem<N, Ty>::~TNLevelSystem()
    {

    }
    
    template<std::size_t N, typename Ty>
    bool TNLevelSystem<N, Ty>::AddTransition(const TTransition<Ty>& trans)
    {
        if (std::find(m_transitions.begin(), m_transitions.end(), trans) != m_transitions.end())
            return false;
        m_transitions.push_back(trans);

        // Recalculate transition frequencies
        m_transitionFreqs.Resize(m_transitions.size(), 1);
        m_dopplerFactors.Resize(m_transitions.size(), 1);
        for (std::size_t i = 0; i < m_transitions.size(); i++)
        {
            Ty lvl1 = m_levels(m_transitions[i].GetLevel1Index(), 0); 
            Ty lvl2 = m_levels(m_transitions[i].GetLevel2Index(), 0); 
            m_transitionFreqs(i, 0) = std::abs(lvl1 - lvl2);
            m_dopplerFactors(i, 0) = 1.0 / SpeedOfLight_v;
        }

        // Rebuild partial hamiltonian
        m_partialHamiltonian.SetZero();

        // Atom hamiltonian
        for (std::size_t i = 0; i < N; i++)
            m_partialHamiltonian(i, i) = m_levels(i, 0);
        
        // Interaction hamiltonian
        for (const auto& tr: m_transitions)
        {
            std::size_t lvlIdx1 = tr.GetLevel1Index();
            std::size_t lvlIdx2 = tr.GetLevel2Index();
            Ty rabi = tr.GetRabi();
            m_partialHamiltonian(lvlIdx1, lvlIdx2) = 0.5 * rabi;
            m_partialHamiltonian(lvlIdx2, lvlIdx1) = 0.5 * rabi;
        }
        
        // Rebuild transition tree and photon basis matrix
        std::vector<TTransitionTree<Ty>> trees;
        std::set<std::size_t> handledLevels;

        m_photonBasis.Resize(N, GetTransitionCount());
        m_photonBasis.SetZero();

        while (handledLevels.size() < N)
        {
            std::size_t head = 0;
            for (; handledLevels.find(head) != handledLevels.end(); head++);

            TTransitionTree<Ty> tree(head);
            if (!tree.BuildTree(m_transitions))
            {
                m_transitions.pop_back();
                return false;
            }
            trees.push_back(tree);
            
            auto treeIndices = tree.GetTreeLevelIndices();
            handledLevels.insert(treeIndices.begin(), treeIndices.end());

            tree.AddPhotonBasis(m_transitions, m_levels, m_photonBasis);
        }

        return true;
    }

    template<std::size_t N, typename Ty>
    std::size_t TNLevelSystem<N, Ty>::GetTransitionIdx(const TTransition<Ty>& trans) const
    {
        auto it = std::find(m_transitions.begin(), m_transitions.end(), trans);
        if (it != m_transitions.end())
            return it - m_transitions.begin();
        return -1;
    }

    template<std::size_t N, typename Ty>
    template<typename VT>
    TStaticMatrix<Ty, N, N> TNLevelSystem<N, Ty>::GetHamiltonian(const TMatrix<VT>& detunings, Ty velocity) const
    {
        assert((~detunings).Rows() == m_transitionFreqs.Rows());

        TStaticMatrix<Ty, N, N> h = m_partialHamiltonian;
        auto laserFreqs = detunings + m_transitionFreqs;
        for (std::size_t i = 0; i < (~laserFreqs).Rows(); i++)
            laserFreqs(i, 0) *= (1 - m_dopplerFactors(i, 0) * velocity);
        auto diag = m_photonBasis * laserFreqs;
        for (std::size_t i = 0; i < N; i++)
            h(i, i) += diag(i, 0); 

        return h; 
    }

    template<std::size_t N, typename Ty>
    template<typename VT>
    TStaticMatrix<std::complex<Ty>, N, N> TNLevelSystem<N, Ty>::GetSteadyState(const TMatrix<VT>& detunings, Ty velocity)
    {
        auto h = GetHamiltonian(detunings, velocity);
        TStaticMatrix<std::complex<Ty>, N*N + 1, N*N> A;

        // von Neumann part of the evolution operator
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                for (std::size_t k = 0; k < N; k++)
                {
                    A(i*N+k, j*N+k) -= std::complex<Ty>(0, 1) * h(i, j);
                    A(k*N+i, k*N+j) += std::complex<Ty>(0, 1) * h(j, i);
                }
            }
        }
        
        // Lindblad part of the evolution operator
        for(const auto& decay: m_decays)
        {
            std::size_t i = decay.GetInitialIndex();
            std::size_t f = decay.GetFinalIndex();
            Ty rate = decay.GetRate();
            A(f*N+f, i*N+i) += rate;
            for (std::size_t j = 0; j < N; j++)
            {
                A(j*N+i, j*N+i) -= 0.5*rate;
                A(i*N+j, i*N+j) -= 0.5*rate;
            }
        }

        // Normalization condition
        for (std::size_t i = 0; i < N; i++)
            A(N*N, i*N+i) += static_cast<Ty>(1);
        
        TStaticMatrix<std::complex<Ty>, N*N+1, 1> b;
        b(N*N, 0) = static_cast<Ty>(1);
        
        auto A_adj = Adjoint(A);
        auto x = LinearSolve(A_adj*A, A_adj*b);

        TStaticMatrix<std::complex<Ty>, N, N> ss;
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
                ss(i, j) = x(i*N+j, 0);
        }
        
        return ss;
    }

    template<std::size_t N, typename Ty>
    template<typename VT>
    Ty TNLevelSystem<N, Ty>::GetAbsorptionCoeff(const TMatrix<VT>& detunings, std::size_t lIdx1, std::size_t lIdx2)
    {
        return std::imag((~GetSteadyState(detunings, 0.0))(lIdx1, lIdx2));
    }
}

#endif
