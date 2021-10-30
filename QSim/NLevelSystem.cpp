// Philipp Neufeld, 2021

#include "NLevelSystem2.h"

#include <algorithm>
#include <set>
#include <cassert>

namespace QSim
{
    NLevelSystem::NLevelSystem()
    {

    }   

    NLevelSystem::~NLevelSystem()
    {

    } 

    void NLevelSystem::SetLevel(const std::string& name, double level)
    {
        // insert (or change) level
        m_levels[name] = level;
    }

    bool NLevelSystem::AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi)
    {
        // check that both levels exist in the system
        if (m_levels.find(lvl1) == m_levels.end() || m_levels.find(lvl2) == m_levels.end())
            return false;

        m_transitions.emplace_back(lvl1, lvl2, rabi);
        
        // try to preapre the calculation
        if (!PrepareCalculation())
        {
            m_transitions.pop_back();
            PrepareCalculation();
            return false;
        }

        return true;
    }
    
    bool NLevelSystem::AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi)
    {
        // check that both levels exist in the system
        if (m_levels.find(lvlFrom) == m_levels.end() || m_levels.find(lvlTo) == m_levels.end())
            return false;

        m_decays.emplace_back(lvlFrom, lvlTo, rabi);

        // try to preapre the calculation
        if (!PrepareCalculation())
        {
            m_decays.pop_back();
            PrepareCalculation();
            return false;
        }

        return true;
    }

    TDynamicMatrix<double> NLevelSystem::GetHamiltonian(const TDynamicMatrix<double>& detunings, double velocity)
    {
        assert((~detunings).Rows() == m_transitionSplittings.Rows());

        // Add light term to the hamiltonian
        TDynamicMatrix<double> hamiltonian = m_hamiltonianNoLight;
        auto laserFreqs = detunings + m_transitionSplittings;
        for (std::size_t i = 0; i < (~laserFreqs).Rows(); i++)
            laserFreqs(i, 0) *= (1 - m_dopplerFactors(i, 0) * velocity);
        auto diag = m_photonBasis * laserFreqs;
        for (std::size_t i = 0; i < (~hamiltonian).Rows(); i++)
            hamiltonian(i, i) += diag(i, 0); 

        return hamiltonian; 
    }

    TDynamicMatrix<std::complex<double>> NLevelSystem::GetSteadyState(const TDynamicMatrix<double>& detunings, double velocity)
    {

        auto h = GetHamiltonian(detunings, velocity);
        std::size_t N = h.Rows();
        TDynamicMatrix<std::complex<double>> A(N*N + 1, N*N);

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
        for(const auto& decay: m_indexedDecays)
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
        
        TDynamicMatrix<std::complex<double>> b(N*N+1, 1);
        b(N*N, 0) = 1.0;
        
        auto A_adj = Adjoint(A);
        auto x = LinearSolve(A_adj*A, A_adj*b);

        TDynamicMatrix<std::complex<double>> ss(N, N);
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
                ss(i, j) = x(i*N+j, 0);
        }
        
        return ss;
    }
    
    double NLevelSystem::GetAbsorptionCoeff(const TDynamicMatrix<double>& detunings, double velocity, const std::string& lvl1, const std::string lvl2)
    {
        std::map<std::string, std::size_t> indexMap;
        for (std::size_t i = 0; i < m_usedLevels.size(); i++)
            indexMap[m_usedLevels[i]] = i;

        return std::imag(GetSteadyState(detunings, velocity)(indexMap[lvl1], indexMap[lvl2]));
    }

    bool NLevelSystem::PrepareCalculation()
    {
        // find used levels
        std::set<std::string> usedLevels;
        for (const auto& transition : m_transitions)
            usedLevels.insert({std::get<0>(transition), std::get<1>(transition)});
        for (const auto& decay : m_decays)
            usedLevels.insert({std::get<0>(decay), std::get<1>(decay)});
        m_usedLevels.assign(usedLevels.begin(), usedLevels.end());

        std::map<std::string, std::size_t> usedLevelsIndices;
        for (std::size_t i = 0; i < m_usedLevels.size(); i++)
            usedLevelsIndices[m_usedLevels[i]] = i;
        
        // 
        m_indexedDecays.resize(m_decays.size());
        for (std::size_t i = 0; i < m_indexedDecays.size(); i++)
            m_indexedDecays[i] = std::make_tuple(usedLevelsIndices[std::get<0>(m_decays[i])],
                usedLevelsIndices[std::get<1>(m_decays[i])], std::get<2>(m_decays[i]));
        

        // get dimensions
        std::size_t dim = m_usedLevels.size();
        std::size_t transCnt = m_transitions.size();

        // calculate transition splittings
        m_transitionSplittings.Resize(transCnt, 1);
        for (std::size_t i = 0; i < transCnt; i++)
        {
            auto lvl1 = m_levels[std::get<0>(m_transitions[i])];
            auto lvl2 = m_levels[std::get<1>(m_transitions[i])];
            m_transitionSplittings(i, 0) = std::abs(lvl1 - lvl2);
        }

        // calulate doppler factors
        m_dopplerFactors.Resize(transCnt, 1);
        for (std::size_t i = 0; i < transCnt; i++)
            m_dopplerFactors(i, 0) = 1.0 / SpeedOfLight2_v;

        // allocate photon basis matrix of the right size and initialize it with zero
        m_photonBasis.Resize(dim, transCnt);
        m_photonBasis.SetZero();

        // Perform the photon basis(deallocate matrix on failure)
        std::vector<std::size_t> trans_path;
        trans_path.reserve(m_transitions.size());

        if (!PreparePhotonBasis(trans_path, m_usedLevels[0]))
        {
            m_usedLevels.clear();
            m_transitionSplittings.Resize(0, 0);
            m_dopplerFactors.Resize(0, 0);
            m_photonBasis.Resize(0, 0);
            m_hamiltonianNoLight.Resize(0, 0);
            return false;
        }

        // Prepare hamiltonian (Atom and interaction part)
        m_hamiltonianNoLight.Resize(dim, dim);
        m_hamiltonianNoLight.SetZero();

        // Atom hamiltonian
        for (std::size_t i = 0; i < dim; i++)
            m_hamiltonianNoLight(i, i) += m_levels[m_usedLevels[i]];
        
        // Interaction hamiltonian
        for (const auto& trans: m_transitions)
        {
            std::size_t l1Idx = usedLevelsIndices[std::get<0>(trans)];
            std::size_t l2Idx = usedLevelsIndices[std::get<1>(trans)];
            auto rabi = std::get<2>(trans);
            m_hamiltonianNoLight(l1Idx, l2Idx) += 0.5*rabi;
            m_hamiltonianNoLight(l2Idx, l1Idx) += 0.5*rabi;
        }

        return true;
    }

    bool NLevelSystem::PreparePhotonBasis(std::vector<std::size_t>& trans_path, const std::string& transFrom)
    {
        for (std::size_t tIdx = 0; tIdx < m_transitions.size(); tIdx++)
        {
            // skip if the transition was already in the current transition path
            if (std::find(trans_path.begin(), trans_path.end(), tIdx) != trans_path.end())
                continue;

            // skip if currLevel is not involved in the current transition
            // otherwise check to which level the transition leads
            std::string transTo;
            if (std::get<0>(m_transitions[tIdx]) == transFrom)
                transTo = std::get<1>(m_transitions[tIdx]);
            else if (std::get<1>(m_transitions[tIdx]) == transFrom)
                transTo = std::get<0>(m_transitions[tIdx]);
            else
                continue;
            
            // calculate the indices in the basis matrix for the levels involved in the transition
            const auto& ul = m_usedLevels;
            std::size_t fromIdx = std::find(ul.begin(), ul.end(), transFrom) - ul.begin();
            std::size_t toIdx = std::find(ul.begin(), ul.end(), transTo) - ul.begin();

            // check whether a photon is absorbed or emitted
            double relPhotonNumber = 1;
            if (m_levels[transTo] > m_levels[transFrom])
                relPhotonNumber = -1;

            // Add(Subtract) one photon to(from) the "from" state to obtain the "to" state
            for (std::size_t j=0; j<m_photonBasis.Cols(); j++)
            {
                // if the row is not zero, the level has already been visited 
                // by another transition (circular transition path)
                if (m_photonBasis(toIdx, j) != 0)
                    return false;
                m_photonBasis(toIdx, j) = m_photonBasis(fromIdx, j);
            }
            m_photonBasis(toIdx, tIdx) += relPhotonNumber;

            // continue recursively
            trans_path.push_back(tIdx);
            auto success = PreparePhotonBasis(trans_path, transTo);
            trans_path.pop_back();

            if (!success)
                return false;
        }

        return true;
    }

} 
