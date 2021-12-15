// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystemQM_H_
#define QSim_NLevelSystemQM_H_

#include <cstdint>
#include <set>
#include <vector>

#include "../Math/Matrix.h"
#include "NLevelSystem.h"

namespace QSim
{
    //
    // Solver for a quantum mechanical Jaynes-Cummings model of a 
    // N-level system that is driven by light fields
    //

    template<std::size_t N>
    class TNLevelSystemQM : public TNLevelSystemCRTP<N, TNLevelSystemQM<N>>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
    public:
        // constructors
        TNLevelSystemQM() : MyParent() { }
        TNLevelSystemQM(const std::array<double, N>& levels) : MyParent(levels) { }
        TNLevelSystemQM(const std::array<std::string, N>& lvlNames) 
            : MyParent(lvlNames) { }
        TNLevelSystemQM(const std::array<std::string, N>& lvlNames, 
            const std::array<double, N>& levels) : MyParent(lvlNames, levels) { }

        // copy operations
        TNLevelSystemQM(const TNLevelSystemQM&) = default;
        TNLevelSystemQM& operator=(const TNLevelSystemQM&) = default;

        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetHamiltonianTI(
            const TColVector<VT>& detunings, double velocity) const;

        // Steady state
        template<typename VT>
        TStaticDensityMatrix<N> GetNaturalDensityMatrixSS(const TColVector<VT>& detunings, double velocity) const;
        template<typename VT>
        TStaticDensityMatrix<N> GetDensityMatrixSS(const TColVector<VT>& detunings) const;

    private:

        // Helper methods to prepare the hamiltonian calculation
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
            std::set<std::size_t>& visitedLevels, std::size_t transFrom);

        bool OnLaserAdded(std::size_t lvl1, std::size_t lvl2, bool counter) { return PrepareCalculation(); }
        void OnLaserRemoved() { PrepareCalculation(); }

        using HamiltonianType = TStaticMatrix<std::complex<double>, N, N>;

        template<typename VT>
        HamiltonianType GetHamiltonianAux(const TColVector<VT>& detunings, double velocity) const;
        HamiltonianType GetHamiltonianFast(const HamiltonianType& auxData, double t) const { return auxData; }

    private:
        TDynamicMatrix<double> m_photonBasis;
        TStaticMatrix<std::complex<double>, N, N> m_hamiltonianNoLight;
    };
 
    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianTI(
        const TColVector<VT>& detunings, double velocity) const
    {
        assert((~detunings).Rows() == this->GetLaserCount());
        auto hamiltonian = m_hamiltonianNoLight;

        // Calculate doppler shifted laser frequencies
        auto laserFreqs = this->GetLaserFrequencies() + detunings;
        const auto& propagation = this->GetLasersCounterPropagation();
        auto doppler = velocity / SpeedOfLight_v;
        for (std::size_t i = 0; i < laserFreqs.Size(); i++)
            laserFreqs[i] *= 1 + propagation[i] * doppler;

        // Add light term to the hamiltonian
        TStaticColVector<double, N> photonTerms;
        MatrixMul(photonTerms, m_photonBasis, laserFreqs);
        for (std::size_t i = 0; i < (~hamiltonian).Rows(); i++)
            hamiltonian(i, i) += TwoPi_v * photonTerms(i, 0); 

        return hamiltonian; 
    }

    template<std::size_t N>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemQM<N>::GetNaturalDensityMatrixSS(
        const TColVector<VT>& detunings, double velocity) const
    {
        const TStaticMatrix<std::complex<double>, N, N>& h = GetHamiltonianTI(detunings, velocity);
        TStaticMatrix<std::complex<double>, N*N + 1, N*N> A;

        // von Neumann part of the evolution operator
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
            {
                for (std::size_t k = 0; k < N; k++)
                {
                    A(i*N+k, j*N+k) -= 1.0i * h(i, j);
                    A(k*N+i, k*N+j) += 1.0i * h(j, i);
                }
            }
        }
        
        // Lindblad part of the evolution operator
        for(const auto& decay: this->GetDecays())
        {
            auto lvls = decay.first;
            std::size_t i = lvls.first;
            std::size_t f = lvls.second;
            auto rate = TwoPi_v * decay.second;
            A(f*N+f, i*N+i) += rate;
            for (std::size_t j = 0; j < N; j++)
            {
                A(j*N+i, j*N+i) -= 0.5 * rate;
                A(i*N+j, i*N+j) -= 0.5 * rate;
            }
        }

        // find order of magnitude
        double oom = 0.0;
        for (std::size_t i = 0; i < A.Rows(); i++)
        {
            for (std::size_t j = 0; j < A.Cols(); j++)
                oom = std::max(oom, std::abs(A(i, j)));
        }       

        // Normalization condition
        for (std::size_t i = 0; i < N; i++)
            A(N*N, i*N+i) += oom;
        
        TStaticMatrix<std::complex<double>, N*N + 1, 1> b;
        b(N*N, 0) = oom;
        
        auto A_adj = A.Adjoint();
        auto x = LinearSolve(A_adj*A, A_adj*b);

        TStaticMatrix<std::complex<double>, N, N> ss;
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < N; j++)
                ss(i, j) = x(i*N+j, 0);
        }
        
        return TStaticDensityMatrix<N>(this->GetLevelNames(), ss);
    }
    
    template<std::size_t N>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemQM<N>::GetDensityMatrixSS(
        const TColVector<VT>& detunings) const
    {
        auto& doppler = this->GetDopplerIntegrator();
        auto ss = doppler.Integrate([&](double vel)
            { return this->GetNaturalDensityMatrixSS(detunings, vel).GetMatrix(); });
        return TStaticDensityMatrix<N>(this->GetLevelNames(), ss);
    }

    template<std::size_t N>
    bool TNLevelSystemQM<N>::PrepareCalculation()
    {
        // calculate transition splittings
        auto laserFrequencies = this->GetLaserFrequencies();

        // create photon basis matrix
        m_photonBasis.Resize(N, laserFrequencies.Size());
        m_photonBasis.SetZero();
        std::vector<std::size_t> trans_path;
        std::set<std::size_t> visited_levels;
        trans_path.reserve(laserFrequencies.Size());

        while (visited_levels.size() < N)
        {
            std::size_t head = 0;
            for(; visited_levels.count(head) != 0; head++);

            if (!PreparePhotonBasis(trans_path, visited_levels, head))
            {
                m_photonBasis.Resize(0, 0);
                return false;
            }
        }

        // Atom hamiltonian
        m_hamiltonianNoLight.SetZero();
        for (std::size_t i = 0; i < N; i++)
            m_hamiltonianNoLight(i, i) += TwoPi_v * this->GetLevel(i);
        
        // Interaction hamiltonian
        for (std::size_t i = 0; i < this->GetLaserCount(); i++)
        {
            auto lvls = this->GetLaserLevels(i);
            auto rabi = 0.5 * this->GetLaserElectricField(i) * 
                this->GetDipoleElement(lvls.first, lvls.second) / ReducedPlanckConstant_v;
            m_hamiltonianNoLight(lvls.first, lvls.second) += rabi;
            m_hamiltonianNoLight(lvls.second, lvls.first) += std::conj(rabi);
        }

        return true;
    }

    template<std::size_t N>
    bool TNLevelSystemQM<N>::PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
        std::set<std::size_t>& visitedLvls, std::size_t transFrom)
    {
        // cehck for circular transition path
        if (visitedLvls.count(transFrom) != 0)
            return false;
        else
            visitedLvls.insert(transFrom);

        // Iterate over all transitions
        for (std::size_t tIdx = 0; tIdx < this->GetLaserCount(); tIdx++)
        {
            // skip if the transition was already in the current transition path
            if (std::find(trans_path.begin(), trans_path.end(), tIdx) != trans_path.end())
                continue;

            // skip if currLevel is not involved in the current transition
            // otherwise check to which level the transition leads
            auto lvls = this->GetLaserLevels(tIdx);
            std::size_t transTo = 0;
            if (lvls.first == transFrom)
                transTo = lvls.second;
            else if (lvls.second == transFrom)
                transTo = lvls.first;
            else
                continue;
            
            // check whether a photon is absorbed or emitted
            double relPhotonNumber = 1;
            if (this->GetLevel(transTo) > this->GetLevel(transFrom))
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

    template<std::size_t N>
    template<typename VT>
    typename TNLevelSystemQM<N>::HamiltonianType TNLevelSystemQM<N>::GetHamiltonianAux(
        const TColVector<VT>& detunings, double velocity) const
    {
        return GetHamiltonianTI(detunings, velocity);
    }

}

#endif
