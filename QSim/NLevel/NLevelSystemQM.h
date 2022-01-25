// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystemQM_H_
#define QSim_NLevelSystemQM_H_

#include <cstdint>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include "NLevelSystem.h"

namespace QSim
{
    //
    // Solver for a quantum mechanical Jaynes-Cummings model of a 
    // N-level system that is driven by light fields
    //

    template<int N>
    class TNLevelSystemQM : public TNLevelSystemCRTP<N, TNLevelSystemQM<N>>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
    public:

        // constructors
        //template<int dummy=N, typename=EnableDefaultCtor_t<dummy>>
        TNLevelSystemQM() : MyParent() { }
        TNLevelSystemQM(unsigned int dims) : MyParent(dims) { }
        
        // copy operations
        TNLevelSystemQM(const TNLevelSystemQM&) = default;
        TNLevelSystemQM& operator=(const TNLevelSystemQM&) = default;

        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianTI(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const;

        // Steady state
        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrixSS(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const;
        
    private:

        // Helper methods to prepare the hamiltonian calculation
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<unsigned int>& trans_path, 
            std::set<unsigned int>& visitedLevels, unsigned int transFrom);

        bool OnLaserAdded(unsigned int lvl1, unsigned int lvl2, bool counter) { return PrepareCalculation(); }
        void OnLaserRemoved() { PrepareCalculation(); }

        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianAux(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const;
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianFast(
            const Eigen::Matrix<std::complex<double>, N, N>& auxData, double t) const;

    private:
        Eigen::Matrix<double, N, Eigen::Dynamic> m_photonBasis;
        Eigen::Matrix<std::complex<double>, N, N> m_hamiltonianNoLight;
    };
 
    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianTI(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const
    {
        auto laserFreqs = this->GetDopplerLaserFreqs(detunings, velocity);
        
        Eigen::Matrix<std::complex<double>, N, N> hamiltonian = m_hamiltonianNoLight;
        hamiltonian.diagonal() += TwoPi_v * (m_photonBasis * laserFreqs);
        return hamiltonian;
    }

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetDensityMatrixSS(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const
    {
        auto h = GetHamiltonianTI(detunings, velocity);
        
        using ATy = std::conditional_t<TNLevelSystemQM<N>::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, N*N>, Eigen::MatrixXcd>;
        using BTy = std::conditional_t<TNLevelSystemQM<N>::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, 1>, Eigen::VectorXcd>;
        using SSTy = Eigen::Matrix<std::complex<double>, N, N>;
        
        unsigned int dims = this->GetDims();
        ATy A(dims*dims + 1, dims*dims);

        // von Neumann part of the evolution operator
        for (unsigned int i = 0; i < dims; i++)
        {
            for (unsigned int j = 0; j < dims; j++)
            {
                for (unsigned int k = 0; k < dims; k++)
                {
                    A(i*dims+k, j*dims+k) -= 1.0i * h(i, j);
                    A(k*dims+i, k*dims+j) += 1.0i * h(j, i);
                }
            }
        }
        
        // Lindblad part of the evolution operator
        for(const auto& decay: this->GetDecays())
        {
            auto [lvls, rate] = decay;
            auto [i, f] = lvls;
            rate *= TwoPi_v;
            
            A(f*dims+f, i*dims+i) += rate;
            for (unsigned int j = 0; j < dims; j++)
            {
                A(j*dims+i, j*dims+i) -= 0.5 * rate;
                A(i*dims+j, i*dims+j) -= 0.5 * rate;
            }
        }

        // Normalization condition
        for (unsigned int i = 0; i < dims; i++)
            A(dims*dims, i*dims+i) += 1.0;

        // bring all columns to the same order of magnitude by
        // dividing every row by its abs-max value
        auto Ascale = (1 / A.array().cwiseAbs().rowwise().maxCoeff()).eval();
        A = Ascale.matrix().asDiagonal() * A;

        // generate right hand side of the equation A*x=b
        // b = (0, 0, ..., 0, 1)
        auto b = BTy::Unit(dims*dims + 1, dims*dims);
        
        // solve and subsequently reshape from vetor to matrix
        auto x = ((A.adjoint()*A).householderQr().solve(A.adjoint()*b)).eval();
        return Eigen::Map<SSTy>(x.data(), dims, dims).transpose();
    }
    
    template<int N>
    bool TNLevelSystemQM<N>::PrepareCalculation()
    {
        // calculate transition splittings
        auto laserFrequencies = this->GetLaserFrequencies();

        // initialize photon basis matrix which can be used for the calculation
        // of the light field contribution to the hamiltonian
        m_photonBasis.setZero(this->GetDims(), laserFrequencies.size());

        std::vector<unsigned int> trans_path;
        std::set<unsigned int> visited_levels;
        trans_path.reserve(laserFrequencies.size());

        while (visited_levels.size() < this->GetDims())
        {
            unsigned int head = 0;
            for(; visited_levels.count(head) != 0; head++);

            if (!PreparePhotonBasis(trans_path, visited_levels, head))
            {
                m_photonBasis.resize(0, 0);
                return false;
            }
        }

        // Atom hamiltonian
        m_hamiltonianNoLight = TwoPi_v * this->GetLevels().asDiagonal();

        // Interaction hamiltonian
        for (unsigned int i = 0; i < this->GetLaserCount(); i++)
        {
            auto lvls = this->GetLaserLevels(i);
            auto rabi = 0.5 * this->GetLaserElectricAmplitude(i) * 
                this->GetDipoleElement(lvls.first, lvls.second) / ReducedPlanckConstant_v;
            m_hamiltonianNoLight(lvls.first, lvls.second) += rabi;
            m_hamiltonianNoLight(lvls.second, lvls.first) += std::conj(rabi);
        }

        return true;
    }

    template<int N>
    bool TNLevelSystemQM<N>::PreparePhotonBasis(std::vector<unsigned int>& trans_path, 
        std::set<unsigned int>& visitedLvls, unsigned int transFrom)
    {
        // m_photonBasis(i, j) contains the (relative) photon number
        // of the j-th laser field in the i-th atomic state

        // check for circular transition path
        if (visitedLvls.count(transFrom) != 0)
            return false;
        
        visitedLvls.insert(transFrom);

        // Iterate over all transitions
        for (unsigned int tIdx = 0; tIdx < this->GetLaserCount(); tIdx++)
        {
            // skip if the transition was already in the current transition path
            if (std::find(trans_path.begin(), trans_path.end(), tIdx) != trans_path.end())
                continue;

            // skip if currLevel is not involved in the current transition
            // otherwise check to which level the transition leads
            auto lvls = this->GetLaserLevels(tIdx);
            unsigned int transTo = 0;
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
            m_photonBasis.row(transTo) = m_photonBasis.row(transFrom);
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

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianAux(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const
    {
        return GetHamiltonianTI(detunings, velocity);
    }

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianFast(
        const Eigen::Matrix<std::complex<double>, N, N>& auxData, double t) const 
    { 
        return auxData; 
    }

}

#endif
