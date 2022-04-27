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
    class TNLevelSystemQM : public TNLevelSystemCRTP<N, TNLevelSystemQM<N>, false>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemQM<N>, false>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemQM<N>, false>;
    public:
        // parent type
        using Laser_t = typename MyParent::Laser_t;

        // constructors
        //template<int dummy=N, typename=EnableDefaultCtor_t<dummy>>
        TNLevelSystemQM() : MyParent() { }
        TNLevelSystemQM(unsigned int dims) : MyParent(dims) { }
        
        // copy operations
        TNLevelSystemQM(const TNLevelSystemQM&) = default;
        TNLevelSystemQM& operator=(const TNLevelSystemQM&) = default;

        // hamiltonian
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianAux(
            const Eigen::VectorXd& laserFreqs) const;
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianFast(
            const Eigen::Matrix<std::complex<double>, N, N>& auxData, double t) const;
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianTI(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;

        // Steady state
        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrixSS(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;
        
    private:
        // Helper methods to prepare the hamiltonian calculation
        void PrepareHamiltonian();
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<unsigned int>& trans_path, 
            std::set<unsigned int>& visitedLevels, unsigned int transFrom);

        bool OnLaserAdded(const Laser_t&) { return PrepareCalculation(); }
        void OnLaserRemoved() { PrepareCalculation(); }
        void OnDipoleOperatorChanged() { PrepareHamiltonian(); }
        
    private:
        Eigen::Matrix<double, N, Eigen::Dynamic> m_photonBasis;
        Eigen::Matrix<std::complex<double>, N, N> m_hamiltonianNoLight;
    };
 

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianAux(
        const Eigen::VectorXd& laserFreqs) const
    {
        return GetHamiltonianTI(laserFreqs);
    }

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianFast(
        const Eigen::Matrix<std::complex<double>, N, N>& auxData, double t) const 
    { 
        return auxData; 
    }

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetHamiltonianTI(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {        
        Eigen::Matrix<std::complex<double>, N, N> hamiltonian = m_hamiltonianNoLight;
        hamiltonian.diagonal() += TwoPi_v * (m_photonBasis * laserFreqs);
        return hamiltonian;
    }

    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemQM<N>::GetDensityMatrixSS(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {
        auto h = GetHamiltonianTI(laserFreqs);
        
        using ATy = std::conditional_t<TNLevelSystemQM<N>::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, N*N>, Eigen::MatrixXcd>;
        using BTy = std::conditional_t<TNLevelSystemQM<N>::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, 1>, Eigen::VectorXcd>;
        using SSTy = Eigen::Matrix<std::complex<double>, N, N>;
        
        unsigned int dims = this->GetDims();
        ATy A = ATy::Zero(dims*dims + 1, dims*dims);

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
            
            double ratePi = Pi_v * rate;
            A(f*dims+f, i*dims+i) += 2 * ratePi;
            for (unsigned int j = 0; j < dims; j++)
            {
                A(j*dims+i, j*dims+i) -= ratePi;
                A(i*dims+j, i*dims+j) -= ratePi;
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
        auto x = ((A.adjoint()*A).ldlt().solve(A.adjoint()*b)).eval();
        return Eigen::Map<SSTy>(x.data(), dims, dims).transpose();
    }
    
    template<int N>
    void TNLevelSystemQM<N>::PrepareHamiltonian()
    {
        // Atom hamiltonian
        m_hamiltonianNoLight = TwoPi_v * this->GetLevels().asDiagonal();

        // Interaction hamiltonian
        for (unsigned int i = 0; i < this->GetLaserCount(); i++)
        {
            constexpr double hbarHalf = 0.5 / ReducedPlanckConstant_v;
            auto& laser = this->GetLaser(i);
            auto [l1, l2] = laser.GetLevels();
            auto rabiHalf = hbarHalf * laser.GetElAmplitude() * this->GetDipoleElement(l1, l2);
            m_hamiltonianNoLight(l1, l2) += rabiHalf;
            m_hamiltonianNoLight(l2, l1) += std::conj(rabiHalf);
        }

        // Set energy zero
        std::size_t dims = this->GetDims();
        m_hamiltonianNoLight -= m_hamiltonianNoLight(0,0) * Eigen::Matrix<std::complex<double>, N, N>::Identity(dims, dims);
    }

    template<int N>
    bool TNLevelSystemQM<N>::PrepareCalculation()
    {
        // initialize photon basis matrix which can be used for the calculation
        // of the light field contribution to the hamiltonian
        m_photonBasis.setZero(this->GetDims(), this->GetLaserCount());

        std::vector<unsigned int> trans_path;
        std::set<unsigned int> visited_levels;
        trans_path.reserve(this->GetLaserCount());

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

        PrepareHamiltonian();

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
            auto [l1, l2] = this->GetLaser(tIdx).GetLevels();
            unsigned int transTo = 0;
            if (l1 == transFrom)
                transTo = l2;
            else if (l2 == transFrom)
                transTo = l1;
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

}

#endif
