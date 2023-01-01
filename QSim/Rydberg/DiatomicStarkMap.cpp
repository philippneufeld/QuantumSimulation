// Philipp Neufeld, 2021-2022

#include "DiatomicStarkMap.h"
#include <iostream>

namespace QSim
{
    DiatomicStarkMap::DiatomicStarkMap(std::unique_ptr<RydbergDiatomic> pRydberg, 
        const std::vector<int>& Rs, int nMax, int mN, double energy, double dE)
        : m_pRydberg(std::move(pRydberg)), m_bPrepared(false)
    {

        dE = std::abs(dE);
        double minEnergy = energy - dE;
        double maxEnergy = energy + dE;

        // generate basis
        for (int R: Rs)
        {
            for (int n = 1; n < nMax; n++)
            {
                bool all_above = true;
                for (int l = 0; l < n; l++)
                {
                    for (int N = std::abs(R-l); N <= R+l; N++)
                    {
                        // only add states with same mN as reference state
                        if (std::abs(mN) <= N)
                        {
                            RydbergDiatomicState_t s(n, l, R, N, mN);
                            double sEnergy = m_pRydberg->GetEnergy(s);

                            if (sEnergy < maxEnergy)
                                all_above = false;

                            // add to basis if in range
                            if (std::abs(sEnergy - energy) <= dE)
                                m_basis.emplace_back(s);
                        }
                    }
                }

                // all states in this manifold were above the desired energy range 
                // -> next manifolds will be even higher in energy
                if (all_above)
                    break;
            }
        }
    }

    void DiatomicStarkMap::PrepareCalculation(const SingleThreaded&)
    {
        int stateCnt = m_basis.size();
        m_hamiltonian0 = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        m_dipoleOperator = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        
        CalcDiagHamiltonianTerms();
        for (int row = 0; row < stateCnt; row++)
            CalcOffDiagHamiltonianTermsRow(row);
        SymmetrizeHamiltonianTerms();
        
        m_bPrepared = true;
    }

    void DiatomicStarkMap::PrepareCalculation(ThreadPool& pool)
    {
        int stateCnt = m_basis.size();
        m_hamiltonian0 = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        m_dipoleOperator = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        
        CalcDiagHamiltonianTerms();
        std::vector<std::future<void>> futures;
        for (int row = 0; row < stateCnt; row++)
        {
            futures.push_back(pool.SubmitWithFuture([&, row=row]()
            {
                CalcOffDiagHamiltonianTermsRow(row);
            }));
        }
        for (auto& future: futures)
            future.wait();

        SymmetrizeHamiltonianTerms();
        
        m_bPrepared = true;
    }

    void DiatomicStarkMap::CalcDiagHamiltonianTerms()
    {
        int stateCnt = m_basis.size();
        m_hamiltonian0 = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        for (int i=0; i<stateCnt; i++)
            m_hamiltonian0(i, i) = m_pRydberg->GetEnergy(m_basis[i]);
    }

    void DiatomicStarkMap::CalcOffDiagHamiltonianTermsRow(int row)
    {
        for (int i2 = 0; i2 <= row; i2++)
        {
            // self dipole and self multielectron interaction
            double selfDip = m_pRydberg->GetSelfDipoleME(m_basis[row], m_basis[i2]);
            double selfEl = m_pRydberg->GetCoreInteractionME(m_basis[row], m_basis[i2]);
            m_hamiltonian0(row, i2) += selfEl + selfDip;
            

            // stark operators
            double dip = m_pRydberg->GetDipoleME(m_basis[row], m_basis[i2]);
            m_dipoleOperator(row, i2) += dip;
        }
    }

    void DiatomicStarkMap::SymmetrizeHamiltonianTerms()
    {
        for (int i1 = 0; i1 < m_basis.size(); i1++)
        {
            for (int i2 = 0; i2 <= i1; i2++)
            {
                m_hamiltonian0(i2, i1) = m_hamiltonian0(i1, i2);
                m_dipoleOperator(i2, i1) = m_dipoleOperator(i1, i2);
            }
        }
    }

    Eigen::VectorXd DiatomicStarkMap::GetEnergies(double electricField)
    {
        return GetEnergiesAndStates(electricField).first;
    }

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> DiatomicStarkMap::GetEnergiesAndStates(double electricField)
    {
        Eigen::MatrixXd hamiltonian = GetHamiltonian(electricField);
        return GetEnergiesAndStatesFromHamiltonian(std::move(hamiltonian));
    }

    Eigen::MatrixXd DiatomicStarkMap::GetHamiltonian(double electricField)
    {
        if (!m_bPrepared)
            PrepareCalculation();

        return m_hamiltonian0 + electricField * m_dipoleOperator;
    }

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> 
        DiatomicStarkMap::GetEnergiesAndStatesFromHamiltonian(const Eigen::MatrixXd& ham)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(ham);
        Eigen::VectorXd energies = solver.eigenvalues();
        Eigen::MatrixXd states = solver.eigenvectors();
        
        // column k of states is the k-th eigen-vector
        return std::make_pair(energies, states);
    }
}
