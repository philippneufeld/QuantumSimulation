// Philipp Neufeld, 2021-2022

#include "DiatomicStarkMap.h"
#include <iostream>

namespace QSim
{
    DiatomicStarkMap::DiatomicStarkMap(
        const RydbergDiatomic& system, 
        int nMin, int nMax, int RMin, int RMax, int mN, 
        double centerEnergy, double maxEnergyDist)
    {
        // generate basis
        for (int n = nMin; n <= nMax; n++)
        {
            for (int l = 0; l < n; l++)
            {
                for (int R = RMin; R <= RMax; R+=1)
                {
                    for (int N = std::abs(R-l); N <= R+l; N++)
                    {
                        // only add states with same mN as reference state
                        if (std::abs(mN) <= N)
                        {
                            RydbergDiatomicState_t s(n, l, R, N, mN);
                            if (std::abs(system.GetEnergy(s) - centerEnergy) <= maxEnergyDist)
                                m_basis.emplace_back(s);
                        }
                    }
                }
            }
        }

        int stateCnt = m_basis.size();

        // get energies
        m_hamiltonian0 = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        for (int i=0; i<stateCnt; i++)
            m_hamiltonian0(i, i) = system.GetEnergy(m_basis[i]);

        // calculate dipole operator (and self-interaction operator)
        m_dipoleOperator = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        for (int i1 = 0; i1 < stateCnt; i1++)
        {
            for (int i2 = 0; i2 <= i1; i2++)
            {
                // self dipole and self multielectron interaction
                double selfDip = system.GetSelfDipoleME(m_basis[i1], m_basis[i2]);
                double selfEl = system.GetCoreInteractionME(m_basis[i1], m_basis[i2]);
                m_hamiltonian0(i1, i2) += selfEl + selfDip;
                m_hamiltonian0(i2, i1) = m_hamiltonian0(i1, i2);

                // stark operators
                double dip = system.GetDipoleME(m_basis[i1], m_basis[i2]);
                m_dipoleOperator(i1, i2) += dip;
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
        Eigen::MatrixXd hamiltonian = m_hamiltonian0 + electricField * m_dipoleOperator;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(hamiltonian);
        Eigen::VectorXd energies = solver.eigenvalues();
        Eigen::MatrixXd states = solver.eigenvectors();
        
        // column k of states is the k-th eigen-vector
        return std::make_pair(energies, states);
    }
}
