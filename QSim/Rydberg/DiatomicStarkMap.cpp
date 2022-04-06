// Philipp Neufeld, 2021-2022

#include "DiatomicStarkMap.h"

namespace QSim
{
    DiatomicStarkMap::DiatomicStarkMap(
        const TRydbergSystem<RydbergDiatomicState_t>& system, 
        const RydbergDiatomicState_t& state, int nMin, int nMax, int RMax, double maxEnergyDist)
    {
        m_referenceStateIdx = -1;

        auto [n, l, R, N, mN] = state;
        double energy = system.GetEnergy(state);

        // generate basis
        for (int n = nMin; n <= nMax; n++)
        {
            for (int l = 0; l <= n; l++)
            {
                for (int R = 0; R <= RMax; R++)
                {
                    for (int N = std::abs(R-l); N <= R+l; N++)
                    {
                        // only add states with same mN as reference state
                        if (std::abs(mN) <= N)
                        {
                            RydbergDiatomicState_t s(n, l, R, N, mN);
                            
                            if (std::abs(system.GetEnergy(s) - energy) <= maxEnergyDist)
                                m_basis.emplace_back(s);
                        }
                    }
                }
            }
        }

        auto it = std::find(m_basis.begin(), m_basis.end(), state);
        if (it == m_basis.end()) m_basis.push_back(state);
        m_referenceStateIdx = it - m_basis.begin();

        int stateCnt = m_basis.size();

        // get energies
        m_energies = Eigen::VectorXd(stateCnt);
        for (int i=0; i<stateCnt; i++)
            m_energies[i] = system.GetEnergy(m_basis[i]);

        // calculate dipole operator
        m_dipoleOperator = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        for (int i1 = 0; i1 < stateCnt; i1++)
        {
            for (int i2 = 0; i2 < i1; i2++)
            {
                double dip = system.GetDipoleME(m_basis[i1], m_basis[i2]);
                m_dipoleOperator(i1, i2) = dip;
                m_dipoleOperator(i2, i1) = dip;
            }
        }
    }

    Eigen::VectorXd DiatomicStarkMap::GetEnergies(double electricField)
    {
        return GetEnergiesAndStates(electricField).first;
    }

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> DiatomicStarkMap::GetEnergiesAndStates(double electricField)
    {
        Eigen::MatrixXd hamiltonian = electricField * m_dipoleOperator;
        hamiltonian += m_energies.asDiagonal();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(hamiltonian);

        Eigen::VectorXd energies = solver.eigenvalues();
        Eigen::MatrixXd states = solver.eigenvectors();
        
        // column k of states is the k-th eigen-vector
        return std::make_pair(energies, states);
    }
}
