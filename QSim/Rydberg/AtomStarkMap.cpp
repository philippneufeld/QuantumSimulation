// Philipp Neufeld, 2021-2022

#include "AtomStarkMap.h"

namespace QSim
{
    AtomStarkMap::AtomStarkMap(
        const TRydbergSystem<RydbergAtomState_t>& system, 
        const RydbergAtomState_t& state, int nMin, int nMax, int lMax)
    {
        auto [n, l, j, mj] = state;

        if (nMin > n || nMax < n || l > n || 
            std::abs(mj) > j || n < 1 || nMin < 1 || l < 0)
            throw std::runtime_error("Invalid quantum numbers");

        // generate basis
        for (int n = nMin; n <= nMax; n++)
        {
            for (int l = 0; l <= n && l <= lMax; l++)
            {
                int jmult = 2;
                for (int i=0; i<jmult; i++)
                {
                    double j = double(l) - 0.5 + i;
                    if (std::abs(mj) - 0.1 < j)
                    {
                        m_basis.emplace_back(n, l, j, mj);
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

    std::pair<Eigen::VectorXd, Eigen::VectorXd> AtomStarkMap::GetEnergies(double electricField)
    {
        Eigen::MatrixXd hamiltonian = electricField * m_dipoleOperator;
        hamiltonian += m_energies.asDiagonal();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(hamiltonian);

        Eigen::VectorXd energies = solver.eigenvalues();
        Eigen::MatrixXd states = solver.eigenvectors();
        Eigen::VectorXd overlaps(energies.size());

        // column k of states is the k-th eigen-vector
        // the l-th element of the eigen vector is the amount of 
        // overlap of the eigen vector and the l-th basis state
        for (int i = 0; i < overlaps.size(); i++)
            overlaps[i] = std::abs(states(m_referenceStateIdx, i));

        return std::make_pair(energies, overlaps);
    }
}
