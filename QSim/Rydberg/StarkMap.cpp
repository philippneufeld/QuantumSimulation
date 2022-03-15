// Philipp Neufeld, 2021-2022

#include "StarkMap.h"

namespace QSim
{
    StarkMap::StarkMap(const HydrogenicSystem& system, int n, int l, 
        double j, double mj, int nMin, int nMax, int lMax)
        : m_n(n), m_l(l), m_j(j), m_mj(mj), m_nMin(nMin), m_nMax(nMax), m_lMax(lMax)
    {
        if (m_nMin > m_n || m_nMax < m_n || m_l > m_n || 
            std::abs(m_mj) > m_j || m_n < 1 || m_nMin < 1 || m_l < 0)
            throw std::runtime_error("Invalid quantum numbers");

        // generate basis
        /*for (int n = m_nMin; n < m_nMax; n++)
        {
            for (int l = 0; l <= n && l <= m_lMax; l++)
            {
                m_basis.emplace_back(n, l, m);
            }
        }*/
        for (int n = m_nMin; n <= m_nMax; n++)
        {
            for (int l = 0; l <= n && l <= m_lMax; l++)
            {
                int jmult = 2;
                for (int i=0; i<jmult; i++)
                {
                    double j = double(l) - 0.5 + i;
                    if (std::abs(mj) - 0.1 < j)
                        m_basis.emplace_back(n, l, j, mj);
                }
                
            }
        }
        int stateCnt = m_basis.size();

        // get energies
        m_energies = Eigen::VectorXd(stateCnt);
        for (int i=0; i<stateCnt; i++)
        {
            auto [n, l, j, m] = m_basis[i];
            m_energies[i] = system.GetEnergy(n);
        }

        // calculate dipole operator
        m_dipoleOperator = Eigen::MatrixXd::Zero(stateCnt, stateCnt);
        for (int i1 = 0; i1 < stateCnt; i1++)
        {
            for (int i2 = 0; i2 < i1; i2++)
            {
                auto [n1, l1, j1, m1] = m_basis[i1];
                auto [n2, l2, j2, m2] = m_basis[i2];

                double dip = system.GetDipoleME2(n1, l1, j1, m1, n2, l2, j2, m2);
                m_dipoleOperator(i1, i2) = dip;
                m_dipoleOperator(i2, i1) = dip;
            }
        }
    }

    Eigen::VectorXd StarkMap::GetEnergies(double electricField)
    {
        Eigen::MatrixXd hamiltonian = electricField * m_dipoleOperator;
        hamiltonian += m_energies.asDiagonal();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(hamiltonian);
        return solver.eigenvalues();
    }
}
