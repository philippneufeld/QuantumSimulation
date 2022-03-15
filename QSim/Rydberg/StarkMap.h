// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_StarkMap_H_
#define QSim_Rydberg_StarkMap_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <exception>
#include <vector>
#include <tuple>
#include "HydrogenicSystem.h"

#include <iostream>

namespace QSim
{

    class StarkMap
    {
    public:
        StarkMap(const HydrogenicSystem& system, int n, int l, double j, 
            double mj, int nMin, int nMax, int lMax);
        // getter
        const std::vector<std::tuple<int, int, double, double>>& GetBasis() const { return m_basis; }
        const Eigen::MatrixXd& GetDipoleOperator() const { return m_dipoleOperator; }

        // eigen-energy calculator
        Eigen::VectorXd GetEnergies(double electricField);

    private:
        int m_n, m_l;
        double m_j, m_mj;
        int m_nMin, m_nMax, m_lMax;
        std::vector<std::tuple<int, int, double, double>> m_basis;

        Eigen::VectorXd m_energies;
        Eigen::MatrixXd m_dipoleOperator;
    };
    
}

#endif
