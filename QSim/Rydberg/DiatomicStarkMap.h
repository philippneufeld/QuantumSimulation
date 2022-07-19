// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_DiatomStarkMap_H_
#define QSim_Rydberg_DiatomStarkMap_H_

#include <Eigen/Dense>

#include <vector>

#include "RydbergDiatomic.h"


namespace QSim
{

    class DiatomicStarkMap
    {
    public:
        DiatomicStarkMap(
            const RydbergDiatomic& system, 
            int nMin, int nMax, int RMax, int mN,
            double centerEnergy, double maxEnergyDist);
        
        // getter
        const std::vector<RydbergDiatomicState_t>& GetBasis() const { return m_basis; }
        const Eigen::MatrixXd& GetDipoleOperator() const { return m_dipoleOperator; }

        // eigen-energy calculator
        Eigen::VectorXd GetEnergies(double electricField);
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> GetEnergiesAndStates(double electricField);
    
    private:
        std::vector<RydbergDiatomicState_t> m_basis;

        Eigen::VectorXd m_energies;
        Eigen::MatrixXd m_hamiltonian0;
        Eigen::MatrixXd m_dipoleOperator;
    };
    
}

#endif
