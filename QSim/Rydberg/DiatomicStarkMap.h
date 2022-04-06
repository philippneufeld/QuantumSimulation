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
            const TRydbergSystem<RydbergDiatomicState_t>& system, 
            const RydbergDiatomicState_t& state, 
            int nMin, int nMax, int RMax, double maxEnergyDist);
        
        // getter
        int GetReferenceStateIndex() const { return m_referenceStateIdx; }
        const std::vector<RydbergDiatomicState_t>& GetBasis() const { return m_basis; }
        const Eigen::MatrixXd& GetDipoleOperator() const { return m_dipoleOperator; }

        // eigen-energy calculator
        Eigen::VectorXd GetEnergies(double electricField);
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> GetEnergiesAndStates(double electricField);
        
    private:
        int m_referenceStateIdx;
        std::vector<RydbergDiatomicState_t> m_basis;

        Eigen::VectorXd m_energies;
        Eigen::MatrixXd m_dipoleOperator;
    };
    
}

#endif
