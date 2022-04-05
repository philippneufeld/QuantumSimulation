// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_DiatomStarkMap_H_
#define QSim_Rydberg_DiatomStarkMap_H_

#include <Eigen/Dense>

#include <vector>

#include "RydbergDiatomic.h"


namespace QSim
{

    /*class AtomStarkMap
    {
    public:
        AtomStarkMap(
            const TRydbergSystem<RydbergAtomState_t>& system, 
            int n, int l, double j, double mj, 
            int nMin, int nMax, int lMax);
        
        // getter
        const std::vector<std::tuple<int, int, double, double>>& GetBasis() const { return m_basis; }
        const Eigen::MatrixXd& GetDipoleOperator() const { return m_dipoleOperator; }

        // eigen-energy calculator
        // Eigen::VectorXd GetEnergies(double electricField);
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetEnergies(double electricField);

    private:
        int m_referenceStateIdx;
        std::vector<RydbergAtomState_t> m_basis;

        Eigen::VectorXd m_energies;
        Eigen::MatrixXd m_dipoleOperator;
    };*/
    
}

#endif
