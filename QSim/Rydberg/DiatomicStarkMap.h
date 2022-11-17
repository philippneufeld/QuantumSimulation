// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_DiatomStarkMap_H_
#define QSim_Rydberg_DiatomStarkMap_H_

#include <Eigen/Dense>

#include <vector>
#include <memory>

#include "RydbergDiatomic.h"

#include "../Execution/SingleThreaded.h"
#include "../Execution/ThreadPool.h"

namespace QSim
{

    class DiatomicStarkMap
    {
    public:
        DiatomicStarkMap(std::unique_ptr<RydbergDiatomic> pRydberg, 
            const std::vector<int>& Rs, int nMax, int mN, double energy, double dE);
        
        // getter
        const std::vector<RydbergDiatomicState_t>& GetBasis() const { return m_basis; }
        const Eigen::MatrixXd& GetDipoleOperator() const { return m_dipoleOperator; }

        void PrepareCalculation(const SingleThreaded& = SingleThreaded{});
        void PrepareCalculation(ThreadPool& pool);

        // eigen-energy calculator
        Eigen::VectorXd GetEnergies(double electricField);
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> GetEnergiesAndStates(double electricField);
    
    private:
        void CalcDiagHamiltonianTerms();
        void CalcOffDiagHamiltonianTermsRow(int row);
        void SymmetrizeHamiltonianTerms();

    private:
        std::unique_ptr<RydbergDiatomic> m_pRydberg;
        std::vector<RydbergDiatomicState_t> m_basis;

        bool m_bPrepared; // indicates if following variables are valid
        Eigen::VectorXd m_energies;
        Eigen::MatrixXd m_hamiltonian0;
        Eigen::MatrixXd m_dipoleOperator;
    };
    
}

#endif
