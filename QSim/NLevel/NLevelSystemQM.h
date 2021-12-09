// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystemQM_H_
#define QSim_NLevelSystemQM_H_

#include <cstdint>
#include <memory>

#include "../Math/Matrix.h"
#include "NLevelSystem.h"

namespace QSim
{
    //
    // Solver for a quantum mechanical Jaynes-Cummings model of a 
    // N-level system that is driven by light fields
    //

    template<std::size_t N>
    class TNLevelSystemQM : public TNLevelSystemCRTP<N, TNLevelSystemQM<N>>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemQM<N>>;
    public:
        // constructors
        TNLevelSystemQM() : MyParent() { }
        TNLevelSystemQM(const std::array<double, N>& levels) : MyParent(levels) { }
        TNLevelSystemQM(const std::array<std::string, N>& lvlNames) 
            : MyParent(lvlNames) { }
        TNLevelSystemQM(const std::array<std::string, N>& lvlNames, 
            const std::array<double, N>& levels) : MyParent(lvlNames, levels) { }

        // copy operations
        TNLevelSystemQM(const TNLevelSystemQM&) = default;
        TNLevelSystemQM& operator=(const TNLevelSystemQM&) = default;  

    private:

        // Helper methods to prepare the hamiltonian calculation
        bool PrepareCalculation();
        bool PreparePhotonBasis(std::vector<std::size_t>& trans_path, 
            std::set<std::size_t>& visitedLevels, std::size_t transFrom);

        bool OnLaserAdded(std::size_t lvl1, std::size_t lvl2, bool counter) { return PrepareCalculation(); }
        void OnLaserRemoved() { PrepareCalculation(); }

        using HamiltonianType = TStaticMatrix<std::complex<double>, N, N>;

        template<typename VT>
        HamiltonianType GetHamiltonianAux(const TColVector<VT>& detunings, double velocity) const;
        HamiltonianType GetHamiltonianFast(const HamiltonianType& auxData, double t) const { return auxData; }   
    };

}

#endif
