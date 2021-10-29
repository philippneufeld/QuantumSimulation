// Philipp Neufeld, 2021

#ifndef QSIM_NLevelSystem2_H_
#define QSIM_NLevelSystem2_H_

#include <cstdint>
#include <vector>
#include <map>
#include <string>
#include <tuple>

#include "Matrix.h"
#include "Transition.h"
#include "Decay.h"

namespace QSim
{

    class NLevelSystem
    {
    public:
        NLevelSystem();
        ~NLevelSystem();

        void SetLevel(const std::string& name, double level);
        double GetLevel(const std::string& name) const { return m_levels.at(name); }

        bool AddTransition(const std::string& lvl1, const std::string& lvl2, double rabi);
        bool AddDecay(const std::string& lvlFrom, const std::string& lvlTo, double rabi);

        TDynamicMatrix<double> GetHamiltonian(const TDynamicMatrix<double>& detunings);

    private:
        bool PrepareCalculation();
        bool PreparePhotonBasis(
            std::vector<std::size_t>& trans_path,
            const std::string& transFrom);

    private:
        std::map<std::string, double> m_levels;

        std::vector<std::tuple<std::string, std::string, double>> m_transitions;
        std::vector<std::tuple<std::string, std::string, double>> m_decays;

        // auxilliary variables for the creation of the hamiltonian
        std::vector<std::string> m_usedLevels;
        TDynamicMatrix<double> m_transitionSplittings;
        TDynamicMatrix<double> m_dopplerFactors;
        TDynamicMatrix<double> m_photonBasis;
        TDynamicMatrix<double> m_hamiltonianNoLight;
    };

}

#endif
