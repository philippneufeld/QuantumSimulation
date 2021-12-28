// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystemSC_H_
#define QSim_NLevelSystemSC_H_

#include <cstdint>
#include <memory>

#include "../Math/Matrix.h"
#include "NLevelSystem.h"

namespace QSim
{
    //
    // Solver for a semi-classical model of a N-level system that is driven by light fields
    //

    template<std::size_t N>
    class TNLevelSystemSC : public TNLevelSystemCRTP<N, TNLevelSystemSC<N>>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemSC<N>>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemSC<N>>;
    public:
        // constructors
        TNLevelSystemSC() : MyParent() { }
        TNLevelSystemSC(const std::array<double, N>& levels) : MyParent(levels) { }
        TNLevelSystemSC(const std::array<std::string, N>& lvlNames) 
            : MyParent(lvlNames) { }
        TNLevelSystemSC(const std::array<std::string, N>& lvlNames, 
            const std::array<double, N>& levels) : MyParent(lvlNames, levels) { }

        // copy operations
        TNLevelSystemSC(const TNLevelSystemSC&) = default;
        TNLevelSystemSC& operator=(const TNLevelSystemSC&) = default;  

    private:
        bool OnLaserAdded(std::size_t lvl1, std::size_t lvl2, bool counter) { return true; }
        void OnLaserRemoved() { }

        using HamiltonianType = TStaticMatrix<std::complex<double>, N, N>;

        template<typename VT>
        TStaticColVector<double, N> CalculateRotatingFrame(
            const TColVector<VT>& laserFreqencies) const;

        // auxilliary data that can be prepared once during a time series and be reused later
        using HAuxData = std::tuple<
            // hamiltonian without time-dependent terms:
            HamiltonianType, 
            // Time dependent laser fields (lvl1, lvl2, rabi, (angular) freq):
            std::vector<std::tuple<std::size_t, std::size_t, std::complex<double>, double>>
            >;

        template<typename VT>
        HAuxData GetHamiltonianAux(const TColVector<VT>& detunings, double velocity) const;
        HamiltonianType GetHamiltonianFast(const HAuxData& auxData, double t) const;   
    };

    template<std::size_t N>
    template<typename VT>
    TStaticColVector<double, N> TNLevelSystemSC<N>::CalculateRotatingFrame(
        const TColVector<VT>& laserFreqencies) const
    {
        auto levels = this->GetLevels();
        TStaticColVector<double, N> frame;
        
        frame[0] = levels[0];
        for (std::size_t i = 1; i < N; i++)
        {
            // find closest matching laser frequency
            double closest = 0.0;
            for (std::size_t j = 0; j < (~laserFreqencies).Size(); j++)
            {
                if (std::abs(levels[i] - (~laserFreqencies)[j]) < std::abs(levels[i] - closest))
                    closest = (~laserFreqencies)[j];
            }
            frame[i] = frame[i-1] + closest;    
        }

        return frame;
    }

    template<std::size_t N>
    template<typename VT>
    typename TNLevelSystemSC<N>::HAuxData TNLevelSystemSC<N>::GetHamiltonianAux(
        const TColVector<VT>& detunings, double velocity) const
    {
        assert((~detunings).Rows() == this->GetLaserCount());
        
        HamiltonianType hamiltonian(N, N);

        // Atom hamiltonian
        auto levels = this->GetLevels();
        for (std::size_t i = 0; i < N; i++)
            hamiltonian(i, i) = TwoPi_v * levels[i];  

        // Calculate doppler shifted laser frequencies
        auto laserFreqs = this->GetLaserFrequencies() + detunings;
        const auto& propagation = this->GetLasersCounterPropagation();
        auto doppler = velocity / SpeedOfLight_v;
        for (std::size_t i = 0; i < laserFreqs.Size(); i++)
            laserFreqs[i] *= 1 + propagation[i] * doppler;
        
        // Rotating frame
        auto frame = this->CalculateRotatingFrame(laserFreqs);
        
        // apply roatating frame to hamiltonian 
        // (additional term that comes from the temporal derivative of rho in the rotating frame)
        for (std::size_t i = 0; i < N; i++)
            hamiltonian(i, i) -= TwoPi_v * frame(i);

        // Calculate electric field and system-light interaction
        // Time dependent laser fields (lvl1, lvl2, rabi, (angular) freq):
        std::vector<std::tuple<std::size_t, std::size_t, std::complex<double>, double>> laserFields;
        for (std::size_t i = 0; i < (~laserFreqs).Size(); i++)
        {
            double E0 = this->GetLaserElectricField(i);    
            for (std::size_t j = 0; j < N; j++)
            {
                for (std::size_t k = j + 1; k < N; k++)
                {
                    std::complex<double> rabi = 0.5 * E0 * this->GetDipoleElement(k, j) / ReducedPlanckConstant_v;
                    if (std::abs(rabi) == 0.0)
                        continue;

                    // calculate frequencies of the electric field components
                    double elFieldFreq = (~laserFreqs)(i);
                    double frame_kj = frame(k) - frame(j);
                    double freqs[] = { frame_kj + elFieldFreq, frame_kj - elFieldFreq };

                    std::complex<double> electricField = 0.0;

                    // Rotating wave approximation: Ignore counter rotating propagating wave
                    // if both frequencies are equal (in absolute terms) then keep both
                    double maxFreq = std::abs(freqs[0]) < std::abs(freqs[1]) ? std::abs(freqs[0]) : std::abs(freqs[1]);
                    for (double freq: freqs)
                    {
                        if (std::abs(freq) > maxFreq)
                            continue;

                        // if the contribution is not time dependent, 
                        // it can directly be added to the hamiltonian
                        if (freq != 0.0)
                        {
                            laserFields.emplace_back(k, j, rabi, TwoPi_v*freq);
                        }
                        else
                        {
                            hamiltonian(k, j) += rabi;
                            hamiltonian(j, k) += std::conj(rabi);
                        }
                    }
                }
            }
        }

        // return auxilliary data
        return std::make_tuple(hamiltonian, std::move(laserFields));
    }
    
    template<std::size_t N>
    typename TNLevelSystemSC<N>::HamiltonianType TNLevelSystemSC<N>::GetHamiltonianFast(
        const HAuxData& auxData, double t) const
    {
        // unpack auxilliary data
        HamiltonianType hamiltonian = std::get<0>(auxData);
        const auto& laserFields = std::get<1>(auxData);

        for (const auto& laserField: laserFields)
        {
            std::size_t lvl1 = std::get<0>(laserField);
            std::size_t lvl2 = std::get<1>(laserField);
            auto rabi = std::get<2>(laserField);
            double angFreq = std::get<3>(laserField);
            
            auto coupling = rabi * std::exp(std::complex<double>(1.0i * angFreq * t));
            
            hamiltonian(lvl1, lvl2) += coupling;
            hamiltonian(lvl2, lvl1) += std::conj(coupling);
        }

        return hamiltonian;
    }

}

#endif
