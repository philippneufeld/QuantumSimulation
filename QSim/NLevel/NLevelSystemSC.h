// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystemQM_H_
#define QSim_NLevelSystemQM_H_

#include <cstdint>

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
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemSC<N>>;
    public:
        // constructors
        TNLevelSystemSC() : TNLevelSystemCRTP<N, TNLevelSystemSC<N>>() { }
        TNLevelSystemSC(const std::array<double, N>& levels) 
            : TNLevelSystemCRTP<N, TNLevelSystemSC<N>>(levels) { }
        TNLevelSystemSC(const std::array<std::string, N>& lvlNames) 
            : TNLevelSystemCRTP<N, TNLevelSystemSC<N>>(lvlNames) { }
        TNLevelSystemSC(const std::array<std::string, N>& lvlNames, const std::array<double, N>& levels) 
            : TNLevelSystemCRTP<N, TNLevelSystemSC<N>>(lvlNames, levels) { }

        // copy operations
        TNLevelSystemSC(const TNLevelSystemSC&) = default;
        TNLevelSystemSC& operator=(const TNLevelSystemSC&) = default;

        template<typename VT>
        TStaticMatrix<std::complex<double>, N, N> GetHamiltonian(
            const TColVector<VT>& detunings, double velocity, double t) const;
    };

    template<std::size_t N>
    template<typename VT>
    TStaticMatrix<std::complex<double>, N, N> TNLevelSystemSC<N>::GetHamiltonian(
        const TColVector<VT>& detunings, double velocity, double t) const
    {
        // The returned hamiltonian is in units of angular frequency
        assert((~detunings).Rows() == this->GetLaserCount());
        
        using HamiltonianType = TStaticMatrix<std::complex<double>, N, N>;
        HamiltonianType hamiltonian(N, N);

        // Atom hamiltonian
        auto angularFreqLevels = TwoPi_v * this->m_levels;
        for (std::size_t i = 0; i < N; i++)
            hamiltonian(i, i) = angularFreqLevels[i];  

        // Calculate doppler shifted laser frequencies
        auto laserFreqs = TwoPi_v * (this->GetLaserFrequencies() + detunings);
        auto laserFreqsDoppler = laserFreqs;
        for (std::size_t i = 0; i < laserFreqsDoppler.Size(); i++)
        {
            auto doppler = velocity / SpeedOfLight_v;
            laserFreqsDoppler[i] *= this->GetLaserCounterPropagation(i) ? 1.0 + doppler : 1.0 - doppler;
        }
        
        // Rotating frame
        TStaticColVector<double, N> frameFrequencies;
        frameFrequencies[0] = angularFreqLevels[0];
        for (std::size_t i = 1; i < N; i++)
        {
            // find closest matching laser frequency
            double closest = 0.0;
            for (std::size_t j = 0; j < (~laserFreqsDoppler).Size(); j++)
            {
                if (std::abs(angularFreqLevels[i] - (~laserFreqsDoppler)[j]) < std::abs(angularFreqLevels[i] - closest))
                    closest = (~laserFreqsDoppler)[j];
            }
            frameFrequencies[i] = frameFrequencies[i-1] + closest;    
        }

        // apply roatating frame to hamiltonian 
        // (additional term that comes from the temporal derivative of rho in the rotating frame)
        for (std::size_t i = 0; i < N; i++)
            hamiltonian(i, i) -= frameFrequencies(i);
        
        // Calculate electric field and system-light interaction
        double zc = velocity * t / SpeedOfLight_v;
        for (std::size_t i = 0; i < (~laserFreqsDoppler).Size(); i++)
        {
            double E0 = this->GetLaserElectricField(i);    
            for (std::size_t j = 0; j < N; j++)
            {
                for (std::size_t k = j + 1; k < N; k++)
                {
                    // calculate frequencies of the electric field components
                    double elFieldFreq = (~laserFreqsDoppler)(i);
                    double frame_kj = frameFrequencies(k) - frameFrequencies(j);
                    double freqs[] = { frame_kj + elFieldFreq, frame_kj - elFieldFreq };

                    std::complex<double> electricField = 0.0;

                    // Rotating wave approximation: Ignore counter rotating propagating wave
                    // if both frequencies are equal (in absolute terms) then keep both
                    double maxFreq = std::abs(freqs[0]) < std::abs(freqs[1]) ? std::abs(freqs[0]) : std::abs(freqs[1]);
                    for (double freq: freqs)
                    {
                        if (std::abs(freq) <= maxFreq) 
                            electricField += 0.5* E0 * (freq != 0.0 ? std::exp(std::complex<double>(1.0i * freq * t)) : 1.0);
                    }

                    electricField /= ReducedPlanckConstant_v;
                    hamiltonian(k, j) += electricField * this->m_dipoleOperator(k, j);
                    hamiltonian(j, k) += std::conj(electricField) * this->m_dipoleOperator(j, k);
                }
            }
        }

        return hamiltonian;
    }

}

#endif
