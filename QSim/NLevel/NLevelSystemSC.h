// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystemSC_H_
#define QSim_NLevelSystemSC_H_

#include <cstdint>
#include <memory>

#include <Eigen/Dense>
#include "NLevelSystem.h"

namespace QSim
{
    //
    // Solver for a semi-classical model of a N-level system that is driven by light fields
    //

    template<int N>
    class TNLevelSystemSC : public TNLevelSystemCRTP<N, TNLevelSystemSC<N>>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemSC<N>>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemSC<N>>;
    public:
        // constructors
        TNLevelSystemSC() : MyParent() { }
        TNLevelSystemSC(unsigned int dims) : MyParent(dims) { }

        // copy operations
        TNLevelSystemSC(const TNLevelSystemSC&) = default;
        TNLevelSystemSC& operator=(const TNLevelSystemSC&) = default;  

    private:
        // Allow every possible laser configuration
        bool OnLaserAdded(std::size_t lvl1, std::size_t lvl2, bool counter) { return true; }
        void OnLaserRemoved() { }

        Eigen::Matrix<double, N, 1> CalculateRotatingFrame(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;

        // auxilliary data that can be prepared once during a time series and be reused later
        using HAuxData = std::tuple<
            // hamiltonian without time-dependent terms:
            Eigen::Matrix<std::complex<double>, N, N>, 
            // Time dependent laser fields (lvl1, lvl2, rabi, (angular) freq):
            std::vector<std::tuple<std::size_t, std::size_t, std::complex<double>, double>>
            >;

        HAuxData GetHamiltonianAux(const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const;
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianFast(const HAuxData& auxData, double t) const;   
    };

    template<int N>
    Eigen::Matrix<double, N, 1> TNLevelSystemSC<N>::CalculateRotatingFrame(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {
        auto levels = this->GetLevels();
        Eigen::Matrix<double, N, 1> frame(this->GetDims());
        
        frame[0] = levels[0];
        for (std::size_t i = 1; i < this->GetDims(); i++)
        {
            // find closest matching laser frequency
            int idx = 0;
            auto freqDiff = (laserFreqs.array() - levels[i]).abs().minCoeff(&idx);
            double offset = freqDiff < std::abs(levels[i]) ? laserFreqs[idx] : 0.0;
            frame[i] = frame[i-1] + offset;    
        }

        return frame;
    }

    template<int N>
    typename TNLevelSystemSC<N>::HAuxData TNLevelSystemSC<N>::GetHamiltonianAux(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const
    {
        assert(detunings.size() == this->GetLaserCount());
        
        unsigned int dims = this->GetDims();
        Eigen::Matrix<std::complex<double>, N, N> h0(dims, dims);

        // Atom hamiltonian
        h0 = TwoPi_v * this->GetLevels().asDiagonal();

        // Calculate doppler shifted laser frequencies and apply 
        // roatating frame to hamiltonian (additional term that 
        // comes from the temporal derivative of rho in the rotating frame)
        auto laserFreqs = this->GetDopplerLaserFreqs(detunings, velocity);
        auto frame = this->CalculateRotatingFrame(laserFreqs);
        h0 -= TwoPi_v * frame.asDiagonal();
        
        // Calculate electric field and system-light interaction
        // Time dependent laser fields (lvl1, lvl2, rabi, (angular) freq):
        std::vector<std::tuple<std::size_t, std::size_t, std::complex<double>, double>> laserFields;
        for (std::size_t i = 0; i < laserFreqs.size(); i++)
        {
            double E0 = this->GetLaserElectricAmplitude(i);    
            for (std::size_t j = 0; j < dims; j++)
            {
                for (std::size_t k = j + 1; k < dims; k++)
                {
                    std::complex<double> rabi = 0.5 * E0 * this->GetDipoleElement(k, j) / ReducedPlanckConstant_v;
                    if (std::abs(rabi) == 0.0)
                        continue;

                    // calculate frequencies of the electric field components
                    double elFieldFreq = laserFreqs(i);
                    double frame_kj = frame(k) - frame(j);
                    double freqs[] = { frame_kj + elFieldFreq, frame_kj - elFieldFreq };

                    // Rotating wave approximation: Ignore counter rotating wave
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
                            h0(k, j) += rabi;
                            h0(j, k) += std::conj(rabi);
                        }
                    }
                }
            }
        }

        // return auxilliary data
        return std::make_tuple(h0, std::move(laserFields));
    }
    
    template<int N>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemSC<N>::GetHamiltonianFast(
        const HAuxData& auxData, double t) const
    {
        // unpack auxilliary data
        const auto& [h0, fields] = auxData;
        Eigen::Matrix<std::complex<double>, N, N> h = h0;

        for (const auto& field: fields)
        {
            auto [l1, l2, rabi, angFreq] = field;
            auto coupling = rabi * std::exp(std::complex<double>(1.0i * angFreq * t));
            h(l1, l2) += coupling;
            h(l2, l1) += std::conj(coupling);
        }

        return h;
    }

}

#endif
