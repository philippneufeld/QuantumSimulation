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

    template<int N, bool AM=false>
    class TNLevelSystemSC : public TNLevelSystemCRTP<N, TNLevelSystemSC<N, AM>, AM>
    {
        friend class TNLevelSystemCRTP<N, TNLevelSystemSC<N, AM>, AM>;
        using MyParent = TNLevelSystemCRTP<N, TNLevelSystemSC<N, AM>, AM>;
    public:
        // parent type
        using Laser_t = typename MyParent::Laser_t;

        // constructors
        TNLevelSystemSC() : MyParent() { }
        TNLevelSystemSC(unsigned int dims) : MyParent(dims) { }

        // copy operations
        TNLevelSystemSC(const TNLevelSystemSC&) = default;
        TNLevelSystemSC& operator=(const TNLevelSystemSC&) = default;  

    private:
        // Allow every possible laser configuration
        bool OnLaserAdded(const Laser_t&) { return true; }
        void OnLaserRemoved() { }
        void OnDipoleOperatorChanged() { }

        Eigen::Matrix<double, N, 1> CalculateRotatingFrame(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;

        // auxilliary data that can be prepared once during a time series and be reused later
        using HAuxData = std::tuple<
            // hamiltonian without time-dependent terms
            Eigen::Matrix<std::complex<double>, N, N>,
            // doppler shifted laser frequencies
            Eigen::Matrix<double, N, 1>,
            // rotating frame frequencies
            Eigen::Matrix<double, N, 1>
            >;

        HAuxData GetHamiltonianAux(const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonianFast(const HAuxData& auxData, double t) const;   
    };

    template<int N, bool AM>
    Eigen::Matrix<double, N, 1> TNLevelSystemSC<N, AM>::CalculateRotatingFrame(
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

    template<int N, bool AM>
    typename TNLevelSystemSC<N, AM>::HAuxData TNLevelSystemSC<N, AM>::GetHamiltonianAux(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {
        assert(laserFreqs.size() == this->GetLaserCount());
        
        unsigned int dims = this->GetDims();
        Eigen::Matrix<std::complex<double>, N, N> h0(dims, dims);

        // Atom hamiltonian
        h0 = TwoPi_v * this->GetLevels().asDiagonal();

        // Apply rotating frame to hamiltonian (additional term that 
        // comes from the temporal derivative of rho in the rotating frame)
        auto frame = this->CalculateRotatingFrame(laserFreqs);
        h0 -= TwoPi_v * frame.asDiagonal();

        // return auxilliary data
        return std::make_tuple(h0, laserFreqs, frame);
    }
    
    template<int N, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemSC<N, AM>::GetHamiltonianFast(
        const HAuxData& auxData, double t) const
    {
        // unpack auxilliary data
        const auto& [h0, laserFreqs, frame] = auxData;
        Eigen::Matrix<std::complex<double>, N, N> h = h0;

        // Calculate time dependent system-light interaction
        for (std::size_t i = 0; i < laserFreqs.size(); i++)
        {
            double elAmp = this->GetLaser(i).GetModElAmplitude(t);
            for (std::size_t j = 0; j < this->GetDims(); j++)
            {
                for (std::size_t k = j + 1; k < this->GetDims(); k++)
                {
                    constexpr double hbarHalf = 0.5 / ReducedPlanckConstant_v;
                    std::complex<double> rabiHalf = elAmp * this->GetDipoleElement(k, j) * hbarHalf;
                    if (std::abs(rabiHalf) == 0.0)
                        continue;

                    // calculate frequencies of the electric field components
                    double elFieldFreq = laserFreqs[i];
                    double frame_kj = frame[k] - frame[j];
                    double freqs[] = { frame_kj + elFieldFreq, frame_kj - elFieldFreq };

                    // Rotating wave approximation: Ignore counter rotating wave
                    // if both frequencies are equal (in absolute terms) then keep both
                    double smallerFreq = std::min(std::abs(freqs[0]), std::abs(freqs[1]));
                    for (double freq: freqs)
                    {
                        if (std::abs(freq) > smallerFreq)
                            continue;

                        // add coupling to the hamiltonian
                        std::complex<double> wave = (freq != 0.0 ? std::exp(1.0i * TwoPi_v * freq * t) : 1.0);
                        std::complex<double> coupling = rabiHalf * wave;
                        
                        h(k, j) += coupling;
                        h(j, k) += std::conj(coupling);
                    }
                }
            }
        }

        return h;
    }

}

#endif
