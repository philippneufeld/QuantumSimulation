// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystem_H_
#define QSim_NLevelSystem_H_

#include <complex>
#include <array>
#include <vector>
#include <map>
#include <functional>
#include <type_traits>
#include <Eigen/Dense>

#include "../Util/CRTP.h"
#include "../Constants.h"
#include "../Math/Ode.h"
#include "Laser.h"

namespace QSim
{
    // Enables the i literal for complex numbers
    using namespace std::complex_literals;

    constexpr int DynamicDim_v = Eigen::Dynamic;

    template<int N, typename MyT, bool AM>
    class TNLevelSystemCRTP : public TCRTP<MyT>
    {
    protected:
        ~TNLevelSystemCRTP() = default;

    public:

        using Laser_t = TNLevelLaser<AM>;

        template<int compileTimeSize>
        using EnableDefaultCtor_t = std::enable_if_t<(compileTimeSize != DynamicDim_v)>;

        //template<int dummy=N, typename=EnableDefaultCtor_t<dummy>>
        TNLevelSystemCRTP() : TNLevelSystemCRTP(N) {}
        TNLevelSystemCRTP(unsigned int dims);

        // copy operations
        TNLevelSystemCRTP(const TNLevelSystemCRTP&) = default;
        TNLevelSystemCRTP(TNLevelSystemCRTP&&) = default;
        TNLevelSystemCRTP& operator=(const TNLevelSystemCRTP&) = default;
        TNLevelSystemCRTP& operator=(TNLevelSystemCRTP&&) = default;

        // dimensionality of the system
        static constexpr bool IsStaticDim() { return (N != DynamicDim_v); }
        unsigned int GetDims() const { return m_levels.size(); }

        // levels
        const Eigen::Matrix<double, N, 1>& GetLevels() const { return m_levels; }
        double GetLevel(unsigned int idx) const;
        bool SetLevel(unsigned int idx, double level);

        // decay rates due to spontaneous emission
        const auto& GetDecays() const { return m_decays; }
        double GetDecay(unsigned int from, unsigned int to) const;
        bool SetDecay(unsigned int from, unsigned int to, double rate);
        bool AddDecay(unsigned int from, unsigned int to, double rate);

        // dipole operator
        const Eigen::Matrix<std::complex<double>, N, N>& GetDipoleOperator() const { return m_dipoleOp; }
        std::complex<double> GetDipoleElement(unsigned int from, unsigned int to) const;
        bool SetDipoleElement(unsigned int from, unsigned int to, std::complex<double> dip);

        // coupling laser
        unsigned int GetLaserCount() const { return m_lasers.size(); }
        Laser_t& GetLaser(unsigned int idx) { return m_lasers.at(idx); }
        const Laser_t& GetLaser(unsigned int idx) const { return m_lasers.at(idx); }
        bool AddLaser(Laser_t laser);
        Eigen::VectorXd GetTransitionFreqs() const;
        Eigen::VectorXd GetLaserDirs() const;

        // hamiltonian
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonian(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const;

        // time evolution of density matrix
        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrix(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const Eigen::Matrix<std::complex<double>, N, N>& rho0,
            double t0, double t, double dt);

        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrixAv(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const Eigen::Matrix<std::complex<double>, N, N>& rho0,
            double t0, double t, double tav, double dt);

        Eigen::VectorXd GetTrajectoryTimeaxis(double t0, double dt, std::size_t n);
        std::vector<Eigen::Matrix<std::complex<double>, N, N>> GetTrajectory(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const Eigen::Matrix<std::complex<double>, N, N>& rho0,
            double t0, double t, double dt);
        
        // create specific density matrices
        Eigen::Matrix<std::complex<double>, N, N> CreateGroundState() const;
        Eigen::Matrix<std::complex<double>, N, N> CreateThermalState(double temperature) const;
        
    private:
        // helper methods for time evolution
        template<typename AuxType>
        Eigen::Matrix<std::complex<double>, N, N> GetDensityOpDerivative(
            const AuxType& auxData, 
            const Eigen::Matrix<std::complex<double>, N, N>& rho,
            double t) const;

    private:
        // Properties of the system
        Eigen::Matrix<double, N, 1> m_levels;
        Eigen::Matrix<std::complex<double>, N, N> m_dipoleOp;
        std::map<std::pair<unsigned int, unsigned int>, double> m_decays;

        // coupling lasers
        std::vector<Laser_t> m_lasers;
    };


    template<int N, typename MyT, bool AM>
    TNLevelSystemCRTP<N, MyT, AM>::TNLevelSystemCRTP(unsigned int dims) 
        : m_levels(dims), m_dipoleOp(dims, dims) 
    {
        m_levels.setZero();
        m_dipoleOp.setZero();
    }

    template<int N, typename MyT, bool AM>
    double TNLevelSystemCRTP<N, MyT, AM>::GetLevel(unsigned int idx) const
    {
        return idx < GetDims() ? m_levels(idx) : 0.0;
    }

    template<int N, typename MyT, bool AM>
    bool TNLevelSystemCRTP<N, MyT, AM>::SetLevel(unsigned int idx, double level)
    {
        if (idx >= GetDims())
            return false;
        m_levels[idx] = level;
        return true;
    }

    template<int N, typename MyT, bool AM>
    double TNLevelSystemCRTP<N, MyT, AM>::GetDecay(
        unsigned int from, unsigned int to) const
    {
        auto it = m_decays.find(std::make_pair(from, to));
        return it != m_decays.end() ? it->second : 0.0;
    }

    template<int N, typename MyT, bool AM>
    bool TNLevelSystemCRTP<N, MyT, AM>::SetDecay(
        unsigned int from, unsigned int to, double rate)
    {
        if (from >= GetDims() || to >= GetDims())
            return false; // index out of bound
        m_decays[std::make_pair(from, to)] = rate;
        return true;
    }

    template<int N, typename MyT, bool AM>
    bool TNLevelSystemCRTP<N, MyT, AM>::AddDecay(
        unsigned int from, unsigned int to, double rate)
    {
        return SetDecay(from, to, GetDecay(from, to) + rate);
    }

    template<int N, typename MyT, bool AM>
    std::complex<double> TNLevelSystemCRTP<N, MyT, AM>::GetDipoleElement(
        unsigned int from, unsigned int to) const
    {
        return (from < GetDims() && to < GetDims()) ? 
            m_dipoleOp(from, to) : 0.0;
    }

    template<int N, typename MyT, bool AM>
    bool TNLevelSystemCRTP<N, MyT, AM>::SetDipoleElement(
        unsigned int from, unsigned int to, std::complex<double> dip)
    {
        if (from >= GetDims() || to >= GetDims())
            return false;  // index out of bound
        m_dipoleOp(from, to) = dip;
        m_dipoleOp(to, from) = std::conj(dip);

        (~(*this)).OnDipoleOperatorChanged();

        return true;
    }

    template<int N, typename MyT, bool AM>
    bool TNLevelSystemCRTP<N, MyT, AM>::AddLaser(TNLevelLaser<AM> laser)
    {
        auto [lvl1, lvl2] = laser.GetLevels();
        if (lvl1 >= GetDims() || lvl2 >= GetDims() || lvl1 == lvl2)
            return false; // invalid levels

        m_lasers.push_back(laser);
        
        // Add laser and remove it again, if OnLaserAdded returns false
        const auto& cLaserRef = laser;
        if (!(~(*this)).OnLaserAdded(cLaserRef))
        {
            m_lasers.pop_back();
            (~(*this)).OnLaserRemoved();
            return false;
        }

        return true;
    }

    template<int N, typename MyT, bool AM>
    Eigen::VectorXd TNLevelSystemCRTP<N, MyT, AM>::GetTransitionFreqs() const
    {
        Eigen::VectorXd freqs(GetLaserCount());
        for (unsigned int i = 0; i < freqs.size(); i++)
        {
            auto [l1, l2] = m_lasers[i].GetLevels();
            freqs[i] = std::abs(m_levels[l1] - m_levels[l2]);
        }
        return freqs;
    }

    template<int N, typename MyT, bool AM>
    Eigen::VectorXd TNLevelSystemCRTP<N, MyT, AM>::GetLaserDirs() const
    {
        Eigen::VectorXd dirs(GetLaserCount());
        for (unsigned int i = 0; i < dirs.size(); i++)
            dirs[i] = m_lasers[i].GetPropDirection();
        return dirs;
    }

    template<int N, typename MyT, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::GetHamiltonian(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(laserFreqs);
        return (~(*this)).GetHamiltonianFast(auxData, t);
    }

    template<int N, typename MyT, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::GetDensityMatrix(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const Eigen::Matrix<std::complex<double>, N, N>& rho0,
        double t0, double t, double dt)
    {
        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil((t-t0) / dt));
        dt = (t-t0) / steps;

        const auto auxData = (~(*this)).GetHamiltonianAux(laserFreqs);

        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54DPPolicy> integrator;
        using YType = Eigen::Matrix<std::complex<double>, N, N>;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(auxData, y, x); };

        return integrator.IntegrateTo(func, rho0, t0, t, dt);
    }

    template<int N, typename MyT, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::GetDensityMatrixAv(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const Eigen::Matrix<std::complex<double>, N, N>& rho0,
        double t0, double t, double tav, double dt)
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(laserFreqs);
        
        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54DPPolicy> integrator;
        using YType = Eigen::Matrix<std::complex<double>, N, N>;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(auxData, y, x); };
        
        // validate averaging time and generate starting and 
        // end time of the averaging process
        tav = tav > t - t0 ? t - t0 : tav;
        double t1 = t - tav / 2;

        // evolve up to the starting point of the averaging
        double dt1 = dt;
        auto rho = integrator.IntegrateTo(func, rho0, t0, t1, dt1);

        // continue evolving while averaging the newly calculated density matrices
        double dtAvStep = dt; // std::min(dt1, dt);
        double dtAv = dt1;
        
        auto rhoAv = rho;
        unsigned int i = 0;
        for (; i*dtAvStep < tav; i++)
        {
            double t0av = t1 + i*dtAvStep;
            double t1av = t1 + (i+1)*dtAvStep;
            rho = integrator.IntegrateTo(func, rho0, t0av, t1av, dtAv);
            rhoAv += rho;
        }
        rhoAv *= 1.0 / (1 + i);

        return rhoAv;
    }
    
    template<int N, typename MyT, bool AM>
    Eigen::VectorXd TNLevelSystemCRTP<N, MyT, AM>::GetTrajectoryTimeaxis(double t0, double dt, std::size_t n)
    {
        if (n==0) return Eigen::VectorXd();
        return Eigen::VectorXd::LinSpaced(n, t0, t0 + (n-1)*dt);
    }

    template<int N, typename MyT, bool AM>
    std::vector<Eigen::Matrix<std::complex<double>, N, N>> TNLevelSystemCRTP<N, MyT, AM>::GetTrajectory(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const Eigen::Matrix<std::complex<double>, N, N>& rho0,
        double t0, double t, double dt)
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(laserFreqs);
        
        // allocate memory for the trajectory data
        std::size_t n = static_cast<std::size_t>(std::ceil((t-t0) / dt));
        std::vector<Eigen::Matrix<std::complex<double>, N, N>> trajectory;
        trajectory.reserve(n + 1);
        trajectory.push_back(rho0);

        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54DPPolicy> integrator;
        using YType = Eigen::Matrix<std::complex<double>, N, N>;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(auxData, y, x); };
        
        auto rho = rho0;
        double dtInner = dt;
        for (std::size_t i = 0; i < n; i++)
        {
            rho = integrator.IntegrateTo(func, rho, t0+i*dt, t0+(i+1)*dt, dtInner);
            trajectory.push_back(rho);
        };

        return trajectory;
    }

    template<int N, typename MyT, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::CreateGroundState() const
    { 
        return CreateThermalState(0.0);
    }

    template<int N, typename MyT, bool AM>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::CreateThermalState(double temperature) const
    { 
        unsigned int dims = GetDims();
        temperature = temperature > 0 ? temperature : 0;

        unsigned int gsIdx;
        m_levels.minCoeff(&gsIdx);
        auto lvlOffset = m_levels.array() - m_levels(gsIdx);
        
        double tmp = PlanckConstant_v / (BoltzmannConstant_v * temperature);
        Eigen::Array<double, N, 1> thermalWeights(dims);
        if (temperature <= 0.0)
            thermalWeights = (lvlOffset > 0.0).select(Eigen::ArrayXd::Zero(dims), 1.0);
        else
            thermalWeights = (-tmp * lvlOffset).exp();
        
        thermalWeights /= thermalWeights.sum();
        return thermalWeights.matrix().asDiagonal();
    }

    template<int N, typename MyT, bool AM>
    template<typename AuxType>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT, AM>::GetDensityOpDerivative(
        const AuxType& auxData, 
        const Eigen::Matrix<std::complex<double>, N, N>& rho, 
        double t) const
    {
        // von Neumann term
        auto h = (~(*this)).GetHamiltonianFast(auxData, t);
        auto rhoPrime = (-1.0i * (h * rho - rho * h)).eval();

        // add lindblad dissipation term
        unsigned int dims = this->GetDims();
        for (const auto& decay: m_decays)
        {
            auto [lvls, rate] = decay;
            auto [i, f] = lvls;
            double ratePi = Pi_v * rate; 

            rhoPrime(f, f) += 2 * ratePi * rho(i, i);
            for (unsigned int j=0; j<dims; j++)
            {
                rhoPrime(i, j) -= ratePi * rho(i, j);
                rhoPrime(j, i) -= ratePi * rho(j, i);
            }
        }

        return rhoPrime;
    }

}

// Include specific implementations
#include "NLevelSystemSC.h"
#include "NLevelSystemQM.h"

#endif
