// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystem2_H_
#define QSim_NLevelSystem2_H_

#include <complex>
#include <array>
#include <vector>
#include <map>
#include <type_traits>
#include <Eigen/Dense>

#include "../Constants.h"
#include "../Math/Ode.h"
#include "Laser.h"

namespace QSim
{
    // Enables the i literal for complex numbers
    using namespace std::complex_literals;

    constexpr int DynamicDim_v = Eigen::Dynamic;

    template<int N, typename MyT>
    class TNLevelSystemCRTP2
    {
        struct CouplingLaser
        {
            std::pair<unsigned int, unsigned int> lvls;
            double E0;
            bool counterProp;
        };

    public:

        template<int compileTimeSize>
        using EnableDefaultCtor_t = std::enable_if_t<(compileTimeSize != DynamicDim_v)>;

        //template<int dummy=N, typename=EnableDefaultCtor_t<dummy>>
        TNLevelSystemCRTP2() : TNLevelSystemCRTP2(N) {}
        TNLevelSystemCRTP2(unsigned int dims);

        // copy operations
        TNLevelSystemCRTP2(const TNLevelSystemCRTP2&) = default;
        TNLevelSystemCRTP2& operator=(const TNLevelSystemCRTP2&) = default;

        static constexpr bool IsStaticDim() { return (N != DynamicDim_v); }
        unsigned int GetDims() const { return m_levels.size(); }

        // CRTP operators
        MyT& operator~() { return static_cast<MyT&>(*this); }
        const MyT& operator~() const { return static_cast<const MyT&>(*this); }

        // levels
        const auto& GetLevels() const { return m_levels; }
        double GetLevel(unsigned int idx) const;
        bool SetLevel(unsigned int idx, double level);

        // decay rates due to spontaneous emission
        const auto& GetDecays() const { return m_decays; }
        double GetDecay(unsigned int from, unsigned int to) const;
        bool SetDecay(unsigned int from, unsigned int to, double rate);

        // Transition dipole operator
        const auto& GetDipoleOperator() const { return m_dipoleOp; }
        std::complex<double> GetDipoleElement(unsigned int from, unsigned int to) const;
        bool SetDipoleElement(unsigned int from, unsigned int to, std::complex<double> dip);

        // coupling laser
        bool AddLaser(unsigned int lvl1, unsigned int lvl2, double intensity, bool counter);

        unsigned int GetLaserCount() const { return m_lasers.size(); }
        const auto& GetLaserFrequencies() const { return m_laserFreqs; }
        Eigen::VectorXd GetDopplerLaserFreqs(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const;

        std::pair<unsigned int, unsigned int> GetLaserLevels(unsigned int idx) const;
        double GetLaserElectricAmplitude(unsigned int idx) const;
        double GetLaserIntensity(unsigned int idx) const;
        bool GetLaserCounterPropagation(unsigned int idx) const;

        bool SetLaserIntensity(unsigned int idx, double intensity);

        // hamiltonian
        Eigen::Matrix<std::complex<double>, N, N> GetHamiltonian(
            const Eigen::Ref<const Eigen::VectorXd>& detunings,
            double velocity, double t) const;

        /*// time evolution of density matrix
        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrix(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, 
            const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
            double velocity, double t0, double t, double dt);

        Eigen::Matrix<std::complex<double>, N, N> GetDensityMatrixAv(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, 
            const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
            double velocity, double t0, double t, double tav, double dt);

        std::pair<Eigen::VectorXd, std::vector<Eigen::Matrix<std::complex<double>, N, N>>>
        GetTrajectory(
            const Eigen::Ref<const Eigen::VectorXd>& detunings, 
            const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
            double velocity, double t0, double t, double dt);
        
        // create stecific density matrices
        Eigen::Matrix<std::complex<double>, N, N> CreateGroundState() const;
        Eigen::Matrix<std::complex<double>, N, N> CreateThermalState(double temperature) const;*/
        
    private:
        // auxilliary laser variable update
        void UpdateLaserFrequencies();

        /*// helper methods for time evolution
        template<typename AuxType>
        Eigen::Matrix<std::complex<double>, N, N> GetDensityOpDerivative(
            const AuxType& auxData, 
            const Eigen::Matrix<std::complex<double>, N, N>& rho,
            double t) const;

        template<typename AuxType>
        Eigen::Matrix<std::complex<double>, N, N> EvolveDensityMatrix(
            const AuxType& auxData, const Eigen::Matrix<std::complex<double>, N, N>& rho0, 
            double t0, double dt, unsigned int steps);*/

    private:

        // Properties of the system
        Eigen::Vector<double, N> m_levels;
        Eigen::Matrix<std::complex<double>, N, N> m_dipoleOp;
        std::map<std::pair<unsigned int, unsigned int>, double> m_decays;

        // coupling lasers
        std::vector<CouplingLaser> m_lasers;

        // laser related performance enhancing auxilliary variables 
        Eigen::VectorXd m_laserFreqs;
    };


    template<int N, typename MyT>
    TNLevelSystemCRTP2<N, MyT>::TNLevelSystemCRTP2(unsigned int dims) 
        : m_levels(dims), m_dipoleOp(dims, dims) 
    {
        m_levels.setZero();
        m_dipoleOp.setZero();
    }

    template<int N, typename MyT>
    double TNLevelSystemCRTP2<N, MyT>::GetLevel(unsigned int idx) const
    {
        return idx < GetDims() ? m_levels(idx) : 0.0;
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::SetLevel(unsigned int idx, double level)
    {
        if (idx >= GetDims())
            return false;
        m_levels[idx] = level;
        UpdateLaserFrequencies();
        return true;
    }

    template<int N, typename MyT>
    double TNLevelSystemCRTP2<N, MyT>::GetDecay(
        unsigned int from, unsigned int to) const
    {
        auto it = m_decays.find(std::make_pair(from, to));
        return it != m_decays.end() ? it->second : 0.0;
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::SetDecay(
        unsigned int from, unsigned int to, double rate)
    {
        if (from >= GetDims() || to >= GetDims())
            return false;  // index out of bound
        m_decays[std::make_pair(from, to)] = rate;
        return true;
    }

    template<int N, typename MyT>
    std::complex<double> TNLevelSystemCRTP2<N, MyT>::GetDipoleElement(
        unsigned int from, unsigned int to) const
    {
        return (from < GetDims() && to < GetDims()) ? 
            m_dipoleOp(from, to) : 0.0;
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::SetDipoleElement(
        unsigned int from, unsigned int to, std::complex<double> dip)
    {
        if (from >= GetDims() || to >= GetDims())
            return false;  // index out of bound
        m_dipoleOp(from, to) = dip;
        m_dipoleOp(to, from) = dip;
        return true;
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::AddLaser(unsigned int lvl1, 
        unsigned int lvl2, double intensity, bool counter)
    {
        if (lvl1 >= GetDims()  || lvl2 >= GetDims() || lvl1 == lvl2)
            return false; // invalid levels

        CouplingLaser laser;
        laser.lvls = std::make_pair(lvl1, lvl2);
        laser.E0 = GetElectricFieldFromIntensity(intensity);
        laser.counterProp = counter;

        m_lasers.push_back(laser);
        UpdateLaserFrequencies();

        // Add laser and remove it again, if OnLaserAdded returns false
        if (!(~(*this)).OnLaserAdded(lvl1, lvl2, counter))
        {
            m_lasers.pop_back();
            UpdateLaserFrequencies();
            (~(*this)).OnLaserRemoved();
            return false;
        }

        return true;
    }

    template<int N, typename MyT>
    Eigen::VectorXd TNLevelSystemCRTP2<N, MyT>::GetDopplerLaserFreqs(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, double velocity) const
    {
        assert(detunings.rows() == this->GetLaserCount()); 
        Eigen::VectorXd laserFreqs = m_laserFreqs + detunings;
        
        double doppler = velocity / SpeedOfLight_v;
        for (unsigned int i = 0; i < m_lasers.size(); i++)
            laserFreqs[i] *= m_lasers[i].counterProp ? 1 + doppler : 1 - doppler;

        return laserFreqs;
    }

    template<int N, typename MyT>
    std::pair<unsigned int, unsigned int> 
        TNLevelSystemCRTP2<N, MyT>::GetLaserLevels(unsigned int idx) const
    {
        if (idx >= m_lasers.size())
            return std::make_pair<unsigned int, unsigned int>(-1, -1);
        return m_lasers[idx].lvls;
    }

    template<int N, typename MyT>
    double TNLevelSystemCRTP2<N, MyT>::GetLaserElectricAmplitude(unsigned int idx) const
    {
        if (idx >= m_lasers.size())
            return 0.0;
        return m_lasers[idx].E0;
    } 

    template<int N, typename MyT>
    double TNLevelSystemCRTP2<N, MyT>::GetLaserIntensity(unsigned int idx) const
    {
        if (idx >= m_lasers.size())
            return 0.0;
        double electricField = m_lasers[idx].E0;
        return GetIntensityFromElectricField(electricField);
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::GetLaserCounterPropagation(unsigned int idx) const
    {
        if (idx >= m_lasers.size())
            return false;
        return m_lasers[idx].counterProp;
    }

    template<int N, typename MyT>
    bool TNLevelSystemCRTP2<N, MyT>::SetLaserIntensity(
        unsigned int idx, double intensity)
    {
        if (idx >= m_lasers.size())
            return false;

        m_lasers[idx].E0 = GetElectricFieldFromIntensity(intensity);
        return true;
    }

    template<int N, typename MyT>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::GetHamiltonian(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, 
        double velocity, double t) const
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        return (~(*this)).GetHamiltonianFast(auxData, t);
    }

    /*template<int N, typename MyT>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::GetDensityMatrix(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, 
        const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
        double velocity, double t0, double t, double dt)
    {
        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil((t-t0) / dt));
        dt = (t-t0) / steps;
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        return this->EvolveDensityMatrix(auxData, rho0, t0, dt, steps);
    }

    template<int N, typename MyT>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::GetDensityMatrixAv(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, 
        const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
        double velocity, double t0, double t, double tav, double dt)
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        
        // validate averaging time and generate starting and 
        // end time of the averaging process
        tav = tav > t - t0 ? t - t0 : tav;
        double t1 = t - tav / 2;
        double t2 = t1 + tav;

        // evolve up to the starting point of the averaging
        unsigned int steps1 = static_cast<unsigned int>(std::ceil((t1-t0) / dt));
        double dt1 = (t1 - t0) / steps1;
        auto rho = this->EvolveDensityMatrix(auxData, rho0, t0, dt1, steps1);

        // continue evolving while averaging the newly calculated density matrices
        unsigned int steps2 = static_cast<unsigned int>(std::ceil((t2-t1) / dt));
        double dt2 = (t2 - t1) / steps2;

        auto rhoAv = rho;
        for (unsigned int i = 0; i < steps2; i++)
        {
            rho = this->EvolveDensityMatrix(
                auxData, rho, t1 + i*dt2, dt2, 1);
            rhoAv += rho;
        }
        rhoAv *= 1.0 / (steps2 + 1);

        return rhoAv;
    }
    
    template<int N, typename MyT>
    std::pair<Eigen::VectorXd, std::vector<Eigen::Matrix<std::complex<double>, N, N>>> TNLevelSystemCRTP2<N, MyT>::GetTrajectory(
        const Eigen::Ref<const Eigen::VectorXd>& detunings, 
        const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0,
        double velocity, double t0, double t, double dt)
    {
        unsigned int steps = static_cast<unsigned int>(std::ceil((t-t0) / dt));

        std::vector<Eigen::Matrix<std::complex<double>, N, N>> trajectory;
        trajectory.reserve(steps + 1);
        trajectory.push_back(rho0);
        
        auto rho = rho0;
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        for (unsigned int i = 0; i < steps; i++)
        {
            rho = this->EvolveDensityMatrix(
                auxData, rho, t0 + i*dt, dt, 1);
            trajectory.push_back(rho);
        }

        return {Eigen::VectorXd::LinSpaced(steps, t0, t0 + steps*dt), trajectory};
    }

    template<int N, typename MyT>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::CreateGroundState() const
    { 
        return CreateThermalState(0.0);
    }

    template<int N, typename MyT>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::CreateThermalState(double temperature) const
    { 
        unsigned int dims = GetDims();
        Eigen::Matrix<std::complex<double>, N, N> res(dims, dims);

        Eigen::Vector<double, N> lvlOffset(dims);
        unsigned int gsIdx = std::min_element(m_levels.Data(), m_levels.Data() + dims) - m_levels.Data();
        for (unsigned int i = 0; i < dims; i++)
            lvlOffset(i) = m_levels(i) - m_levels(gsIdx);

        Eigen::Vector<double, N> thermalWeights(dims);
        if (temperature <= 0.0)
        {
            for (unsigned int i = 0; i < dims; i++)
            {
                if (lvlOffset(i) == 0.0)
                    thermalWeights(i) = 1.0;
            }
        }
        else
        {
            double tmp = PlanckConstant_v / (BoltzmannConstant_v * temperature);
            for (unsigned int i = 0; i < dims; i++)
                thermalWeights(i) = std::exp(-tmp * lvlOffset(i));
        }

        double norm = 0.0;
        for (unsigned int i = 0; i < dims; i++)
            norm += thermalWeights(i);
        norm = 1.0 / norm;
        thermalWeights *= norm;

        for (unsigned int i = 0; i < dims; i++)
            res(i, i) = thermalWeights(i);
        return res;
    }*/

    template<int N, typename MyT>
    void TNLevelSystemCRTP2<N, MyT>::UpdateLaserFrequencies()
    {
        m_laserFreqs.resize(GetLaserCount());
        for (unsigned int i = 0; i < m_laserFreqs.size(); i++)
        {
            auto [l1, l2] = m_lasers[i].lvls;
            m_laserFreqs[i] = std::abs(m_levels[l1] - m_levels[l2]);
        }
    }

    /*template<int N, typename MyT>
    template<typename AuxType>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::GetDensityOpDerivative(
        const AuxType& auxData, 
        const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho, 
        double t) const
    {
        // von Neumann term
        auto h = (~(*this)).GetHamiltonianFast(auxData, t);
        Eigen::Matrix<std::complex<double>, N, N> rhoPrime = -1.0i * (h * rho - rho * h);

        // add lindblad dissipation term
        for (const auto& decay: this->m_decays)
        {
            auto idxPair = decay.first; // from, to
            auto rate = TwoPi_v * decay.second;
            std::complex<double> popDecayRate = rate * rho(idxPair.first, idxPair.first);
            rhoPrime(idxPair.first, idxPair.first) -= popDecayRate;
            rhoPrime(idxPair.second, idxPair.second) += popDecayRate;
            rhoPrime(idxPair.first, idxPair.second) -= 0.5 * rate * rho(idxPair.first, idxPair.second);
            rhoPrime(idxPair.second, idxPair.first) -= 0.5 * rate * rho(idxPair.second, idxPair.first);
        }

        return rhoPrime;
    }

    template<int N, typename MyT>
    template<typename AuxType>
    Eigen::Matrix<std::complex<double>, N, N> TNLevelSystemCRTP2<N, MyT>::EvolveDensityMatrix(
        const AuxType& auxData, const Eigen::Ref<const Eigen::Matrix<std::complex<double>, N, N>>& rho0, 
        double t0, double dt, unsigned int steps)
    {
        using YType = Eigen::Matrix<std::complex<double>, N, N>;
        ODERK4<double, YType> integrator;
        YType rho = rho0;

        for (unsigned int i = 0; i < steps; i++)
        {
            auto func = [&](double x, const YType& y) 
            { 
                return this->GetDensityOpDerivative(auxData, y, x);
            };
            rho += integrator.Step(func, rho, t0 + i * dt, dt);       
        }

        return rho;
    }*/
}

// Include specific implementations
#include "NLevelSystemQM2.h"

#endif
