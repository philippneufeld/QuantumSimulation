// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystem_H_
#define QSim_NLevelSystem_H_

#include <cstdint>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <algorithm>
#include <complex>

#include "../Math/Matrix.h"
#include "../Math/Integrator.h"
#include "../Constants.h"
#include "../Doppler.h"
#include "../DensityMatrix.h"

namespace QSim
{
    // Enables the i literal for complex numbers
    using namespace std::complex_literals;


    //
    // N-level quantum system solver (CRTP base class)
    //

    template<std::size_t N, typename MyT>
    class TNLevelSystemCRTP
    {
        using IndexPair = std::pair<std::size_t, std::size_t>;

        // name, lvl1, lvl2, electric field amplitude, counter-propagating
        using CouplingLaser = std::tuple<std::string, std::size_t, std::size_t, double, bool>;
    public:
        // constructors
        TNLevelSystemCRTP();
        TNLevelSystemCRTP(const std::array<double, N>& levels);
        TNLevelSystemCRTP(const std::array<std::string, N>& lvlNames);
        TNLevelSystemCRTP(const std::array<std::string, N>& lvlNames, const std::array<double, N>& levels);

        // copy operations
        TNLevelSystemCRTP(const TNLevelSystemCRTP&) = default;
        TNLevelSystemCRTP& operator=(const TNLevelSystemCRTP&) = default;

        // CRTP operators
        MyT& operator~() { return static_cast<MyT&>(*this); }
        const MyT& operator~() const { return static_cast<const MyT&>(*this); }

        // level names
        std::string GetLevelNameByIndex(std::size_t idx) const;
        std::size_t GetLevelIndexByName(const std::string& name) const;
        bool SetLevelName(std::size_t idx, const std::string& newName);

        // levels
        const TStaticColVector<double, N>& GetLevels() const { return m_levels; }
        double GetLevel(std::size_t idx) const;
        double GetLevelByName(const std::string& name) const;
        bool SetLevel(std::size_t idx, double level);
        bool SetLevelByName(const std::string& name, double level);

        // decay rates due to spontaneous emission
        double GetDecay(std::size_t from, std::size_t to) const;
        double GetDecayByName(const std::string& from, std::string& to) const;
        bool SetDecay(std::size_t from, std::size_t to, double rate);
        bool SetDecayByName(const std::string& from, const std::string& to, double rate);

        // Transition dipole operator
        const TStaticMatrix<std::complex<double>, N, N>& GetDipoleOperator() const { return m_dipoleOperator; }
        std::complex<double> GetDipoleElement(std::size_t from, std::size_t to) const;
        std::complex<double> GetDipoleElementByName(const std::string& from, const std::string& to) const;
        bool SetDipoleElement(std::size_t from, std::size_t to, std::complex<double> dip);
        bool SetDipoleElementByName(const std::string& from, const std::string& to, std::complex<double> rate);

        // coupling laser
        std::size_t GetLaserIdxByName(const std::string& name) const;
        std::size_t GetLaserCount() const { return m_couplingLasers.size(); }
        IndexPair GetLaserLevels(std::size_t idx) const;
        double GetLaserIntensity(std::size_t idx) const;
        double GetLaserElectricField(std::size_t idx) const;
        bool GetLaserCounterPropagation(std::size_t idx) const;
        IndexPair GetLaserLevelsByName(const std::string& name) const;
        double GetLaserIntensityByName(const std::string& name) const;
        double GetLaserElectricFieldByName(const std::string& name) const;
        bool GetLaserCounterPropagationByName(const std::string& name) const;
        TDynamicColVector<double> GetLaserFrequencies() const;
        bool AddLaser(const std::string& name, std::size_t lvl1, 
            std::size_t lvl2, double intensity, bool counter);
        bool AddLaserByName(const std::string& name, const std::string& lvl1, 
            const std::string& lvl2, double intensity, bool counter);
        bool RemoveLaser(const std::string& name);
        bool SetLaserIntensity(const std::string& name, double intensity);

        // thermal environment and properties needed for the doppler integration
        double GetMass() const { return m_doppler.GetMass(); }
        double GetTemperature() const { return m_doppler.GetTemperature(); }
        void SetMass(double mass) { m_doppler.SetMass(mass); }
        void SetTemperature(double temp) { m_doppler.SetTemperature(temp); }

        double GetDopplerWidth(std::size_t from, std::size_t to) const;
        double GetDopplerWidthByName(const std::string& from, const std::string& to) const;

        // create stecific density matrices
        TStaticDensityMatrix<N> CreateGroundState() const;
        TStaticDensityMatrix<N> CreateThermalState() const;

        // time evolution of density matrix
        template<typename VT>
        TStaticDensityMatrix<N> GetNaturalDensityMatrix(
            const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
            double velocity, double t0, double t, double dt);

        template<typename VT>
        TStaticDensityMatrix<N> GetNaturalDensityMatrixAv(
            const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
            double velocity, double t0, double t, double tav, double dt);

        template<typename VT> 
        std::pair<TDynamicColVector<double>, std::vector<TStaticDensityMatrix<N>>> 
        GetNaturalTrajectory(
            const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
            double velocity, double t0, double t, double dt);
        
        // time evolution of doppler broadened density matrix
        template<typename VT>
        TStaticDensityMatrix<N> GetDensityMatrix(
            const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
            double t0, double t, double dt);

        template<typename VT>
        TStaticDensityMatrix<N> GetDensityMatrixAv(
            const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
            double t0, double t, double tav, double dt);

    private:
        // helper method for thermal state creation
        TStaticDensityMatrix<N> CreateThermalStateHelper(double temperature) const;

        // helper methods for default initialization
        std::array<double, N> GenerateDefaultLevels() const;
        std::array<std::string, N> GenerateDefaultLevelNames() const;

        // helper methods for time evolution
        template<typename AuxType>
        TStaticMatrix<std::complex<double>, N, N> GetDensityOpDerivative(
            const AuxType& auxData, 
            const TStaticMatrix<std::complex<double>, N, N>& rho,
            double t) const;

        template<typename AuxType>
        TStaticDensityMatrix<N> EvolveNaturalDensityMatrix(
            const AuxType& auxData, const TStaticDensityMatrix<N>& initial, 
            double t0, double dt, std::size_t steps);

    private:
        // Properties of the system
        std::array<std::string, N> m_levelNames;
        TStaticColVector<double, N> m_levels;
        std::map<IndexPair, double> m_decays;
        TStaticMatrix<std::complex<double>, N, N> m_dipoleOperator;

        // coupling lasers
        std::vector<CouplingLaser> m_couplingLasers;

        // thermal environment
        TDopplerIntegrator<double> m_doppler;
    };

    template<std::size_t N, typename MyT>
    TNLevelSystemCRTP<N, MyT>::TNLevelSystemCRTP()
        : TNLevelSystemCRTP(GenerateDefaultLevelNames(), GenerateDefaultLevels()) { }

    template<std::size_t N, typename MyT>
    TNLevelSystemCRTP<N, MyT>::TNLevelSystemCRTP(const std::array<double, N>& levels)
        : TNLevelSystemCRTP(GenerateDefaultLevelNames(), levels) { }

    template<std::size_t N, typename MyT>
    TNLevelSystemCRTP<N, MyT>::TNLevelSystemCRTP(const std::array<std::string, N>& lvlNames)
        : TNLevelSystemCRTP(lvlNames, GenerateDefaultLevels()) { }

    template<std::size_t N, typename MyT>
    TNLevelSystemCRTP<N, MyT>::TNLevelSystemCRTP(const std::array<std::string, N>& lvlNames, 
        const std::array<double, N>& levels)
        : m_levelNames(lvlNames), m_doppler(AtomicMassUnit_v, 300.0)
    {
        for (size_t i = 0; i < N; i++)
            m_levels(i) = levels[i];
    }

    template<std::size_t N, typename MyT>
    std::string TNLevelSystemCRTP<N, MyT>::GetLevelNameByIndex(std::size_t idx) const
    {
        return idx < N ? m_levelNames[idx] : std::string();
    }

    template<std::size_t N, typename MyT>
    std::size_t TNLevelSystemCRTP<N, MyT>::GetLevelIndexByName(
        const std::string& name) const
    {
        auto it = std::find(m_levelNames.begin(), m_levelNames.end(), name);
        return it != m_levelNames.end() ? it - m_levelNames.begin() : -1;
    }
    
    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetLevelName(
        std::size_t idx, const std::string& newName)
    {
        if (GetLevelIndexByName(newName) < N)
            return false;  // level name already present
        m_levelNames[idx] = newName;
        return true;
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLevel(std::size_t idx) const
    {
        return idx < N ? m_levels[idx] : 0.0;
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLevelByName(const std::string& name) const
    {
        return GetLevel(GetLevelIndexByName(name));
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetLevel(std::size_t idx, double level)
    {
        if (idx >= N)
            return false;
        m_levels[idx] = level;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetLevelByName(
        const std::string& name, double level)
    {
        return SetLevel(GetLevelIndexByName(name), level);
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetDecay(
        std::size_t from, std::size_t to) const
    {
        auto it = m_decays.find(std::make_pair(from, to));
        return it != m_decays.end() ? it->second : 0.0;
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetDecayByName(
        const std::string& from, std::string& to) const
    {
        return GetDecay(GetLevelIndexByName(from), 
            GetLevelIndexByName(to));
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetDecay(
        std::size_t from, std::size_t to, double rate)
    {
        if (from >= N || to >= N)
            return false;  // index out of bound
        m_decays[std::make_pair(from, to)] = rate;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetDecayByName(
        const std::string& from, const std::string& to, double rate)
    {
        return SetDecay(GetLevelIndexByName(from), 
            GetLevelIndexByName(to), rate);
    }

    template<std::size_t N, typename MyT>
    std::complex<double> TNLevelSystemCRTP<N, MyT>::GetDipoleElement(
        std::size_t from, std::size_t to) const
    {
        return (from < N && to < N) ? m_dipoleOperator(from, to) : 0.0;
    }

    template<std::size_t N, typename MyT>
    std::complex<double> TNLevelSystemCRTP<N, MyT>::GetDipoleElementByName(
        const std::string& from, const std::string& to) const
    {
        return GetDipoleElement(GetLevelIndexByName(from), 
            GetLevelIndexByName(to));
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetDipoleElement(
        std::size_t from, std::size_t to, std::complex<double> dip)
    {
        if (from >= N || to >= N)
            return false;  // index out of bound
        m_dipoleOperator(from, to) = dip;
        m_dipoleOperator(to, from) = dip;
        return true;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetDipoleElementByName(
        const std::string& from, const std::string& to, std::complex<double> dip)
    {
        return SetDipoleElement(GetLevelIndexByName(from), 
            GetLevelIndexByName(to), dip);
    }

    template<std::size_t N, typename MyT>
    std::size_t TNLevelSystemCRTP<N, MyT>::GetLaserIdxByName(const std::string& name) const
    {
        auto it = m_couplingLasers.begin();
        for (;it != m_couplingLasers.end() && std::get<0>(*it) != name; it++);
        return it != m_couplingLasers.end() ? it - m_couplingLasers.begin() : -1;
    }

    template<std::size_t N, typename MyT>
    typename TNLevelSystemCRTP<N, MyT>::IndexPair 
        TNLevelSystemCRTP<N, MyT>::GetLaserLevels(std::size_t idx) const
    {
        return (idx < m_couplingLasers.size()) ? 
            std::make_pair(std::get<1>(m_couplingLasers[idx]), std::get<2>(m_couplingLasers[idx])) 
            : std::make_pair<std::size_t, std::size_t>(-1, -1);
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLaserIntensity(std::size_t idx) const
    {
        return GetIntensityFromElectricField(GetLaserElectricField());
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLaserElectricField(std::size_t idx) const
    {
        return (idx < m_couplingLasers.size()) ? std::get<3>(m_couplingLasers[idx]) : 0.0;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::GetLaserCounterPropagation(std::size_t idx) const
    {
        return (idx < m_couplingLasers.size()) ? std::get<4>(m_couplingLasers[idx]) : false;
    }

    template<std::size_t N, typename MyT>
    typename TNLevelSystemCRTP<N, MyT>::IndexPair 
        TNLevelSystemCRTP<N, MyT>::GetLaserLevelsByName(const std::string& name) const
    {
        return GetLaserLevels(GetLaserIdxByName(name));
    }
    
    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLaserIntensityByName(
        const std::string& name) const 
    {
        return GetLaserIntensity(GetLaserIdxByName(name));
    }
    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetLaserElectricFieldByName(
        const std::string& name) const
    {
        return GetLaserElectricField(GetLaserIdxByName(name));
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::GetLaserCounterPropagationByName(
        const std::string& name) const
    {
        return GetLaserCounterPropagation(GetLaserIdxByName(name));
    }

    template<std::size_t N, typename MyT>
    TDynamicColVector<double> TNLevelSystemCRTP<N, MyT>::GetLaserFrequencies() const
    {
        TDynamicColVector<double> frequencies(GetLaserCount());
        for (size_t i = 0; i < frequencies.Size(); i++)
        {
            auto lvls = GetLaserLevels(i);
            frequencies[i] = abs(GetLevel(lvls.first) - GetLevel(lvls.second));
        }
        return frequencies;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::AddLaser(
        const std::string& name, std::size_t lvl1, 
        std::size_t lvl2, double intensity, bool counter)
    {
        if (GetLaserIdxByName(name) < m_couplingLasers.size())
            return false; // coupling laser with the given name already exists

        if (lvl1 >= N  || lvl2 >= N || lvl1 == lvl2)
            return false; // invalid levels

        m_couplingLasers.emplace_back(name, lvl1, lvl2, 0.0, counter);
        SetLaserIntensity(name, intensity);

        if (!(~(*this)).OnLaserAdded(lvl1, lvl2, counter))
        {
            m_couplingLasers.pop_back();
            return false;
        }

        return true;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::AddLaserByName(
        const std::string& name, const std::string& lvl1, 
        const std::string& lvl2, double intensity, bool counter)
    {
        return AddLaser(name, GetLevelIndexByName(lvl1), 
            GetLevelIndexByName(lvl2), intensity, counter);
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::RemoveLaser(const std::string& name)
    {
        auto idx = GetLaserIdxByName(name);
        if (idx >= m_couplingLasers.size())
            return false;

        m_couplingLasers.erase(m_couplingLasers.begin() + idx);

        (~(*this)).OnLaserRemoved();

        return true;
    }

    template<std::size_t N, typename MyT>
    bool TNLevelSystemCRTP<N, MyT>::SetLaserIntensity(
        const std::string& name, double intensity)
    {
        std::size_t idx = GetLaserIdxByName(name);
        if (idx >= m_couplingLasers.size())
            return false;

        double conv = 2 / (SpeedOfLight_v * VacuumPermittivity_v);
        double electricField = std::sqrt(conv * std::abs(intensity));
        std::get<3>(m_couplingLasers[idx]) = GetElectricFieldFromIntensity(intensity);

        return true;
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetDopplerWidth(std::size_t from, std::size_t to) const
    {
        double frequency = std::abs(GetLevel(to) - GetLevel(from));
        return m_doppler.GetDopplerWidth(frequency);
    }

    template<std::size_t N, typename MyT>
    double TNLevelSystemCRTP<N, MyT>::GetDopplerWidthByName(const std::string& from, const std::string& to) const
    {
        return GetDopplerWidth(GetLevelIndexByName(from), GetLevelIndexByName(to));
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::CreateGroundState() const
    { 
        return CreateThermalStateHelper(0.0);
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::CreateThermalState() const
    { 
        return CreateThermalStateHelper(GetTemperature()); 
    }

    template<std::size_t N, typename MyT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::CreateThermalStateHelper(double temperature) const
    { 
        TStaticDensityMatrix<N> res(m_levelNames);

        TStaticColVector<double, N> lvlOffset;
        std::size_t gsIdx = std::min_element(m_levels.Data(), m_levels.Data() + N) - m_levels.Data();
        for (size_t i = 0; i < N; i++)
            lvlOffset(i) = m_levels(i) - m_levels(gsIdx);

        TStaticColVector<double, N> thermalWeights;
        if (temperature <= 0.0)
        {
            for (size_t i = 0; i < N; i++)
            {
                if (lvlOffset(i) == 0.0)
                    thermalWeights(i) = 1.0;
            }
        }
        else
        {
            double tmp = PlanckConstant_v / (BoltzmannConstant_v * temperature);
            for (size_t i = 0; i < N; i++)
                thermalWeights(i) = std::exp(-tmp * lvlOffset(i));
        }

        double norm = 0.0;
        for (size_t i = 0; i < N; i++)
            norm += thermalWeights(i);
        norm = 1.0 / norm;
        thermalWeights *= norm;

        for (size_t i = 0; i < N; i++)
            res(i, i) = thermalWeights(i);
        return res;
    }

    template<std::size_t N, typename MyT>
    std::array<double, N> TNLevelSystemCRTP<N, MyT>::GenerateDefaultLevels() const
    {
        std::array<double, N> lvls;
        std::fill(lvls.begin(), lvls.end(), 0.0);
        return lvls;
    }

    template<std::size_t N, typename MyT>
    std::array<std::string, N> TNLevelSystemCRTP<N, MyT>::GenerateDefaultLevelNames() const
    {
        std::array<std::string, N> names;
        for (size_t i = 0; i < N; i++)
            names[i] = std::to_string(i);
        return names;
    }

    template<std::size_t N, typename MyT>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::GetNaturalDensityMatrix(
        const TColVector<VT>& detunings,
        const TStaticDensityMatrix<N>& initial, 
        double velocity,
        double t0, double t, double dt)
    {
        // calculate appropriate amount of steps to obtain a dt near to the required dt
        std::size_t steps = static_cast<std::size_t>(std::ceil((t-t0) / dt));
        dt = (t-t0) / steps;
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        return this->EvolveNaturalDensityMatrix(
            auxData, initial, t0, dt, steps);
    }

    template<std::size_t N, typename MyT>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::GetNaturalDensityMatrixAv(
        const TColVector<VT>& detunings,
        const TStaticDensityMatrix<N>& initial, 
        double velocity,
        double t0, double t, double tav, double dt)
    {
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        
        // validate averaging time and generate starting and 
        // end time of the averaging process
        tav = tav > t - t0 ? t - t0 : tav;
        double t1 = t - tav / 2;
        double t2 = t1 + tav;

        // evolve up to the starting point of the averaging
        auto rho = this->GetNaturalDensityMatrix(
            detunings, initial, velocity, t0, t1, dt);

        std::size_t steps = static_cast<std::size_t>(std::ceil((t2-t1) / dt));
        dt = (t2 - t1) / steps;

        auto rhoAv = rho;
        for (std::size_t i = 0; i < steps; i++)
        {
            rho = this->EvolveNaturalDensityMatrix(
                auxData, rho, t1 + i*dt, dt, 1);
            rhoAv += rho;
        }
        rhoAv *= 1.0 / (steps + 1);

        return rhoAv;
    }
    
    template<std::size_t N, typename MyT>
    template<typename VT>
    std::pair<TDynamicColVector<double>, std::vector<TStaticDensityMatrix<N>>> TNLevelSystemCRTP<N, MyT>::GetNaturalTrajectory(
        const TColVector<VT>& detunings, const TStaticDensityMatrix<N>& initial, 
        double velocity, double t0, double t, double dt)
    {
        std::size_t steps = static_cast<std::size_t>(std::ceil((t-t0) / dt));

        std::vector<TStaticDensityMatrix<N>> trajectory;
        trajectory.reserve(steps + 1);
        trajectory.push_back(initial);
        
        auto rho = initial;
        const auto auxData = (~(*this)).GetHamiltonianAux(detunings, velocity);
        for (std::size_t i = 0; i < steps; i++)
        {
            rho = this->EvolveNaturalDensityMatrix(
                auxData, rho, t0 + i*dt, dt, 1);
            trajectory.push_back(rho);
        }

        return {QSim::CreateLinspaceCol(t0, t0+steps*dt, steps + 1), trajectory};
    }

    template<std::size_t N, typename MyT>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::GetDensityMatrix(
        const TColVector<VT>& detunings,
        const TStaticDensityMatrix<N>& initial, 
        double t0, double t, double dt)
    {
        auto func = [&](double vel)
        {
            return this->GetNaturalDensityMatrix(
                detunings, initial, vel, t0, t, dt).GetMatrix();
        };
        auto rho = this->m_doppler.Integrate(func);
        return TStaticDensityMatrix<N>(this->m_levelNames, rho);
    }

    template<std::size_t N, typename MyT>
    template<typename VT>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::GetDensityMatrixAv(
        const TColVector<VT>& detunings,
        const TStaticDensityMatrix<N>& initial, 
        double t0, double t, double tav, double dt)
    {
        auto func = [&](double vel)
        {
            return this->GetNaturalDensityMatrixAv(
                detunings, initial, vel, t0, t, tav, dt).GetMatrix();
        };
        auto rho = this->m_doppler.Integrate(func);
        return TStaticDensityMatrix<N>(this->m_levelNames, rho);
    }

    template<std::size_t N, typename MyT>
    template<typename AuxType>
    TStaticMatrix<std::complex<double>, N, N> TNLevelSystemCRTP<N, MyT>::GetDensityOpDerivative(
        const AuxType& auxData, const TStaticMatrix<std::complex<double>, N, N>& rho, double t) const
    {
        // von Neumann term
        auto h = (~(*this)).GetHamiltonianFast(auxData, t);
        TStaticMatrix<std::complex<double>, N, N> rhoPrime = -1.0i * (h * rho - rho * h);

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

    template<std::size_t N, typename MyT>
    template<typename AuxType>
    TStaticDensityMatrix<N> TNLevelSystemCRTP<N, MyT>::EvolveNaturalDensityMatrix(
        const AuxType& auxData, const TStaticDensityMatrix<N>& initial, 
        double t0, double dt, std::size_t steps)
    {
        using YType = TStaticMatrix<std::complex<double>, N, N>;
        RK4Integrator<double, YType> integrator;
        YType rho = initial;

        for (std::size_t i = 0; i < steps; i++)
        {
            auto func = [&](double x, const YType& y) 
            { 
                return this->GetDensityOpDerivative(auxData, y, x);
            };
            rho += integrator.Step(rho, t0 + i * dt, dt, func);       
        }

        return TStaticDensityMatrix<N>(this->m_levelNames, rho);
    }

}

// Include specific implementations
#include "NLevelSystemSC.h"
#include "NLevelSystemQM.h"

#endif
