// Philipp Neufeld, 2021-2022

#ifndef QSim_NLevelSystem_H_
#define QSim_NLevelSystem_H_

#include <complex>
#include <array>
#include <vector>
#include <map>
#include <set>
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

    template<int N, bool AM=false>
    class TNLevelSystem
    {
        using MyT = TNLevelSystem<N, AM>;
    public:

        using MatOp_t = Eigen::Matrix<std::complex<double>, N, N>;
        using Rabi_t = std::conditional_t<AM, std::function<double(double)>, double>;

        template<int compileTimeSize>
        using EnableDefaultCtor_t = std::enable_if_t<(compileTimeSize != DynamicDim_v)>;

        //template<int dummy=N, typename=EnableDefaultCtor_t<dummy>>
        TNLevelSystem() : TNLevelSystem(N) {}
        TNLevelSystem(unsigned int dims);

        // copy operations
        TNLevelSystem(const TNLevelSystem&) = default;
        TNLevelSystem(TNLevelSystem&&) = default;
        TNLevelSystem& operator=(const TNLevelSystem&) = default;
        TNLevelSystem& operator=(TNLevelSystem&&) = default;

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

        // level coupulings
        Rabi_t& GetCoupling(unsigned int cidx) { return std::get<2>(m_couplings[cidx]); }
        const Rabi_t& GetCoupling(unsigned int cidx) const { return std::get<2>(m_couplings[cidx]); }
        bool AddCoupling(unsigned int l1, unsigned int l2, Rabi_t rabi);
        void ClearCouplings();
        Eigen::VectorXd GetCouplingResonanceFreqs() const;

        // hamiltonian
        template<bool dummy=AM>
        std::enable_if_t<!dummy, MatOp_t> GetHamiltonian(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;
        template<bool dummy=AM>
        std::enable_if_t<!dummy, MatOp_t> GetHamiltonian(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const;
        template<bool dummy=AM>
        std::enable_if_t<dummy, MatOp_t> GetHamiltonian(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const;

        // Steady state
        template<bool dummy=AM>
        std::enable_if_t<!dummy, MatOp_t> GetDensityMatrixSS(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const;

        // time evolution of density matrix
        MatOp_t GetDensityMatrix(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const MatOp_t& rho0,
            double t0, double t, double dt);

        MatOp_t GetDensityMatrixAv(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const MatOp_t& rho0,
            double t0, double t, double tav, double dt);

        Eigen::VectorXd GetTrajectoryTimeaxis(double t0, double dt, std::size_t n);
        std::vector<MatOp_t> GetTrajectory(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const MatOp_t& rho0,
            double t0, double t, double dt);
        
        // create specific density matrices
        MatOp_t CreateGroundState() const;
        MatOp_t CreateThermalState(double temperature) const;
       
        // get temporal derivative of the density operator
        MatOp_t GetDensityOpDerivative(
            const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
            const MatOp_t& rho, double t) const;

    private:
        //
        // Helper methods to prepare the hamiltonian calculation
        //

        bool GeneratePhotonBasis(std::vector<unsigned int>& trans_path, 
            std::set<unsigned int>& visitedLevels, unsigned int transFrom);
        
        template<bool dummy=AM>
        std::enable_if_t<!dummy, void> GeneratePartialHamiltonian();
        template<bool dummy=AM>
        std::enable_if_t<dummy, void> GeneratePartialHamiltonian();

        bool UpdateAuxilliaries();

    private:
        // Properties of the system
        Eigen::Matrix<double, N, 1> m_levels;
        std::map<std::pair<unsigned int, unsigned int>, double> m_decays;
        std::vector<std::tuple<unsigned int, unsigned int, Rabi_t>> m_couplings;

        // auxilliaries
        Eigen::Matrix<double, N, Eigen::Dynamic> m_photonBasis;
        MatOp_t m_partialHamiltonian;
    };

    template<int N, bool AM>
    using TNLevelSystemMatOp_t = typename TNLevelSystem<N, AM>::MatOp_t;


    template<int N, bool AM>
    TNLevelSystem<N, AM>::TNLevelSystem(unsigned int dims) 
        : m_levels(dims)
    {
        m_levels.setZero();
    }

    template<int N, bool AM>
    double TNLevelSystem<N, AM>::GetLevel(unsigned int idx) const
    {
        return idx < GetDims() ? m_levels(idx) : 0.0;
    }

    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::SetLevel(unsigned int idx, double level)
    {
        if (idx >= GetDims())
            return false;
        m_levels[idx] = level;
        return true;
    }

    template<int N, bool AM>
    double TNLevelSystem<N, AM>::GetDecay(
        unsigned int from, unsigned int to) const
    {
        auto it = m_decays.find(std::make_pair(from, to));
        return it != m_decays.end() ? it->second : 0.0;
    }

    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::SetDecay(
        unsigned int from, unsigned int to, double rate)
    {
        if (from >= GetDims() || to >= GetDims())
            return false; // index out of bound
        m_decays[std::make_pair(from, to)] = rate;
        return true;
    }

    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::AddDecay(
        unsigned int from, unsigned int to, double rate)
    {
        return SetDecay(from, to, GetDecay(from, to) + rate);
    }

    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::AddCoupling(unsigned int l1, unsigned int l2, Rabi_t rabi)
    {
        m_couplings.emplace_back(l1, l2, rabi);
        
        if (!UpdateAuxilliaries())
        {
            m_couplings.pop_back();
            UpdateAuxilliaries();
            return false;
        }

        return true;
    }

    template<int N, bool AM>
    void TNLevelSystem<N, AM>::ClearCouplings()
    {
        m_couplings.clear();
        UpdateAuxilliaries();
    }

    template<int N, bool AM>
    Eigen::VectorXd TNLevelSystem<N, AM>::GetCouplingResonanceFreqs() const
    {
        Eigen::VectorXd freqs(m_couplings.size());
        for (int i=0; i<m_couplings.size(); i++)
        {
            auto [l1, l2, rabi] = m_couplings[i];
            freqs(i) = std::abs(this->GetLevel(l2) - this->GetLevel(l1));
        }
        return freqs;
    }

    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<!dummy, typename TNLevelSystem<N, AM>::MatOp_t> TNLevelSystem<N, AM>::GetHamiltonian(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {        
        // add light field hamiltonian
        MatOp_t hamiltonian = m_partialHamiltonian;
        hamiltonian.diagonal() += TwoPi_v * (m_photonBasis * laserFreqs);
        return hamiltonian;
    }

    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<!dummy, typename TNLevelSystem<N, AM>::MatOp_t> TNLevelSystem<N, AM>::GetHamiltonian(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const
    {
        // system is time independent -> drop t
        return GetHamiltonian(laserFreqs);
    }

    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<dummy, typename TNLevelSystem<N, AM>::MatOp_t> TNLevelSystem<N, AM>::GetHamiltonian(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, double t) const
    {
        MatOp_t hamiltonian = m_partialHamiltonian;
        
        // add light field hamiltonian
        hamiltonian.diagonal() += TwoPi_v * (m_photonBasis * laserFreqs);

        // add time-dependent coupling term
        for (auto [l1, l2, rabi]: m_couplings)
        {
            double rabiT = std::invoke(rabi, t);
            hamiltonian(l1, l2) += Pi_v * rabiT;
            hamiltonian(l2, l1) += Pi_v * std::conj(rabiT);
        }

        return hamiltonian;
    }

    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<!dummy, typename TNLevelSystem<N, AM>::MatOp_t> TNLevelSystem<N, AM>::GetDensityMatrixSS(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs) const
    {
        auto h = GetHamiltonian(laserFreqs);
        
        using ATy = std::conditional_t<MyT::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, N*N>, Eigen::MatrixXcd>;
        using BTy = std::conditional_t<MyT::IsStaticDim(),
            Eigen::Matrix<std::complex<double>, N*N + 1, 1>, Eigen::VectorXcd>;
        using SSTy = MatOp_t;
        
        unsigned int dims = this->GetDims();
        ATy A = ATy::Zero(dims*dims + 1, dims*dims);

        // von Neumann part of the evolution operator
        for (unsigned int i = 0; i < dims; i++)
        {
            for (unsigned int j = 0; j < dims; j++)
            {
                for (unsigned int k = 0; k < dims; k++)
                {
                    A(i*dims+k, j*dims+k) -= 1.0i * h(i, j);
                    A(k*dims+i, k*dims+j) += 1.0i * h(j, i);
                }
            }
        }
        
        // Lindblad part of the evolution operator
        for(const auto& decay: this->GetDecays())
        {
            auto [lvls, rate] = decay;
            auto [i, f] = lvls;
            
            double ratePi = Pi_v * rate;
            A(f*dims+f, i*dims+i) += 2 * ratePi;
            for (unsigned int j = 0; j < dims; j++)
            {
                A(j*dims+i, j*dims+i) -= ratePi;
                A(i*dims+j, i*dims+j) -= ratePi;
            }
        }

        // Normalization condition
        for (unsigned int i = 0; i < dims; i++)
            A(dims*dims, i*dims+i) += 1.0;

        // bring all columns to the same order of magnitude by
        // dividing every row by its abs-max value
        auto Ascale = (1 / A.array().cwiseAbs().rowwise().maxCoeff()).eval();
        A = Ascale.matrix().asDiagonal() * A;

        // generate right hand side of the equation A*x=b
        // b = (0, 0, ..., 0, 1)
        auto b = BTy::Unit(dims*dims + 1, dims*dims);
        
        // solve and subsequently reshape from vetor to matrix
        auto x = ((A.adjoint()*A).ldlt().solve(A.adjoint()*b)).eval();
        return Eigen::Map<SSTy>(x.data(), dims, dims).transpose();
    }
    
    template<int N, bool AM>
    typename TNLevelSystem<N, AM>::MatOp_t TNLevelSystem<N, AM>::GetDensityMatrix(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const MatOp_t& rho0,
        double t0, double t, double dt)
    {
        // calculate appropriate amount of steps to obtain a dt near to the required dt
        unsigned int steps = static_cast<unsigned int>(std::ceil((t-t0) / dt));
        dt = (t-t0) / steps;

        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54CKPolicy> integrator;
        using YType = MatOp_t;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(laserFreqs, y, x); };

        return integrator.IntegrateTo(func, rho0, t0, t, dt);
    }

    template<int N, bool AM>
    typename TNLevelSystem<N, AM>::MatOp_t TNLevelSystem<N, AM>::GetDensityMatrixAv(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const MatOp_t& rho0,
        double t0, double t, double tav, double dt)
    {
        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54CKPolicy> integrator;
        using YType = MatOp_t;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(laserFreqs, y, x); };
        
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
    
    template<int N, bool AM>
    Eigen::VectorXd TNLevelSystem<N, AM>::GetTrajectoryTimeaxis(double t0, double dt, std::size_t n)
    {
        if (n==0) return Eigen::VectorXd();
        return Eigen::VectorXd::LinSpaced(n, t0, t0 + (n-1)*dt);
    }

    template<int N, bool AM>
    std::vector<typename TNLevelSystem<N, AM>::MatOp_t> TNLevelSystem<N, AM>::GetTrajectory(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const MatOp_t& rho0,
        double t0, double t, double dt)
    {
        // allocate memory for the trajectory data
        std::size_t n = static_cast<std::size_t>(std::ceil((t-t0) / dt));
        std::vector<MatOp_t> trajectory;
        trajectory.reserve(n + 1);
        trajectory.push_back(rho0);

        // define integrator and function to be integrated
        TODEIntegrator<ODEAd54CKPolicy> integrator;
        using YType = MatOp_t;
        auto func = [&](double x, const YType& y) { return GetDensityOpDerivative(laserFreqs, y, x); };
        
        auto rho = rho0;
        double dtInner = dt;
        for (std::size_t i = 0; i < n; i++)
        {
            rho = integrator.IntegrateTo(func, rho, t0+i*dt, t0+(i+1)*dt, dtInner);
            trajectory.push_back(rho);
        };

        return trajectory;
    }

    template<int N, bool AM>
    typename TNLevelSystem<N, AM>::MatOp_t TNLevelSystem<N, AM>::CreateGroundState() const
    { 
        return CreateThermalState(0.0);
    }

    template<int N, bool AM>
    typename TNLevelSystem<N, AM>::MatOp_t TNLevelSystem<N, AM>::CreateThermalState(double temperature) const
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

    template<int N, bool AM>
    typename TNLevelSystem<N, AM>::MatOp_t TNLevelSystem<N, AM>::GetDensityOpDerivative(
        const Eigen::Ref<const Eigen::VectorXd>& laserFreqs, 
        const MatOp_t& rho, 
        double t) const
    {
        // von Neumann term
        auto h = this->GetHamiltonian(laserFreqs, t);
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

    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::GeneratePhotonBasis(std::vector<unsigned int>& trans_path, 
        std::set<unsigned int>& visitedLvls, unsigned int transFrom)
    {
        // m_photonBasis(i, j) contains the (relative) photon number
        // of the j-th laser field in the i-th atomic state

        // check for circular transition path
        if (visitedLvls.count(transFrom) != 0)
            return false;
        
        visitedLvls.insert(transFrom);

        // Iterate over all transitions
        for (unsigned int tIdx = 0; tIdx < m_couplings.size(); tIdx++)
        {
            // skip if the transition was already in the current transition path
            if (std::find(trans_path.begin(), trans_path.end(), tIdx) != trans_path.end())
                continue;

            // skip if currLevel is not involved in the current transition
            // otherwise check to which level the transition leads
            auto [l1, l2, rabi] = m_couplings[tIdx];
            unsigned int transTo = 0;
            if (l1 == transFrom)
                transTo = l2;
            else if (l2 == transFrom)
                transTo = l1;
            else
                continue;
            
            // check whether a photon is absorbed or emitted
            double relPhotonNumber = 1;
            if (this->GetLevel(transTo) > this->GetLevel(transFrom))
                relPhotonNumber = -1;

            // Add(Subtract) one photon to(from) the "from" state to obtain the "to" state
            m_photonBasis.row(transTo) = m_photonBasis.row(transFrom);
            m_photonBasis(transTo, tIdx) += relPhotonNumber;

            // continue recursively
            trans_path.push_back(tIdx);
            auto success = GeneratePhotonBasis(trans_path, visitedLvls, transTo);
            trans_path.pop_back();

            if (!success)
                return false;
        }

        return true;
    }

    // Atom and interaction hamiltonian
    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<!dummy, void> TNLevelSystem<N, AM>::GeneratePartialHamiltonian()
    {
        // Atom hamiltonian
        m_partialHamiltonian = TwoPi_v * (m_levels.array() - m_levels(0)).eval().matrix().asDiagonal();
        
        // Interaction hamiltonian
        for (auto [l1, l2, rabi]: m_couplings)
        {
            m_partialHamiltonian(l1, l2) += Pi_v * rabi;
            m_partialHamiltonian(l2, l1) += Pi_v * std::conj(rabi);
        }
    }

    template<int N, bool AM>
    template<bool dummy>
    std::enable_if_t<dummy, void> TNLevelSystem<N, AM>::GeneratePartialHamiltonian()
    {
        // Atom hamiltonian
        m_partialHamiltonian = TwoPi_v * (m_levels.array() - m_levels(0)).eval().matrix().asDiagonal();
    }


    template<int N, bool AM>
    bool TNLevelSystem<N, AM>::UpdateAuxilliaries()
    {
        // initialize photon basis matrix which can be used for the calculation
        // of the light field contribution to the hamiltonian
        m_photonBasis.setZero(this->GetDims(), m_couplings.size());

        // allocate memory for helper variables
        std::vector<unsigned int> trans_path;
        std::set<unsigned int> visited_levels;
        trans_path.reserve(m_couplings.size());

        while (visited_levels.size() < this->GetDims())
        {
            unsigned int head = 0;
            for(; visited_levels.count(head) != 0; head++);

            if (!GeneratePhotonBasis(trans_path, visited_levels, head))
            {
                m_photonBasis.resize(this->GetDims(), 0);
                return false;
            }
        }

        GeneratePartialHamiltonian();

        return true;
    }

}

#endif
