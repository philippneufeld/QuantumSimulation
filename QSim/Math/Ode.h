// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Ode_H_
#define QSim_Math_Ode_H_

#include "../Platform.h"
#include "MatrixTraits.h"

namespace QSim
{
    // ODEStepper policy
    // Implements at least the following functions:
    // 1) Calculate evolution step of func from x to x+dx:
    //      auto Step(func, y, x, dx)
    // Additionally a static constexpr boolean member 
    // "IsAdaptive_v" must be present

    namespace Internal
    {
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        struct TODEResultType
        {
            using type = TMatrixEvalType_t<decltype(std::declval<DXTy>() *
                std::declval<std::invoke_result_t<Func, XTy, YTy>>())>;
        };
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        using TODEResultType_t = typename TODEResultType<Func, YTy, XTy, DXTy>::type;
    }

    // Euler ode integrator
    class ODEEulerPolicy
    {
    protected:
        ~ODEEulerPolicy() = default;

    public:
        constexpr static bool IsAdaptive_v = false;

        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            return func(x, y) * dx;
        }
    };

    // Runge-Kutta 4th order ode integrator
    class ODERK4Policy
    {
    protected:
        ~ODERK4Policy() = default;

    public:
        constexpr static bool IsAdaptive_v = false;

        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            auto dx2 = dx/2;
            auto k1 = func(x, y);
            auto k2 = func(x + dx2, y + dx2*k1);
            auto k3 = func(x + dx2, y + dx2*k2);
            auto k4 = func(x + dx, y + dx*k3);
            return dx/6 * (k1 + 2*k2 + 2*k3 + k4);
        }
    };

    // Heun-Euler adaptive integrator
    class ODEAd21HEPolicy
    {
    protected:
        ~ODEAd21HEPolicy() = default;

    public:
        constexpr static bool IsAdaptive_v = true;
        constexpr static std::size_t OrderLower_v = 1;
        constexpr static std::size_t OrderHigher_v = 2;
    
        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            return StepWithErrorEst(func, y, x, dx).first;
        }

        template<typename Func, typename YTy>
        static auto StepWithErrorEst(Func&& func, YTy&& y, double x, double dx)
        {
            auto k1 = func(x, y);
            auto k2 = func(x + dx, y + dx*k1);

            using RTy = TMatrixMulResultFP_t<YTy, double>;
            RTy dy1 = dx*k1;
            RTy dy2 = (dx/2) * (k1 + k2);
            RTy err = dy2 - dy1;

            return std::make_pair(dy2, err);
        }

    };
    
    // Bogacki–Shampine ode integrator
    // https://doi.org/10.1016/0893-9659(89)90079-7
    class ODEAd32BSPolicy
    {
    protected:
        ~ODEAd32BSPolicy() = default;

    public:
        constexpr static bool IsAdaptive_v = true;
        constexpr static std::size_t OrderLower_v = 2;
        constexpr static std::size_t OrderHigher_v = 3;
    
        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            return StepWithErrorEst(func, y, x, dx).first;
        }

        template<typename Func, typename YTy>
        static auto StepWithErrorEst(Func&& func, YTy&& y, double x, double dx)
        {
            using RTy = TMatrixMulResultFP_t<YTy, double>;
            
            auto k1 = func(x, y);
            auto k2 = func(x + (dx/2), y + (dx/2)*k1);
            auto k3 = func(x + (3*dx/4), y + (3*dx/4)*k2); 
            RTy dy1 = (2*dx/9)*k1 + (1*dx/3)*k2 + (4*dx/9)*k3;

            auto k4 = func(x + dx, y + dy1);
            RTy dy2 = (7*dx/24)*k1 + (dx/4)*k2 + (dx/3)*k3 + (dx/8)*k4;  
            RTy err = dy2 - dy1;

            return std::make_pair(dy2, err);
        }
    };

    // Cash-Karp ode integrator
    // https://doi.org/10.1145/79505.79507
    class ODEAd54CKPolicy
    {
    protected:
        ~ODEAd54CKPolicy() = default;

    public:
        constexpr static bool IsAdaptive_v = true;
        constexpr static std::size_t OrderLower_v = 4;
        constexpr static std::size_t OrderHigher_v = 5;
    
        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            return StepWithErrorEst(func, y, x, dx).first;
        }

        template<typename Func, typename YTy>
        static auto StepWithErrorEst(Func&& func, YTy&& y, double x, double dx)
        {
            using RTy = TMatrixMulResultFP_t<YTy, double>;
            
            auto k1 = func(x, y);
            auto k2 = func(x + (dx/5), y + (dx/5)*k1);
            auto k3 = func(x + (3*dx/10), y + (3*dx/40)*k1 + (9*dx/40)*k2);
            auto k4 = func(x + (3*dx/5), y + (3*dx/10)*k1 - (9*dx/10)*k2 + (6*dx/5)*k3); 
            auto k5 = func(x + dx, y - (11*dx/54)*k1 + (5*dx/2)*k2 - (70*dx/27)*k3 + (35*dx/27)*k4); 
            auto k6 = func(x + (7*dx/8), y + (1631*dx/55296)*k1 + (175*dx/512)*k2 + (575*dx/13824)*k3 + (44275*dx/110592)*k4 + (253*dx/4096)*k5); 
            RTy dy1 = (2825*dx/27648)*k1 + (18575*dx/48384)*k3 + (13525*dx/55296)*k4 + (277*dx/14336)*k5 + (dx/4)*k6;  
            RTy dy2 = (37*dx/378)*k1 + (250*dx/621)*k3 + (125*dx/594)*k4 + (512*dx/1771)*k6;
            RTy err = dy2 - dy1;

            return std::make_pair(dy2, err);
        }
    };

    // Dormand–Prince ode integrator
    // see https://doi.org/10.1016/0771-050X(80)90013-3
    class ODEAd54DPPolicy
    {
    protected:
        ~ODEAd54DPPolicy() = default;

    public:
        constexpr static bool IsAdaptive_v = true;
        constexpr static std::size_t OrderLower_v = 4;
        constexpr static std::size_t OrderHigher_v = 5;
    
        template<typename Func, typename YTy>
        static auto Step(Func&& func, YTy&& y, double x, double dx) -> 
            TMatrixMulResultFP_t<YTy, double>
        {
            return StepWithErrorEst(func, y, x, dx).first;
        }

        template<typename Func, typename YTy>
        static auto StepWithErrorEst(Func&& func, YTy&& y, double x, double dx)
        {
            using RTy = TMatrixMulResultFP_t<YTy, double>;
            
            auto k1 = func(x, y);
            auto k2 = func(x + (dx/5), y + (dx/5)*k1);
            auto k3 = func(x + (3*dx/10), y + (3*dx/40)*k1 + (9*dx/40)*k2); 
            auto k4 = func(x + (4*dx/5), y + (44*dx/45)*k1 - (56*dx/15)*k2 + (32*dx/9)*k3); 
            auto k5 = func(x + (8*dx/9), y + (19372*dx/6561)*k1 - (25360*dx/2187)*k2 + (64448*dx/6561)*k3 - (212*dx/729)*k4); 
            auto k6 = func(x + dx, y + (9017*dx/3168)*k1 - (355*dx/33)*k2 + (46732*dx/5247)*k3 + (49*dx/176)*k4 - (5103*dx/18656)*k5); 
            RTy dy1 = (35*dx/384)*k1 + (500*dx/1113)*k3 + (125*dx/192)*k4 - (2187*dx/6784)*k5 + (11*dx/84)*k6;

            auto k7 = func(x + dx, y + dy1); 
            RTy dy2 = (5179*dx/57600)*k1 + (7571*dx/16695)*k3 + (393*dx/640)*k4 - (92097*dx/339200)*k5 + (187*dx/2100)*k6 + (1*dx/40)*k6;  
            RTy err = dy2 - dy1;

            return std::make_pair(dy2, err);
        }
    };

    template<typename ODEStepperPolicy>
    class TODEStepperPolicyIsAdaptive : public std::integral_constant<bool, ODEStepperPolicy::IsAdaptive_v> {};
    template<typename ODEStepperPolicy>
    constexpr static bool TODEStepperPolicyIsAdaptive_v 
        = TODEStepperPolicyIsAdaptive<ODEStepperPolicy>::value;

    // 
    namespace Internal
    {
        // Non-adaptive implementation
        template<typename ODEStepperPolicy>
        class TODEIntegratorHelperFixed : public ODEStepperPolicy
        {
        public:
            template<typename Func, typename Y0Ty>
            auto IntegrateToFixed(Func&& func, Y0Ty&& y0, double x, double x1, double dx)
            {
                using YTy = TMatrixMulResultFP_t<Y0Ty, double>;
                YTy y = std::forward<Y0Ty>(y0);
                
                while (x < x1)
                {
                    double dxEff = std::min(x1 - x, dx);
                    x += dxEff;
                    y += ODEStepperPolicy::Step(func, y, x, dxEff);
                }
                
                return y;
            }
        };

        // Non-adaptive implementation
        template<typename ODEStepperPolicy, bool isAdaptive>
        class TODEIntegratorHelper 
            : public TODEIntegratorHelperFixed<ODEStepperPolicy>
        {
        public:
            template<typename Func, typename Y0Ty>
            auto IntegrateTo(Func&& func, Y0Ty&& y0, double x, double x1, double dx)
            {
                return TODEIntegratorHelperFixed<ODEStepperPolicy>::IntegrateToFixed(func, y0, x, x1, dx);
            }
        };

        // Adaptive implementation
        template<typename ODEStepperPolicy>
        class TODEIntegratorHelper<ODEStepperPolicy, true> 
            : public TODEIntegratorHelperFixed<ODEStepperPolicy>
        {
        public:
            TODEIntegratorHelper() : m_precision(1e-10), m_safetyFactor(0.95) {}

            void SetPrecision(double precision) { m_precision = precision; }
            void SetSafetyFactor(double safety) { m_safetyFactor = safety; }

            double SetPrecision(double precision) const { return m_precision; }
            double SetSafetyFactor(double safety) const { return m_safetyFactor; }

            template<typename Func, typename Y0Ty>
            auto IntegrateTo(Func&& func, Y0Ty&& y0, double x, double x1, double& dx)
            {
                using YTy = TMatrixMulResultFP_t<Y0Ty, double>;                
                using YTyAbsMax = TMatrixMulResultFP_t<decltype(TMatrixAbsMax<YTy>::Get(std::declval<YTy>())), double>;
                constexpr double exponentInc = 1.0 / (ODEStepperPolicy::OrderLower_v);
                constexpr double exponentDec = 1.0 / (ODEStepperPolicy::OrderLower_v + 1);

                YTy y = std::forward<Y0Ty>(y0);
                
                while (x < x1)
                {
                    double dxEff = std::min(x1 - x, dx);
                    auto [dy, err] = ODEStepperPolicy::StepWithErrorEst(func, y, x, dxEff);

                    YTyAbsMax maxErr = TMatrixAbsMax<YTy>::Get(err);
                    
                    double adjust = 1.0;
                    if (maxErr > m_precision)
                    {
                        // don't accept and decrease step size
                        adjust = std::max(0.1, m_safetyFactor * std::pow(m_precision / maxErr, exponentDec));    
                    }
                    else
                    {
                        // increase step size if needed
                        if (!(dxEff < dx)) // ignore remainder step
                            adjust = std::min(10.0, m_safetyFactor * std::pow(m_precision / maxErr, exponentInc));
                        
                        // accept
                        x += dxEff;
                        y += dy;
                    }

                    dx = dxEff * adjust;
                }

                return y;
            }

        private:
            double m_precision;
            double m_safetyFactor;
        };
    }

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename ODEStepperPolicy=ODERK4Policy>
    class TODEIntegrator : public Internal::TODEIntegratorHelper<
        ODEStepperPolicy, TODEStepperPolicyIsAdaptive_v<ODEStepperPolicy>> {};


}

#endif
