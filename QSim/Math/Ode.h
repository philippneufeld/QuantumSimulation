// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Ode_H_
#define QSim_Math_Ode_H_

#include "../Platform.h"

namespace QSim
{
    // ODEStepper policy
    // Implements at least the following functions:
    // 1) Calculate evolution step of func from x to x+dx:
    //      auto Step(func, y, x, dx)

    // Euler ode integrator
    class ODEEulerPolicy
    {
    protected:
        ~ODEEulerPolicy() = default;

    public:
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        static auto Step(Func&& func, YTy&& y, XTy&& x, DXTy&& dx)
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
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        static auto Step(Func&& func, YTy&& y, XTy&& x, DXTy&& dx)
        {
            auto dx2 = dx/2;
            auto k1 = func(x, y);
            auto k2 = func(x + dx2, y + dx2*k1);
            auto k3 = func(x + dx2, y + dx2*k2);
            auto k4 = func(x + dx, y + dx*k3);
            return dx/6 * (k1 + 2*k2 + 2*k3 + k4);
        }
    };

    
    // Runge-Kutta 6th order ode integrator
    class ODERK6Policy
    {
    protected:
        ~ODERK6Policy() = default;

    public:
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        static auto Step(Func&& func, YTy&& y, XTy&& x, DXTy&& dx)
        {
            auto k1 = func(x, y);
            auto k2 = func(x + dx/3, y + (dx/3)*k1);
            auto k3 = func(x + 2*dx/3, y + (2*dx/3)*k2);
            auto k4 = func(x + dx/3, y + (dx/12)*k1 + (dx/3)*k2 - (dx/12)*k3);
            auto k5 = func(x + 5*dx/6, y + (25*dx/48)*k1 - (55*dx/24)*k2 + (35*dx/48)*k3 + (15*dx/8)*k4);
            auto k6 = func(x + dx/6, y + (3*dx/20)*k1 - (11*dx/20)*k2 - (1*dx/8)*k3 + (dx/2)*k4 + (dx/10)*k5);
            auto k7 = func(x + dx, y - (261*dx/260)*k1 + (33*dx/13)*k2 + (43*dx/156)*k3 - (118*dx/39)*k4 + (32*dx/195)*k5 + (80*dx/39)*k6);
            return (13*dx/200)*k1 + (11*dx/40)*k3 + (11*dx/40)*k4 + (4*dx/25)*k5 + (4*dx/25)*k6 + (13*dx/200)*k7;
        }
    };

    // Bogacki–Shampine ode integrator
    class ODEBS32Policy
    {
    protected:
        ~ODEBS32Policy() = default;

    public:      
        
        template<typename Func, typename YTy, typename XTy, typename DXTy>
        static auto Step(Func&& func, YTy&& y, XTy&& x, DXTy&& dx)
        {
            return StepWithErrorEst(func, y, x, dx).first;
        }

       template<typename Func, typename YTy, typename XTy, typename DXTy>
        static auto StepWithErrorEst(Func&& func, YTy&& y, XTy&& x, DXTy&& dx)
        {
            auto k1 = func(x, y);
            auto k2 = func(x + (dx/2), y + (dx/2)*k1);
            auto k3 = func(x + (3*dx/4), y + (3*dx/4)*k2);
            auto dy1 = (2*dx/9)*k1 + (1*dx/9)*k2 + (4*dx/9)*k3;
            auto k4 = func(x + dx, y + dy1);
            auto dy2 = (7*dx/24)*k1 + (dx/4)*k2 + (dx/3)*k3 + (dx/8)*k4;  
            auto err = dy2 - dy1;
            return std::make_pair(dy2, err);
        }
    };

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename ODEStepperPolicy=ODERK4Policy>
    class TODEStepper : public ODEStepperPolicy {};

}

#endif
