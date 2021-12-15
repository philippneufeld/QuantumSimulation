// Philipp Neufeld, 2021

#ifndef QSim_Math_Quadrature_H
#define QSim_Math_Quadrature_H

#include <cstdint>

namespace QSim
{

    // Rectangle integrator
    class QuadMidpoint
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t cntFev)
        {
            using YTy = decltype(func(std::declval<XTy>()));
            if (cntFev == 0)
                return YTy{};
            
            XTy dx = (x1 - x0) / cntFev;

            // use first iteration for initialization of the result (no addition needed)
            x0 = x0 + dx/2;
            YTy result = func(x0);

            // execute integration
            for (std::size_t i = 1; i < cntFev; i++)
                result += func(x0 + i*dx);
            
            return result * dx;
        }
    };


    // Trapezoidal integrator
    class QuadTrapezoidal
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t cntFev)
        {
            if (cntFev < 2)
                return QuadMidpoint{}.Integrate(func, x0, x1, cntFev);
            
            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / cntFev;
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) / 2;
            
            // execute integration
            for (std::size_t i = 1; i < cntFev - 1; i++)
                result += func(x0 + i*dx);
            
            return result * dx;
        }
    };

    // Simpson integrator
    class QuadSimpson
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t cntFev)
        {
            std::size_t n = cntFev / 2;
            if (n < 1)
                return QuadTrapezoidal{}.Integrate(func, x0, x1, cntFev);
            
            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / cntFev;
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) / 2;

            // execute first part of the integration
            YTy result1 = func(x0 + dx);
            for (std::size_t i = 1; i < n; i++)
                result1 += func(x0 + (2*i + 1)*dx);
            result += 2*result1;
            
            // execute second part of the integration
            for (std::size_t i = 1; i < n; i++)
                result += func(x0 + (2*i)*dx);

            return result * (dx / 6);
        }
    };

}

#endif
