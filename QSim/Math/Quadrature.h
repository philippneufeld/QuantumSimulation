// Philipp Neufeld, 2021

#ifndef QSim_Math_Quadrature_H
#define QSim_Math_Quadrature_H

#include <cstdint>
#include <random> // for monte carlo

namespace QSim
{

    // Rectangle integrator
    class QuadMidpoint
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            using YTy = decltype(func(std::declval<XTy>()));
            if (n == 0)
                return YTy{};
            
            XTy dx = (x1 - x0) / n;

            // use first iteration for initialization of the result (no addition needed)
            x0 = x0 + dx/2;
            YTy result = func(x0);

            // execute integration
            for (std::size_t i = 1; i < n; i++)
                result += func(x0 + i*dx);
            
            return result * dx;
        }
    };


    // Trapezoidal integrator
    class QuadTrapezoidal
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only n >= 3 allowed
            if (n < 3)
                return QuadMidpoint{}.Integrate(func, x0, x1, n);
            
            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) / 2;
            
            // execute integration
            for (std::size_t i = 1; i < n - 1; i++)
                result += func(x0 + i*dx);
            
            return result * dx;
        }
    };

    // Simpson integrator
    class QuadSimpson
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only odd n >= 5 allowed
            if (n < 5)
                return QuadMidpoint{}.Integrate(func, x0, x1, n);     
            n -= (1 - n % 2);

            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) / 2;

            // execute first part of the integration
            YTy tmp = func(x0 + dx);
            for (std::size_t i = 3; i < n; i+=2)
                tmp += func(x0 + i*dx);
            result += 2*tmp;
            
            // execute second part of the integration
            for (std::size_t i = 2; i < n - 1; i+=2)
                result += func(x0 + i*dx);

            return result * (2 * dx / 3);
        }
    };

    // Alternative Simpson integrator that is suitable for Narrow peaks
    // https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson%27s_rule
    class QuadSimpsonAlt
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only n >= 6 allowed
            if (n < 5)
                return QuadMidpoint{}.Integrate(func, x0, x1, n);     

            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) * (9.0 / 24.0);
            result += (func(x0 + dx) + func(x1 - dx)) * (28.0 / 24.0);
            result += (func(x0 + 2*dx) + func(x1 - 2*dx)) * (23.0 / 24.0);

            // execute first part of the integration
            for (std::size_t i = 3; i < n - 3; i++)
                result += func(x0 + i*dx);

            return result * dx;
        }
    };

    // Boole integrator
    class QuadBoole
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only n = {5, 9, 13, 17, ...}
            if (n < 5)
               return QuadMidpoint{}.Integrate(func, x0, x1, n);
            n -= (n - 1) % 4;
            
            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) * 7;

            // execute first part of the integration
            YTy tmp = func(x0 + dx);
            for (std::size_t i = 3; i < n; i+= 2)
                tmp += func(x0 + i*dx);
            result += 32*tmp;

            tmp = func(x0 + 2*dx);
            for (std::size_t i = 6; i < n; i+= 4)
                tmp += func(x0 + i*dx);
            result += 12*tmp;

            if (n > 5)
            {
                YTy tmp = func(x0 + 4*dx);
                for (std::size_t i = 8; i <= n - 4; i+= 4)
                    tmp += func(x0 + i*dx);
                result += 14*tmp;
            }
            
            return result * (2 * dx / 45);
        }
    };

    // Weddle integrator
    class QuadWeddle
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only n = {7, 13, 19, 25, ...}
            if (n < 7)
               return QuadMidpoint{}.Integrate(func, x0, x1, n);
            n -= (n - 1) % 6;
            
            using YTy = decltype(func(std::declval<XTy>()));
            XTy dx = (x1 - x0) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(x0) + func(x1)) * 41;

            YTy tmp = func(x0 + dx);
            for (std::size_t i = 7; i < n; i+= 6)
                tmp += func(x0 + i*dx);
            result += 216*tmp;

            tmp = func(x0 + 2*dx);
            for (std::size_t i = 8; i < n; i+= 6)
                tmp += func(x0 + i*dx);
            result += 27*tmp;

            tmp = func(x0 + 3*dx);
            for (std::size_t i = 9; i < n; i+= 6)
                tmp += func(x0 + i*dx);
            result += 272*tmp;

            tmp = func(x0 + 4*dx);
            for (std::size_t i = 10; i < n; i+= 6)
                tmp += func(x0 + i*dx);
            result += 27*tmp;

            tmp = func(x0 + 5*dx);
            for (std::size_t i = 11; i < n; i+= 6)
                tmp += func(x0 + i*dx);
            result += 216*tmp;

            if (n > 7)
            {
                YTy tmp = func(x0 + 6*dx);
                for (std::size_t i = 12; i <= n - 6; i+= 6)
                    tmp += func(x0 + i*dx);
                result += 82*tmp;
            }
            
            return result * (dx / 140);
        }
    };


    // Monte carlo integrator (uniform distribution)
    class QuadMC
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<XTy> dist(x0, x1);
            
            using YTy = decltype(func(std::declval<XTy>()));
            if (n == 0)
                return YTy{};
            
            YTy result = func(dist(gen));
            for (std::size_t i = 1; i < n; i++)
                result += func(dist(gen));        

            result *= (x1 - x0) / n;
            return result;
        }
    };

}

#endif
