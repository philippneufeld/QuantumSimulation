// Philipp Neufeld, 2021

#ifndef QSim_Math_Quadrature_H
#define QSim_Math_Quadrature_H

#include <cstdint>
#include <random> // for monte carlo

#include "../Util/ConstList.h"
#include "../Util/ConstexprFor.h"

namespace QSim
{
    // Midpoint integrator
    template<typename XTy>
    class TQuadMidpoint
    {
    public:
        template<typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
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

    template<typename XTy, int Denominator, int W1, int W2, int... Ws>
    class TQuadNewtonCotes
    {
    public:
        template<typename Func>
        auto Integrate(const Func& func, XTy a, XTy b, std::size_t n)
        {
            using YTy = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using WeightList = TConstList<int, W1, W2, Ws...>;
            constexpr std::size_t wcnt = TConstListSizeof_v<WeightList>;

            // adjust n to fit the method used
            if (n < wcnt)
                n = wcnt;
            n -= (n - 1) % (wcnt - 1);

            XTy dist = b - a;
            XTy dx = dist / (n - 1);

            // handle borders
            constexpr auto wFirst = TConstListFront_v<WeightList>;
            constexpr auto wLast = TConstListBack_v<WeightList>;
            YTy result = wFirst * func(a) +  wLast * func(b);

            using CLHelper = TConstListPopFront_t<TConstListPopBack_t<WeightList>>;
            using CL = TConstListAppend_t<int, wFirst + wLast, CLHelper>;
            constexpr std::size_t ccnt = TConstListSizeof_v<CL>;

            auto helper = [=, &result, &func](auto i)
            {
                constexpr std::size_t j0 = i + 1;
                YTy tmp = func(a + j0*dx);
                for (std::size_t j = j0 + ccnt; j < n - 1; j+=ccnt)
                    tmp += func(a + j*dx);
                result += TConstListGet_v<i, CL> * tmp;
            };
            ConstexprFor<std::size_t, 0, ccnt, 1>(helper);

            constexpr XTy den = static_cast<XTy>(Denominator) / ccnt;
            return result * (dx / den);
        }
    };

    // see https://de.wikipedia.org/wiki/Newton-Cotes-Formeln#Abgeschlossene_Newton-Cotes-Formeln
    template<typename XTy>
    using TQuadTrapezoidal = TQuadNewtonCotes<XTy, 2, 1, 1>;
    template<typename XTy>
    using TQuadSimpson = TQuadNewtonCotes<XTy, 6, 1, 4, 1>;
    template<typename XTy>
    using TQuadSimpson38 = TQuadNewtonCotes<XTy, 8, 1, 3, 3, 1>;
    template<typename XTy>
    using TQuadBoole = TQuadNewtonCotes<XTy, 90, 7, 32, 12, 32, 7>;
    template<typename XTy>
    using TQuadNC6Point = TQuadNewtonCotes<XTy, 288, 19, 75, 50, 50, 75, 19>;
    template<typename XTy>
    using TQuadWeddle = TQuadNewtonCotes<XTy, 840, 41, 216, 27, 272, 27, 216, 41>;
    template<typename XTy>
    using TQuadNC8Point = TQuadNewtonCotes<XTy, 17280, 751, 3577, 1323, 2989, 2989, 1323, 3577, 751>;

    // Alternative Simpson integrator that is suitable for Narrow peaks
    // https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson%27s_rule
    template<typename XTy>
    class TQuadSimpsonAlt
    {
    public:
        template<typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            // only n >= 6 allowed
            if (n < 5)
                return TQuadMidpoint<XTy>{}.Integrate(func, x0, x1, n);     

            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
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

    // Monte carlo integrator (uniform distribution)
    template<typename XTy>
    class TQuadMC
    {
    public:
        template<typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t n)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<XTy> dist(x0, x1);
            
            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
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
