// Philipp Neufeld, 2021-2022

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
        auto Integrate(const Func& func, XTy a, XTy b, std::size_t n)
        {
            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
            if (n == 0)
                return YTy{};
            
            XTy dx = (b - a) / n;

            // use first iteration for initialization of the result (no addition needed)
            a = a + dx/2;
            YTy result = func(a);

            // execute integration
            for (std::size_t i = 1; i < n; i++)
                result += func(a + i*dx);
            
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
        auto Integrate(const Func& func, XTy a, XTy b, std::size_t n)
        {
            // only n >= 6 allowed
            if (n < 5)
                return TQuadMidpoint<XTy>{}.Integrate(func, a, b, n);     

            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
            XTy dx = (b - a) / (n - 1);
            
            // use first iteration for initialization of the result (no addition needed)
            YTy result = (func(a) + func(b)) * (9.0 / 24.0);
            result += (func(a + dx) + func(b - dx)) * (28.0 / 24.0);
            result += (func(a + 2*dx) + func(b - 2*dx)) * (23.0 / 24.0);

            // execute first part of the integration
            for (std::size_t i = 3; i < n - 3; i++)
                result += func(a + i*dx);

            return result * dx;
        }
    };

    // Adaptive step size integrator
    template<typename XTy>
    class TQuadAdaptive
    {
        template<typename Func>
        using YTy = decltype(std::declval<Func>()(std::declval<XTy>()) * std::declval<XTy>());
    public:

        template<typename Func>
        std::pair<YTy<Func>, std::size_t> IntegrateFevs(const Func& func, XTy a, XTy b, 
            std::size_t n, YTy<Func> rtol, YTy<Func> atol, std::size_t depth)
        {
            // adjust n to fit the simpson method used
            if (n < 5)
                n = 5;
            n -= (n - 1) % 4;
            std::size_t sections = (n - 1) / 4;

            
            XTy dist = b - a;
            XTy dx = dist / sections;
            atol /= sections;

            // first iteration
            YTy<Func> f0 = func(a);
            YTy<Func> f1 = func(a + 0.5*dx);
            YTy<Func> f2 = func(a + dx);
            auto res = IntegrateHelper(func, a, dx, 0, f0, f1, f2, rtol, atol, depth);
            YTy<Func> result = res.first;
            std::size_t fevs = 3 + res.second;

            // rest of the iterations
            for (size_t i = 1; i < sections; i++)
            {
                f0 = f2;
                f1 = func(a + (i+0.5)*dx);
                f2 = func(a + (i+1)*dx);
                res = IntegrateHelper(func, a+i*dx, dx, 0, f0, f1, f2, rtol, atol, depth);

                result += res.first;
                fevs += 2 + res.second;
            }
            
            return std::make_pair(result, fevs);
        }

        template<typename Func>
        auto Integrate(const Func& func, XTy a, XTy b, 
            std::size_t n, YTy<Func> rtol, YTy<Func> atol, std::size_t depth)
        {
            return IntegrateFevs(func, a, b, n, rtol, atol, depth).first;
        }

    private:
        template<typename Func>
        std::pair<YTy<Func>, std::size_t> IntegrateHelper(const Func& func, XTy a, XTy dx, std::size_t d,
            YTy<Func> f0, YTy<Func> f2, YTy<Func> f4, YTy<Func> rtol, YTy<Func> atol, std::size_t depth)
        {
            std::size_t fevs = 2;

            // calculate area via Simpson's rule
            YTy<Func> I1 = (dx / 6) * (f0 + 4*f2 + f4);
            
            // calculate function value at intermediary points
            YTy<Func> f1 = func(a + 0.25*dx);
            YTy<Func> f3 = func(a + 0.75*dx);

            // calculate area via two Simpson's rule steps with half step size
            YTy<Func> I2 = (dx / 12) * (f0 + 4*(f1 + f3) + 2*f2 + f4);
            
            YTy<Func> Ierr = (I2 - I1) / 15;
            YTy<Func> I = I2 + Ierr; // add error to I2 (I is equivalent to Boole's rule)

            // check if accuracy is sufficient
            if (std::abs(Ierr) > rtol*std::abs(I) + atol && d < depth)
            {
                XTy dx2 = dx / 2;
                auto res1 = IntegrateHelper(func, a, dx2, d+1, f0, f1, f2, rtol, atol / 2, depth);
                auto res2 = IntegrateHelper(func, a+dx2, dx2, d+1, f2, f3, f4, rtol, atol / 2, depth);
                
                I = res1.first + res2.first;
                fevs += res1.second + res2.second;
            }

            return std::make_pair(I, fevs);
        }

    };


    // Monte carlo integrator (uniform distribution)
    template<typename XTy>
    class TQuadMC
    {
    public:
        template<typename Func>
        auto Integrate(const Func& func, XTy a, XTy b, std::size_t n)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<XTy> dist(a, b);
            
            using FuncRet = decltype(func(std::declval<XTy>()) * std::declval<XTy>()); 
            using YTy = std::conditional_t<std::is_integral<FuncRet>::value, double, FuncRet>;
            if (n == 0)
                return YTy{};
            
            YTy result = func(dist(gen));
            for (std::size_t i = 1; i < n; i++)
                result += func(dist(gen));        

            result *= (b - a) / n;
            return result;
        }
    };

}

#endif
