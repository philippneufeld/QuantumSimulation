// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Quad_H
#define QSim_Math_Quad_H

#include <cstdint>
#include <limits>

#include "MatrixTraits.h"
#include "../Util/ConstList.h"
#include "../Util/ConstexprFor.h"

namespace QSim
{

    // Quadrature policy
    // Implements at least the following functions:
    // 1) Integrate func over the interval [a, b] in (at least) n steps:
    //      auto Integrate(func, a, b, n)

    // Helper type traits
    namespace Internal
    {
        template<typename XTyA, typename XTyB>
        struct TQuadXType
        {
        private:
            using AddRes_t = TMatrixEvalType_t<decltype(std::declval<XTyB>() + std::declval<XTyA>())>;
        
        public:
            using type = std::conditional_t<std::is_integral_v<AddRes_t>, double, AddRes_t>;
        };
        template<typename XTyA, typename XTyB>
        using TQuadXType_t = typename TQuadXType<XTyA, XTyB>::type;

        template<typename Func, typename XTyA, typename XTyB>
        struct TQuadResultType
        {
            using type = TMatrixEvalType_t<decltype(
                std::declval<std::invoke_result_t<Func, TQuadXType_t<XTyA, XTyB>>>() * 
                std::declval<TDxLength_t<TQuadXType_t<XTyA, XTyB>>>())>;
        };
        template<typename Func, typename XTyA, typename XTyB>
        using TQuadResultType_t = typename TQuadResultType<Func, XTyA, XTyB>::type;
    }

    // Midpoint integrator
    class QuadMidpointPolicy
    {
    protected:
        ~QuadMidpointPolicy() = default;

    public:
        template<typename Func, typename XTyA, typename XTyB>
        static auto Integrate(Func&& func, XTyA&& a, XTyB&& b, std::size_t n)
        {
            // define types
            using XTy = Internal::TQuadXType_t<XTyA, XTyB>;
            using YTy = Internal::TQuadResultType_t<Func, XTyA, XTyB>;

            if (n == 0)
                return YTy{};
            
            XTy dx = static_cast<XTy>(b - a) / n;

            // use first iteration for initialization of the result (no addition needed)
            XTy x0 = a + dx/2;
            YTy result = std::invoke(func, x0);

            // execute integration
            for (std::size_t i = 1; i < n; i++)
                result += std::invoke(func, x0 + i*dx);
            
            result *= TDxLength<XTy>::Get(dx);
            return result;
        }
    };

    template<int Den, int W1, int W2, int... Ws>
    class TQuadNewtonCotesPolicy
    {
    protected:
        ~TQuadNewtonCotesPolicy() = default;

    public:
        template<typename Func, typename XTyA, typename XTyB>
        static auto Integrate(Func&& func, XTyA a, XTyB b, std::size_t n)
        {
            // define types
            using XTy = Internal::TQuadXType_t<XTyA, XTyB>;
            using YTy = Internal::TQuadResultType_t<Func, XTyA, XTyB>;
            using WeightList = TConstList<int, W1, W2, Ws...>;
            constexpr std::size_t wcnt = TConstListSizeof_v<WeightList>;

            // adjust n to fit the method used
            if (n < wcnt)
                n = wcnt;
            n -= (n - 1) % (wcnt - 1);

            XTy dx = static_cast<XTy>(b - a) / (n - 1);

            // handle borders
            constexpr auto wFirst = TConstListFront_v<WeightList>;
            constexpr auto wLast = TConstListBack_v<WeightList>;
            YTy result = wFirst*std::invoke(func, a) + wLast*std::invoke(func, b);

            using CLHelper = TConstListPopFront_t<TConstListPopBack_t<WeightList>>;
            using CL = TConstListAppend_t<int, wFirst + wLast, CLHelper>;
            constexpr std::size_t ccnt = TConstListSizeof_v<CL>;

            auto helper = [=, &result, &func](auto i)
            {
                constexpr std::size_t j0 = i + 1;
                YTy tmp = std::invoke(func, a + j0*dx);
                for (std::size_t j = j0 + ccnt; j < n - 1; j+=ccnt)
                    tmp += std::invoke(func, a + j*dx);
                result += TConstListGet_v<i, CL> * tmp;
            };
            ConstexprFor<std::size_t, 0, ccnt, 1>(helper);

            constexpr double den = static_cast<double>(Den) / ccnt;
            result *= TDxLength<XTy>::Get(dx) / den;
            return result;
        }
    };

    // see https://de.wikipedia.org/wiki/Newton-Cotes-Formeln#Abgeschlossene_Newton-Cotes-Formeln
    using QuadTrapezoidalPolicy = TQuadNewtonCotesPolicy<2, 1, 1>;
    using QuadSimpsonPolicy = TQuadNewtonCotesPolicy<6, 1, 4, 1>;
    using QuadSimpson38Policy = TQuadNewtonCotesPolicy<8, 1, 3, 3, 1>;
    using QuadBoolePolicy = TQuadNewtonCotesPolicy<90, 7, 32, 12, 32, 7>;
    using QuadNC6PointPolicy = TQuadNewtonCotesPolicy<288, 19, 75, 50, 50, 75, 19>;
    using QuadWeddlePolicy = TQuadNewtonCotesPolicy<840, 41, 216, 27, 272, 27, 216, 41>;
    using QuadNC8PointPolicy = TQuadNewtonCotesPolicy<17280, 751, 3577, 1323, 2989, 2989, 1323, 3577, 751>;

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename QuadPolicy=QuadSimpsonPolicy>
    class TQuadrature : public QuadPolicy {};

}

#endif
