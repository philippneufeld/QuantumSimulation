// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Quad_H
#define QSim_Math_Quad_H

#include <cstdint>
#include <limits>
#include <Eigen/Dense>

#include "MatrixTraits.h"
#include "../Util/ConstList.h"
#include "../Util/ConstexprFor.h"

namespace QSim
{

    // Quadrature policy
    // Implements at least the following functions:
    // 1) Integrate func over the interval [a, b] in (at least) n steps:
    //      auto Integrate(func, a, b, n)

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
            using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
            using YTy = TMatrixMulResultFP_t<std::invoke_result_t<Func, XTy>, XTy>;

            if (n == 0)
                return YTy{};
            
            XTy dx = static_cast<XTy>(b - a) / n;

            // use first iteration for initialization of the result (no addition needed)
            XTy x0 = a + dx/2;
            YTy result = std::invoke(func, x0);

            // execute integration
            for (std::size_t i = 1; i < n; i++)
                result += std::invoke(func, x0 + i*dx);
            
            result *= TMatrixNorm<XTy>::Get(dx);
            return result;
        }
    };

    template<int Den1, typename W1s, int Den2, typename W2s>
    class TQuadNewtonCotesPolicy;
    template<int Den1, int W11, int W12, int... W1s, int Den2, int W21, int W22, int... W2s>
    class TQuadNewtonCotesPolicy<Den1, TConstList<int, W11, W12, W1s...>, Den2, TConstList<int, W21, W22, W2s...>>
    {
        using WeightList1_t = TConstList<int, W11, W12, W1s...>;
        using WeightList2_t = TConstList<int, W21, W22, W2s...>;
        static_assert(((TConstListSizeof_v<WeightList1_t> + 1) == TConstListSizeof_v<WeightList2_t>) ||
            (TConstListSizeof_v<WeightList1_t> == 2 && TConstListSizeof_v<WeightList2_t> == 2) && 
            TConstListSizeof_v<WeightList1_t> > 1);

        template<typename YTy, typename XTy>
        using RType_t = TMatrixEvalType_t<decltype(std::declval<YTy>()[0]*std::declval<XTy>())>;

        template<typename YTy>
        constexpr static bool IsValidData_v = TIsMatrix_v<YTy> && 
            (TMatrixCompileTimeRows_v<YTy> == 1 || TMatrixCompileTimeCols_v<YTy> == 1);
        template<typename Func, typename XTyA, typename XTyB>
        constexpr static bool IsValidFunc_v = !IsValidData_v<Func> && 
            std::is_invocable_v<Func, TMatrixAddResultFP_t<XTyA, XTyB>>;

    protected:
        ~TQuadNewtonCotesPolicy() = default;

    public:

        template<typename YTy, typename XTy, typename=std::enable_if_t<IsValidData_v<YTy>>>
        static auto Integrate(const YTy& y, XTy dx)
        {
            constexpr std::size_t wcnt = TConstListSizeof_v<WeightList1_t>;
            
            // check for zero or one size
            if (TMatrixCompileTimeSize_v<YTy> < 2 && y.size() < 2)
            {
                if (y.size() == 0) return RType_t<YTy, XTy>{};
                else return y[0] * TMatrixNorm<XTy>::Get(dx);
            }
            
            // use composite rule
            auto rem = (y.size() - 1) % (wcnt - 1);
            if (rem == 0)
            {
                return IntegrateHelper<Den1, WeightList1_t>(y) * TMatrixNorm<XTy>::Get(dx);
            }
            else if(y.size() - 1 > wcnt*rem)
            {
                return (IntegrateHelper<Den1, WeightList1_t>(y(Eigen::seqN(Eigen::fix<0>, y.size() - rem*wcnt, Eigen::fix<1>))) 
                    + IntegrateHelper<Den2, WeightList2_t>(y(Eigen::lastN(rem*wcnt + Eigen::fix<1>, Eigen::fix<1>)))) 
                    * TMatrixNorm<XTy>::Get(dx);
            }
            else if(y.size() - 1 == wcnt*rem)
            {
                return IntegrateHelper<Den2, WeightList2_t>(y) * TMatrixNorm<XTy>::Get(dx);
            }
            else
            {
                // fallback to trapezoidal rule
                return IntegrateHelper<2, TConstList<int, 1, 1>>(y) * TMatrixNorm<XTy>::Get(dx);
            }
        }

        template<typename Func, typename XTyA, typename XTyB, 
            typename=std::enable_if_t<IsValidFunc_v<Func, XTyA, XTyB>>>
        static auto Integrate(Func&& func, XTyA a, XTyB b, std::size_t n)
        {
            using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
            if (n < 0)
                return RType_t<Eigen::Matrix<std::invoke_result_t<Func, XTy>, 1, 1>, XTy>{};

            auto x = Eigen::Matrix<XTy, Eigen::Dynamic, 1>::LinSpaced(n, a, b);
            return Integrate(x.unaryExpr(func), TMatrixNorm<XTy>::Get((b-a) / (n-1)));
        }

    private:
        template<int Den, typename WeightList, typename YTy>
        static auto IntegrateHelper(YTy&& y)
        {
            // compile time assertion
            constexpr std::size_t wcnt = TConstListSizeof_v<WeightList>;
            constexpr int ysize = TMatrixCompileTimeSize_v<YTy>;
            static_assert(ysize < 0 || (ysize >= wcnt && ((ysize - 1) % (wcnt - 1)) == 0));

            // compile time actions
            constexpr auto wFirst = TConstListFront_v<WeightList>;
            constexpr auto wLast = TConstListBack_v<WeightList>;
            
            using CLHelper = TConstListPopFront_t<TConstListPopBack_t<WeightList>>;
            using CL = TConstListAppend_t<int, wFirst + wLast, CLHelper>;
            constexpr std::size_t ccnt = TConstListSizeof_v<CL>;

            constexpr double denRcp = static_cast<double>(ccnt) / Den;
            
            // runtime assertion
            assert((y.size() >= wcnt && ((y.size() - 1) % (wcnt - 1)) == 0));

            // carry out integration
            TMatrixElementType_t<YTy> result = wFirst*y[0] + wLast*y[y.size()-1];
            ConstexprFor<std::size_t, 0, ccnt, 1>([=, &result, &y](auto i)
            {
                result += TConstListGet_v<i, CL> * y(Eigen::seq(Eigen::fix<i+1>, y.size()-2, Eigen::fix<ccnt>)).sum();
            });
            return result * denRcp;
        }

    };

    // see https://de.wikipedia.org/wiki/Newton-Cotes-Formeln#Abgeschlossene_Newton-Cotes-Formeln
    using QuadTrapezoidalPolicy = TQuadNewtonCotesPolicy<2, TConstList<int, 1, 1>, 2, TConstList<int, 1, 1>>;
    using QuadSimpsonPolicy = TQuadNewtonCotesPolicy<6, TConstList<int, 1, 4, 1>, 8, TConstList<int, 1, 3, 3, 1>>;
    using QuadBoolePolicy = TQuadNewtonCotesPolicy<90, TConstList<int, 7, 32, 12, 32, 7>, 288, TConstList<int, 19, 75, 50, 50, 75, 19>>;
    using QuadWeddlePolicy = TQuadNewtonCotesPolicy<840, TConstList<int, 41, 216, 27, 272, 27, 216, 41>, 17280, TConstList<int, 751, 3577, 1323, 2989, 2989, 1323, 3577, 751>>;

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename QuadPolicy=QuadSimpsonPolicy>
    class TQuadrature : public QuadPolicy {};

}

#endif
