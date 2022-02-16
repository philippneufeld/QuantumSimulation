// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Quadrature_H
#define QSim_Math_Quadrature_H

#include <cstdint>

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
            using type = TMatrixEvalType_t<decltype(std::declval<XTyB>() + std::declval<XTyA>())>;
        };
        template<typename XTyA, typename XTyB>
        using TQuadXType_t = typename TQuadXType<XTyA, XTyB>::type;

        template<typename Func, typename XTyA, typename XTyB>
        struct TQuadResultType
        {
            using type = TMatrixEvalType_t<decltype(
                std::declval<std::invoke_result_t<Func, TQuadXType_t<XTyA, XTyB>>>() * 
                std::declval<TQuadXType_t<XTyA, XTyB>>())>;
        };
        template<typename Func, typename XTyA, typename XTyB>
        using TQuadResultType_t = typename TQuadResultType<Func, XTyA, XTyB>::type;

        // Implements a unifying interface for creating scalar constants 
        // and matrix constants in which all elements have the same value
        template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
        struct TQuadValueLike
        {
            static auto Generate(TMatrixElementType_t<Ty> val, const Ty& like) 
            { 
                return TMatrixEvalType_t<Ty>::Ones(like.rows(), like.cols()) * val; 
            }
        };
        template<typename Ty>
        struct TQuadValueLike<Ty, false>
        {
            static auto Generate(TMatrixElementType_t<Ty> val, const Ty&) { return val; }
        };
        
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
            
            XTy dx = (b - a) / n;

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

            XTy dx = (b - a) / (n - 1);

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
            result *= dx / den;
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

    // Adaptive step size integrator
    class QuadAdaptivePolicy
    {
    protected:
        ~QuadAdaptivePolicy() = default;

    public:
        QuadAdaptivePolicy()
            : m_rtol(1e-8), m_atol(1e-12), m_depth(10) {}

        void SetIntegrationRTol(double rtol) { m_rtol = rtol; }
        void SetIntegrationATol(double atol) { m_atol = atol; }
        void SetIntegrationDepth(double depth) { m_depth = depth; }

        template<typename Func, typename XTyA, typename XTyB>
        auto Integrate(Func&& func, XTyA&& a, XTyB&& b, std::size_t n)
        {
            return Integrate(func, a, b, n, m_rtol, m_atol, m_depth);
        }

        template<typename Func, typename XTyA, typename XTyB>
        static auto IntegrateFevs(
            Func&& func, XTyA&& a, XTyB&& b, std::size_t n, 
            double rtol, double atol, std::size_t depth)
        {
            // define types
            using XTy = Internal::TQuadXType_t<XTyA, XTyB>;
            using YTy = Internal::TQuadResultType_t<Func, XTyA, XTyB>;

            // adjust n to fit the simpson method used
            if (n < 5)
                n = 5;
            n -= (n - 1) % 4;
            std::size_t sections = (n - 1) / 4;

            XTy dist = b - a;
            XTy dx = dist / sections;
            atol /= sections;

            // first iteration
            YTy f0 = std::invoke(func, a);
            YTy f1 = std::invoke(func, a + 0.5*dx);
            YTy f2 = std::invoke(func, a + dx);
            auto [res, sub_fevs] = IntegrateHelper(func, a, dx, 0, f0, f1, f2, rtol, atol, depth);
            YTy result = std::move(res);
            std::size_t fevs = 3 + sub_fevs;

            // rest of the iterations
            for (size_t i = 1; i < sections; i++)
            {
                f0 = f2;
                f1 = std::invoke(func, a + (i+0.5)*dx);
                f2 = std::invoke(func, a + (i+1)*dx);
                auto [res, sub_fevs] = IntegrateHelper(func, a+i*dx, dx, 0, f0, f1, f2, rtol, atol, depth);

                result += res;
                fevs += 2 + sub_fevs;
            }
            
            return std::make_pair(result, fevs);
        }

        template<typename Func, typename XTyA, typename XTyB>
        static auto Integrate(
            Func&& func, XTyA&& a, XTyB&& b, std::size_t n, 
            double rtol, double atol, std::size_t depth)
        {
            return IntegrateFevs(func, a, b, n, rtol, atol, depth).first;
        }

    private:
        template<typename Func, typename XTyA, typename DXTy>
        static std::pair<Internal::TQuadResultType_t<Func, XTyA, DXTy>, std::size_t> IntegrateHelper(
            Func& func, XTyA&& a, DXTy&& dx, std::size_t d,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f0,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f2,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f4,
            double rtol, double atol,
            std::size_t depth)
        {
            // define types
            using XTy = Internal::TQuadXType_t<XTyA, DXTy>;
            using YTy = Internal::TQuadResultType_t<Func, XTyA, DXTy>;
            
            std::size_t fevs = 2;

            // calculate area via Simpson's rule
            YTy I1 = (dx / 6) * (f0 + 4*f2 + f4);
            
            // calculate function value at intermediary points
            YTy f1 = std::invoke(func, a + 0.25*dx);
            YTy f3 = std::invoke(func, a + 0.75*dx);

            // calculate area via two Simpson's rule steps with half step size
            YTy I2 = (dx / 12) * (f0 + 4*(f1 + f3) + 2*f2 + f4);
            
            YTy Ierr = (I2 - I1) / 15;
            YTy I = I2 + Ierr; // add error to I2 (I is equivalent to Boole's rule)

            // check if accuracy is sufficient
            if (std::abs(Ierr) - Internal::TQuadValueLike<YTy>::Generate(atol, Ierr) > rtol*std::abs(I) && d < depth)
            {
                XTy dx2 = dx / 2;
                auto [res1, sub_fevs1] = IntegrateHelper(func, a, dx2, d+1, f0, f1, f2, rtol, atol / 2, depth);
                auto [res2, sub_fevs2] = IntegrateHelper(func, a+dx2, dx2, d+1, f2, f3, f4, rtol, atol / 2, depth);
                
                I = res1 + res2;
                fevs += sub_fevs1 + sub_fevs2;
            }

            return std::make_pair(I, fevs);
        }

    private:
        double m_rtol;
        double m_atol;
        std::size_t m_depth;
    };

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename QuadPolicy=QuadSimpsonPolicy>
    class TQuadrature : public QuadPolicy {};

}

#endif
