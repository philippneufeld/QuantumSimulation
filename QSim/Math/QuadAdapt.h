// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_QuadAdapt_H
#define QSim_Math_QuadAdapt_H

#include <cstdint>
#include <limits>

#include "MatrixTraits.h"
#include "Quad.h"

namespace QSim
{

    // Quadrature policy (see Math/Quad.h)
    
    // Helper type traits
    namespace Internal
    {
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


        // Implements a unifying interface for calculating the 
        // componentwise absolute value (works for matrices and scalars)
        template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
        struct TQuadCwiseAbs
        {
            static auto Abs(const Ty& mat) { return mat.cwiseAbs().eval(); }
        };
        template<typename Ty>
        struct TQuadCwiseAbs<Ty, false>
        {
            static auto Abs(const Ty& val) { return std::abs(val); }
        };

        // Implements a unifying interface for calculating the 
        // componentwise absolute maximum value (works for matrices and scalars)
        template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
        struct TQuadCwiseAbsMax
        {
            static auto AbsMax(const Ty& mat1, const Ty& mat2) 
            { 
                return (mat1.cwiseAbs()).cwiseMax((mat2.cwiseAbs())).eval(); 
            }
        };
        template<typename Ty>
        struct TQuadCwiseAbsMax<Ty, false>
        {
            static auto AbsMax(const Ty& val1, const Ty& val2) 
            { 
                return std::max(std::abs(val1), std::abs(val2)); 
            }
        };

        // Implements a unifying interface for calculating the 
        // componentwise absolute maximum value (works for matrices and scalars)
        template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
        struct TQuadAnyCwiseLess
        {
            static bool Compare(const Ty& mat1, const Ty& mat2) 
            { 
                return (mat1.array() < mat2.array()).any(); 
            }
        };
        template<typename Ty>
        struct TQuadAnyCwiseLess<Ty, false>
        {
            static bool Compare(const Ty& val1, const Ty& val2) 
            { 
                return val1 < val2; 
            }
        };  
    }

    // Adaptive step size integrator
    class QuadAdaptivePolicy
    {
    protected:
        ~QuadAdaptivePolicy() = default;

    public:
        QuadAdaptivePolicy()
            : m_rtol(std::numeric_limits<double>::epsilon()), 
            m_atol(std::numeric_limits<double>::min()), m_depth(10) {}

        void SetIntegrationRTol(double rtol) { m_rtol = rtol; }
        void SetIntegrationATol(double atol) { m_atol = atol; }
        void SetIntegrationDepth(double depth) { m_depth = depth; }

        template<typename Func, typename XTyA, typename XTyB>
        auto Integrate(Func&& func, XTyA&& a, XTyB&& b, std::size_t n) const
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

            XTy dx = static_cast<XTy>(b - a) / sections;
            atol /= sections;

            auto dxlen = TDxLength<XTy>::Get(dx);

            std::size_t fsCnt = 1+2*sections;
            Eigen::Matrix<YTy, -1, 1> fs(fsCnt);
            for (size_t i = 0; i < fsCnt; i++)
                fs[i] = std::invoke(func, a+i*0.5*dx);
            std::size_t fevs = fsCnt;

            // estimate using the trapezoid rule (divided by section count)
            YTy IestSec = (fs.sum() - 0.5*(fs[0]-fs[fsCnt-1])) * (0.5* dxlen / sections);

            // first iteration
            auto [res, sub_fevs] = IntegrateHelper(func, a, dx, dxlen, 
                1, fs[0], fs[1], fs[2], IestSec, rtol, atol, depth);
            YTy result = std::move(res);
            fevs += sub_fevs;

            // rest of the iterations
            for (size_t i = 1; i < sections; i++)
            {
                auto [res, sub_fevs] = IntegrateHelper(func, a+i*dx, dx, dxlen, 
                    1, fs[2*i], fs[2*i+1], fs[2*i+2], IestSec, rtol, atol, depth);
                result += res;
                fevs += sub_fevs;
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
            Func& func, XTyA&& a, DXTy&& dx, TDxLength_t<DXTy> dxlen, std::size_t d,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f0,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f2,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& f4,
            const Internal::TQuadResultType_t<Func, XTyA, DXTy>& Iest,
            double rtol, double atol, std::size_t depth)
        {
            // define types
            using XTy = Internal::TQuadXType_t<XTyA, DXTy>;
            using YTy = Internal::TQuadResultType_t<Func, XTyA, DXTy>;
            
            std::size_t fevs = 2;

            // calculate area via Simpson's rule
            YTy I1 = (dxlen / 6) * (f0 + 4*f2 + f4);
            
            // calculate function value at intermediary points
            YTy f1 = std::invoke(func, a + 0.25*dx);
            YTy f3 = std::invoke(func, a + 0.75*dx);

            // calculate area via two Simpson's rule steps with half step size
            YTy I2 = (dxlen / 12) * (f0 + 4*(f1 + f3) + 2*f2 + f4);
            
            YTy Ierr = (I2 - I1) / 15;
            YTy I = I2 + Ierr; // add error to I2 (I is equivalent to Boole's rule)

            // check if accuracy is sufficient (otherwise recurse)
            YTy IerrAbs = Internal::TQuadCwiseAbs<YTy>::Abs(Ierr);
            YTy atolY = Internal::TQuadValueLike<YTy>::Generate(atol, Ierr);
            YTy IestAdj = Internal::TQuadCwiseAbsMax<YTy>::AbsMax(Iest, I);
            if (d < depth && Internal::TQuadAnyCwiseLess<YTy>::Compare(rtol*IestAdj, IerrAbs-atolY))
            {
                // NOTE:
                // static_cast<XTy>(a+dx2) in IntegrateHelper is necessary in order to
                // prevent the expression template to be passed directly and thus
                // exceeding the template recursion limit
                XTy dx2 = dx / 2;
                YTy Iest2 = IestAdj / 2;
                auto dxlen2 = dxlen / 2;
                auto [res1, sub_fevs1] = IntegrateHelper(func, a, dx2, 
                    dxlen2, d+1, f0, f1,  f2, Iest2, rtol, atol / 2, depth);
                auto [res2, sub_fevs2] = IntegrateHelper(func, static_cast<XTy>(a+dx2), 
                    dx2, dxlen2, d+1, f2, f3, f4, Iest2, rtol, atol / 2, depth);
                
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

}

#endif
