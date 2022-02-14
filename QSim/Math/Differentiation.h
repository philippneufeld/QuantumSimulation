// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Differentiation_H
#define QSim_Math_Differentiation_H

#include <functional>
#include "MatrixTraits.h"

namespace QSim
{
    // Differentiation policy
    // Implements at least two functions:
    // 1) Differentiate func at position x with a step dx:
    //      auto Differentiate(func, x, dx)
    // 2) Same as 1) but with an auxilliary parameter specifying
    //    the length of dx:
    //      auto Differentiate(func, x, dx, dxLen)

    // Central difference differentiation policy
    // Differentiator that calculates the first derivative with a 
    // precision up to second order in dx
    // template<typename XTy>
    class Diff1O2Policy
    {
    protected:
        ~Diff1O2Policy() = default;
    
    public:
        template<typename Func, typename XTy, typename DXTy>
        static auto Differentiate(Func&& func, XTy&& x, DXTy&& dx)
        {
            return Differentiate(std::forward<Func>(func), std::forward<XTy>(x), 
                std::forward<DXTy>(dx), TDxLength<XTy>::Get(dx));
        }

        template<typename Func, typename XTy, typename DXTy>
        static TMatrixEvalType_t<std::invoke_result_t<Func, XTy>>
            Differentiate(Func&& func, XTy&& x, DXTy&& dx, TDxLength_t<XTy> dxLen)
        {
            return ((std::invoke(func, x + dx) - std::invoke(func, x - dx)) / (2*dxLen));
        }
    };

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename DiffPolicy=Diff1O2Policy>
    class TDifferentiator : public DiffPolicy {};


    namespace Internal
    {
        template<typename DiffPolicy=Diff1O2Policy>
        class TJacobianHelper : public DiffPolicy
        {
        private:
            // Unwraps a type that has been wrapped inside of the Jacobian function
            template<bool wrapped>
            struct UnwrapHelper { template<typename T> static auto Do(T&& x) { return x[0]; }};
            template<>
            struct UnwrapHelper<false> { template<typename T> static auto Do(T&& x) { return x; }};
                
        public:
            // Public API access point
            // Forwards parameters to JacobianHelper which only accepts 
            // Eigen matrix types, so if scalar types are involved, they are
            // wrapped like Eigen::Matrix<Scalar, 1, 1> and forwarded this way
            // to JacobianHelper
            template<typename Func, typename XTy, typename DXTy>
            static auto Jacobian(Func&& func, XTy&& x, DXTy&& dx)
            {
                using WrappedXTy = std::conditional_t<TIsMatrix_v<XTy>, 
                    XTy&&, Eigen::Matrix<std::decay_t<XTy>, 1, 1>>;
                using WrappedDXTy = std::conditional_t<TIsMatrix_v<DXTy>, 
                    DXTy&&, Eigen::Matrix<std::decay_t<DXTy>, 1, 1>>;

                using FTy = std::invoke_result_t<Func, XTy>;
                using VectorFTy = std::conditional_t<TIsMatrix_v<FTy>, 
                    FTy, Eigen::Matrix<std::decay_t<FTy>, 1, 1>>;

                auto adjFunc = [&](auto t) 
                { 
                    auto unwrappedT = UnwrapHelper<!TIsMatrix_v<XTy>>::Do(t);
                    return static_cast<VectorFTy>(std::invoke(func, unwrappedT)); 
                };

                return JacobianHelper(adjFunc, 
                    static_cast<WrappedXTy>(x), static_cast<WrappedDXTy>(dx));
            }

        private:

            template<typename FTy, typename XTy, typename DXTy>
            static constexpr bool IsValidExpr_v = 
                // all matrix types
                (TIsMatrix_v<FTy> && TIsMatrix_v<XTy> && TIsMatrix_v<DXTy>) &&
                // same number of elements (or dynamic)
                (TMatrixSizeAtCompileTime_v<XTy> == TMatrixSizeAtCompileTime_v<DXTy> 
                || TMatrixSizeAtCompileTime_v<XTy> < 0 || TMatrixSizeAtCompileTime_v<DXTy> < 0) &&
                // XTy is a vector type
                (TMatrixRowsAtCompileTime_v<XTy> == 1 || TMatrixColsAtCompileTime_v<XTy> == 1) &&
                // DXTy is a vector type
                (TMatrixRowsAtCompileTime_v<DXTy> == 1 || TMatrixColsAtCompileTime_v<DXTy> == 1) &&
                // FTy is a vector type
                (TMatrixRowsAtCompileTime_v<FTy> == 1 || TMatrixColsAtCompileTime_v<FTy> == 1);
                
            template<typename Func, typename XTy, typename DXTy>
            using EnableJacobian_t = std::enable_if_t<IsValidExpr_v<std::invoke_result_t<Func, XTy>, XTy, DXTy>>;

            template<typename Func, typename XTy, typename DXTy, typename=EnableJacobian_t<Func, XTy, DXTy>>
            static auto JacobianHelper(Func&& func, XTy&& x, DXTy&& dx)
            {
                // J_ij = df_i / dx_j
                assert(x.rows() == 1 || x.cols() == 1); // x is vector
                assert(dx.rows() == 1 || dx.cols() == 1); // dx is vector
                assert(x.size() == dx.size()); // x and dx have same length
        
                auto df_dx0 = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), 0) * dx[0], dx[0]); 
                // return df_dx1;

                using ElementType = std::decay_t<decltype(df_dx0[0])>;
                constexpr int JRows = TMatrixSizeAtCompileTime_v<decltype(df_dx0)>;
                constexpr int JCols = TMatrixSizeAtCompileTime_v<XTy>;

                Eigen::Matrix<ElementType, JRows, JCols> J(df_dx0.size(), x.size());

                for (int i = 0; i < J.rows(); i++)
                    J(i, 0) = df_dx0[i];
                
                for (int j = 1; j < J.cols(); j++)
                {
                    auto df_dxj = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), j) * dx[j], dx[j]);
                    for (int i = 0; i < J.rows(); i++)
                        J(i, j) = df_dxj[i];
                }

                return J;
            }

        };
    }
    
}

#endif
