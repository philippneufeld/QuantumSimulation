// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Jacobian_H
#define QSim_Math_Jacobian_H

#include <functional>
#include "MatrixTraits.h"
#include "Differentiation.h"

namespace QSim
{
    template<typename DiffPolicy=Diff1O2Policy>
    class TJacobian : public DiffPolicy
    {
    private:
        template<typename Ty>
        using WrappedType_t = std::conditional_t<TIsMatrix_v<Ty>, 
            Ty, Eigen::Matrix<std::decay_t<Ty>, 1, 1>>;

        // Unwraps a type that has been wrapped inside of the Jacobian function
        template<bool wrapped>
        struct UnwrapHelper { template<typename T> static auto Do(T&& x) { return x[0]; }};
        template<>
        struct UnwrapHelper<false> { template<typename T> static auto Do(T&& x) { return x; }};

        template<bool strict, typename FTy, typename XTy, typename DXTy>
        static constexpr bool IsValidExpr_v = 
            // all matrix types
            (!strict || (TIsMatrix_v<FTy> && TIsMatrix_v<XTy> && TIsMatrix_v<DXTy>)) &&
            // same number of elements (or dynamic)
            (TMatrixSizeAtCompileTime_v<XTy> == TMatrixSizeAtCompileTime_v<DXTy> || 
             TMatrixSizeAtCompileTime_v<XTy> < 0 || TMatrixSizeAtCompileTime_v<DXTy> < 0) &&
            // XTy is a vector type
            (TMatrixRowsAtCompileTime_v<XTy> == 1 || TMatrixColsAtCompileTime_v<XTy> == 1) &&
            // DXTy is a vector type
            (TMatrixRowsAtCompileTime_v<DXTy> == 1 || TMatrixColsAtCompileTime_v<DXTy> == 1) &&
            // FTy is a vector type
            (TMatrixRowsAtCompileTime_v<FTy> == 1 || TMatrixColsAtCompileTime_v<FTy> == 1);

        template<bool strict, typename JTy, typename FTy, typename XTy, typename DXTy>
        static constexpr bool IsValidExprInPlace_v = IsValidExpr_v<strict, FTy, XTy, DXTy> &&
            // Is JTy capable of holding right amount of rows
            (TMatrixSizeAtCompileTime_v<FTy> == TMatrixRowsAtCompileTime_v<JTy> || 
             TMatrixSizeAtCompileTime_v<FTy> < 0 || TMatrixRowsAtCompileTime_v<JTy> < 0) &&
            // Is JTy capable of holding right amount of columns
            (TMatrixSizeAtCompileTime_v<XTy> == TMatrixColsAtCompileTime_v<JTy> || 
             TMatrixSizeAtCompileTime_v<XTy> < 0 || TMatrixColsAtCompileTime_v<JTy> < 0) &&
            (TMatrixSizeAtCompileTime_v<DXTy> == TMatrixColsAtCompileTime_v<JTy> || 
             TMatrixSizeAtCompileTime_v<DXTy> < 0 || TMatrixColsAtCompileTime_v<JTy> < 0);

        template<bool strict, typename Func, typename XTy, typename DXTy>
        using EnableJac_t = std::enable_if_t<IsValidExpr_v<strict, std::invoke_result_t<Func, XTy>, XTy, DXTy>>;
        template<bool strict, typename JTy, typename Func, typename XTy, typename DXTy>
        using EnableJacIP_t = std::enable_if_t<IsValidExprInPlace_v<strict, JTy, std::invoke_result_t<Func, XTy>, XTy, DXTy>>;

    public:
        // Calculates jacobian
        // Forwards parameters to JacobianHelper which only accepts 
        // Eigen matrix types, so if scalar types are involved, they are
        // wrapped like Eigen::Matrix<Scalar, 1, 1> and forwarded this way
        // to JacobianHelper
        template<typename Func, typename XTy, typename DXTy, typename=EnableJac_t<false, Func, XTy, DXTy>>
        static auto Jacobian(Func&& func, XTy&& x, DXTy&& dx)
        {
            return JacobianHelper(WrapFunction<XTy>(func), 
                static_cast<WrappedType_t<XTy>>(x), 
                static_cast<WrappedType_t<DXTy>>(dx));
        }

        // Calculates jacobian
        // Just like function Jacobian but instead of returning the
        // result it is stored into jac
        template<typename JTy, typename Func, typename XTy, typename DXTy, typename=EnableJacIP_t<false, JTy, Func, XTy, DXTy>>
        static auto JacobianInPlace(JTy& jac, Func&& func, XTy&& x, DXTy&& dx)
        {
            return JacobianInPlaceHelper(jac, WrapFunction<XTy>(func), 
                static_cast<WrappedType_t<XTy>>(x), 
                static_cast<WrappedType_t<DXTy>>(dx));
        }

    private:

        template<typename XTy, typename Func>
        static auto WrapFunction(Func& func)
        {
            using FTy = std::invoke_result_t<Func, XTy>;
            return [&](auto t) 
            { 
                auto unwrappedT = UnwrapHelper<!TIsMatrix_v<XTy>>::Do(t);
                return static_cast<WrappedType_t<FTy>>(std::invoke(func, unwrappedT)); 
            };
        }
         
        template<typename Func, typename XTy, typename DXTy, typename=EnableJac_t<true, Func, XTy, DXTy>>
        static auto JacobianHelper(Func&& func, XTy&& x, DXTy&& dx)
        {
            // J_ij = df_i / dx_j
            assert(x.rows() == 1 || x.cols() == 1); // x is vector
            assert(dx.rows() == 1 || dx.cols() == 1); // dx is vector
            assert(x.size() == dx.size()); // x and dx have same length
    
            auto df_dx0 = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), 0) * dx[0], dx[0]); 

            // setup jacobian matrix
            using ElementType = std::decay_t<decltype(df_dx0[0])>;
            constexpr int JRows = TMatrixSizeAtCompileTime_v<decltype(df_dx0)>;
            constexpr int JCols = TMatrixSizeAtCompileTime_v<XTy>;
            Eigen::Matrix<ElementType, JRows, JCols> jac(df_dx0.size(), x.size());

            for (int i = 0; i < jac.rows(); i++)
                jac(i, 0) = df_dx0[i];
            
            for (int j = 1; j < jac.cols(); j++)
            {
                auto df_dxj = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), j) * dx[j], dx[j]);
                for (int i = 0; i < jac.rows(); i++)
                    jac(i, j) = df_dxj[i];
            }

            return jac;
        }

        template<typename JTy, typename Func, typename XTy, typename DXTy, typename=EnableJacIP_t<true, JTy, Func, XTy, DXTy>>
        static auto JacobianInPlaceHelper(JTy& jac, Func&& func, XTy&& x, DXTy&& dx)
        {
            // J_ij = df_i / dx_j
            assert(x.rows() == 1 || x.cols() == 1); // x is vector
            assert(dx.rows() == 1 || dx.cols() == 1); // dx is vector
            assert(x.size() == dx.size()); // x and dx have same length
    
            auto df_dx0 = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), 0) * dx[0], dx[0]); 
            assert(jac.rows() == df_dx0.size());
            assert(jac.cols() == x.size());
            for (int i = 0; i < jac.rows(); i++)
                jac(i, 0) = df_dx0[i];
            
            for (int j = 1; j < jac.cols(); j++)
            {
                auto df_dxj = DiffPolicy::Differentiate(func, x, XTy::Unit(x.size(), j) * dx[j], dx[j]);
                for (int i = 0; i < jac.rows(); i++)
                    jac(i, j) = df_dxj[i];
            }
        }

    };
    
}

#endif
