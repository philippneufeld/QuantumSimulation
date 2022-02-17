// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_LevenbergMarquard_H
#define QSim_Math_LevenbergMarquard_H

#include <type_traits>
#include <functional>
#include "MatrixTraits.h"
#include "Jacobian.h"

namespace QSim
{
    
    template<typename DiffPolicy=Diff1O2Policy>
    class LevenbergMarquard
    {
    private:

        template<typename Func, typename FTy, typename YTy, int N>
        static constexpr bool IsValidExpVHelper_v = 
            // all matrix types
            (TIsMatrix_v<FTy> && TIsMatrix_v<YTy>) &&
            // same number of elements (or dynamic)
            (TMatrixSizeAtCompileTime_v<FTy> == TMatrixSizeAtCompileTime_v<YTy> || 
             TMatrixSizeAtCompileTime_v<FTy> < 0 || TMatrixSizeAtCompileTime_v<YTy> < 0) &&
            // FTy is a vector type
            (TMatrixRowsAtCompileTime_v<FTy> == 1 || TMatrixColsAtCompileTime_v<FTy> == 1) &&
            // YTy is a vector type
            (TMatrixRowsAtCompileTime_v<YTy> == 1 || TMatrixColsAtCompileTime_v<YTy> == 1) &&
            // Function invocable
            std::is_invocable_r_v<FTy, Func, Eigen::Matrix<double, N, 1>>;

        template<typename Func, typename XTy, typename YTy, int N>
        static constexpr bool IsValidExpV_v = 
            // Helper
            IsValidExpVHelper_v<YTy(Eigen::Matrix<double, N, 1>), YTy, YTy, N> &&
            // all matrix types
            TIsMatrix_v<XTy> &&
            // same number of elements (or dynamic)
            (TMatrixSizeAtCompileTime_v<XTy> == TMatrixSizeAtCompileTime_v<YTy> || 
             TMatrixSizeAtCompileTime_v<XTy> < 0 || TMatrixSizeAtCompileTime_v<YTy> < 0) &&
            // XTy is a vector type
            (TMatrixRowsAtCompileTime_v<XTy> == 1 || TMatrixColsAtCompileTime_v<XTy> == 1) &&
            // Function check
            std::is_invocable_r_v<double, Func, TMatrixElementType_t<XTy>, Eigen::Matrix<double, N, 1>>;

        template<typename Func, typename YTy, int N>
        using EnableCFVH_t = std::enable_if_t<IsValidExpVHelper_v<Func, std::invoke_result_t<Func, Eigen::Matrix<double, N, 1>>, YTy, N>>;  
        template<typename Func, typename XTy, typename YTy, int N>
        using EnableCFV_t = std::enable_if_t<IsValidExpV_v<Func, XTy, YTy, N>>;   

    public:
        template<typename Func, typename XTy, typename YTy, int N, typename=EnableCFV_t<Func, XTy, YTy, N>>
        Eigen::Matrix<double, N, 1> CurveFitV(
            Func&& func, XTy&& x, YTy&& y,
            const Eigen::Matrix<double, N, 1>& beta,
            const Eigen::Matrix<double, N, 1>& step)
            // const Eigen::Matrix<double, N, 1>& order)
        {
            // compile time operations
            using FTy = Eigen::VectorXd;
            using VecTy = TMatrixEvalType_t<XTy>; // TODO: Improve

            // run-time checks
            assert(x.size() == y.size());

            const TMatrixEvalType_t<XTy>& xeval = x;
            auto funcV = [&](const Eigen::Matrix<double, N, 1>& beta)
            {
                FTy fs(xeval.size());
                for (int i = 0; i < xeval.size(); i++)
                    fs[i] = std::invoke(func, xeval[i], beta);
                return fs;
            };

            return CurveFitVHelper(funcV, std::forward<YTy>(y), beta, step);
        }

    private:
        template<typename Func, typename YTy, int N, typename=EnableCFVH_t<Func, YTy, N>>
        Eigen::Matrix<double, N, 1> CurveFitVHelper(
            Func&& func, YTy&& y,
            Eigen::Matrix<double, N, 1> beta,
            const Eigen::Matrix<double, N, 1>& step)
            // const Eigen::Matrix<double, N, 1>& order)
        {
            // define types
            using FuncRet = std::invoke_result_t<Func, Eigen::Matrix<double, N, 1>>;
            using FTy = TMatrixEvalType_t<FuncRet>; // Function values type
            using PTy = Eigen::Matrix<double, N, 1>; // Parameter type
            using JCalc = TJacobian<DiffPolicy>;
            using JTy = TMatrixEvalType_t<typename JCalc::template Jacobian_t<Func, PTy>>;
            using ATy = Eigen::Matrix<double, N, N>;

            double eps1 = 1e-6;
            double eps2 = 1e-6;
            double eps1sq = eps1*eps1;
            double eps2sq = eps2*eps2;

            const auto id = Eigen::Matrix<double, N, N>::Identity();
            const int kmax = 500;
            double lambda = 1e-3;
            
            // model deviation from the data
            FTy g = y - std::invoke(func, beta);

            // variables associated with the jacobian
            bool updateJacobian = true;
            JTy jac(y.size(), N); ATy a; ATy diag;

            for(int k = 0; k < kmax; k++)
            {
                // recalculate the jacobian since the value for beta has changed
                if (updateJacobian)
                {
                    JCalc::JacobianInPlace(jac, func, beta, step);
                    a = jac.transpose()*jac;
                    diag = a.diagonal().asDiagonal(); 
                    updateJacobian = false;
                }

                // calculate step
                PTy delta = (a+lambda*diag).colPivHouseholderQr().solve(jac.transpose()*g);
                
                
                
                FTy gnew = y - std::invoke(func, beta + delta);
                if (gnew.squaredNorm() > g.squaredNorm())
                {
                    // new model worse than old one -> increase lambda
                    // this increases the gradiant descent contribution
                    lambda *= 10;
                    continue;
                }
                else
                {
                    // new model better than old one -> decrease lambda and update beta
                    // this increases the Gauss-Newton contribution
                    lambda /= 10;
                    beta += delta;
                    g = std::move(gnew);
                    updateJacobian = true;
                }

                if (delta.norm() <= eps2*beta.norm()+eps2)
                    break;
            }

            return beta;
        }
    };
    
}

#endif
