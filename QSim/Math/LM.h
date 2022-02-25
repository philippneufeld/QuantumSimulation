// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_LM_H
#define QSim_Math_LM_H

#include <type_traits>
#include <functional>
#include <cmath>
#include "MatrixTraits.h"
#include "Jacobian.h"
#include "Gamma.h"

namespace QSim
{
    
    template<typename DiffPolicy=Diff1O2Policy>
    class LevenbergMarquardt
    {
    private:

        template<typename Func, typename FTy, typename YTy, int N>
        static constexpr bool IsValidExpVHelper_v = 
            // all matrix types
            (TIsMatrix_v<FTy> && TIsMatrix_v<YTy>) &&
            // same number of elements (or dynamic)
            (TMatrixCompileTimeSize_v<FTy> == TMatrixCompileTimeSize_v<YTy> || 
             TMatrixCompileTimeSize_v<FTy> < 0 || TMatrixCompileTimeSize_v<YTy> < 0) &&
            // FTy is a vector type
            (TMatrixCompileTimeRows_v<FTy> == 1 || TMatrixCompileTimeCols_v<FTy> == 1) &&
            // YTy is a vector type
            (TMatrixCompileTimeRows_v<YTy> == 1 || TMatrixCompileTimeCols_v<YTy> == 1) &&
            // Function invocable
            std::is_invocable_r_v<FTy, Func, Eigen::Matrix<double, N, 1>>;

        template<typename Func, typename XTy, typename YTy, int N>
        static constexpr bool IsValidExpV_v = 
            // Helper
            IsValidExpVHelper_v<YTy(Eigen::Matrix<double, N, 1>), YTy, YTy, N> &&
            // all matrix types
            TIsMatrix_v<XTy> &&
            // same number of elements (or dynamic)
            (TMatrixCompileTimeSize_v<XTy> == TMatrixCompileTimeSize_v<YTy> || 
             TMatrixCompileTimeSize_v<XTy> < 0 || TMatrixCompileTimeSize_v<YTy> < 0) &&
            // XTy is a vector type
            (TMatrixCompileTimeRows_v<XTy> == 1 || TMatrixCompileTimeCols_v<XTy> == 1) &&
            // Function check
            std::is_invocable_r_v<double, Func, TMatrixElementType_t<XTy>, Eigen::Matrix<double, N, 1>>;

        template<typename Func, typename YTy, int N>
        using EnableCFVH_t = std::enable_if_t<IsValidExpVHelper_v<Func, std::invoke_result_t<Func, Eigen::Matrix<double, N, 1>>, YTy, N>>;  
        template<typename Func, typename XTy, typename YTy, int N>
        using EnableCFV_t = std::enable_if_t<IsValidExpV_v<Func, XTy, YTy, N>>;   

    public:
        template<typename Func, typename XTy, typename YTy, int N, typename=EnableCFV_t<Func, XTy, YTy, N>>
        std::pair<Eigen::Matrix<double, N, 1>, Eigen::Matrix<double, N, N>> CurveFitV(
            Func&& func, XTy&& x, YTy&& y,
            const Eigen::Matrix<double, N, 1>& beta)
        {
            // no errors for the y values are known, use yerr=1 for the time being

            // compile time operations
            using FTy = Eigen::VectorXd;
            using VecTy = TMatrixEvalType_t<XTy>; // TODO: Improve

            // run-time checks
            assert(x.size() > N);
            assert(x.size() == y.size());

            const TMatrixEvalType_t<XTy>& xeval = x;
            auto funcV = [&](const Eigen::Matrix<double, N, 1>& beta)
            {
                FTy fs(xeval.size());
                for (int i = 0; i < xeval.size(); i++)
                    fs[i] = std::invoke(func, xeval[i], beta);
                return fs;
            };

            auto [params, cov, chi2] = CurveFitVHelper(funcV, std::forward<YTy>(y), beta/*, step*/);

            // estimate the errors
            Eigen::Matrix<double, N, N> adjCov = cov * chi2 / (2*(x.size() - N));

            // here, no independent "goodness-of-fit" metric is available
            return std::pair(params, adjCov);
        }

        template<typename Func, typename XTy, typename YTy, typename YErrTy, int N, typename=EnableCFV_t<Func, XTy, YTy, N>>
        std::tuple<Eigen::Matrix<double, N, 1>, Eigen::Matrix<double, N, N>, double> CurveFitV(
            Func&& func, XTy&& x, YTy&& y, YErrTy&& yerr,
            const Eigen::Matrix<double, N, 1>& beta)
        {
            // Weights are added to the function and the y values, such that 
            // chi^2 = sum_i ((y_i - f(x_i | beta)) / yerr_i)^2 = sum_i (y'_i - f'(x_i | beta))

            // compile time operations
            using FTy = Eigen::VectorXd;
            using XEvalTy = TMatrixEvalType_t<XTy>;
            using YEvalTy = TMatrixEvalType_t<YTy>;
            using YErrEvalTy = TMatrixEvalType_t<YErrTy>;

            // run-time checks
            assert(x.size() > N);
            assert(x.size() == y.size());

            const XEvalTy& xEval = x;
            const YErrEvalTy& weigths = YErrEvalTy::Ones(yerr.size()).cwiseQuotient(yerr);

            // calculate weighted function values
            auto funcV = [&](const Eigen::Matrix<double, N, 1>& beta) -> FTy
            {
                FTy fs(xEval.size());
                for (int i = 0; i < xEval.size(); i++)
                    fs[i] = std::invoke(func, xEval[i], beta);
                return fs.cwiseProduct(weigths).eval();
            };

            auto [params, cov, chi2] =  CurveFitVHelper(funcV, y.cwiseProduct(weigths), beta/*, step*/);

            // calculate the probability that the fit is good, given the data 
            // approximate via the chi-squared distribution (exact for linear models)
            auto nu = x.size() - N;
            double prob = 1 - GammaFunction::GammaP(0.5*nu, 0.5*chi2);

            return std::make_tuple(params, cov, prob);
        }

    private:
        template<typename Func, typename YTy, int N, typename=EnableCFVH_t<Func, YTy, N>>
        std::tuple<Eigen::Matrix<double, N, 1>, Eigen::Matrix<double, N, N>, double> CurveFitVHelper(
            Func&& func, YTy&& y, Eigen::Matrix<double, N, 1> beta)
        {
            // define types
            using FuncRet = std::invoke_result_t<Func, Eigen::Matrix<double, N, 1>>;
            using YEvalTy = TMatrixEvalType_t<YTy>;
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
            const int kmax = 100;
            const int ndone = 5;
            double rtol = 1e-3;
            double lambda = 1e-3;
            
            // model deviation from the data
            const YEvalTy& yEval = y;
            FTy chi = yEval - std::invoke(func, beta);
            double chi2 = chi.squaredNorm();

            // define jacobian and the curvature matrix (alpha)
            bool updateJacobian = true;
            JTy jac(y.size(), N); ATy alpha; ATy diag;

            // used to estimate the optimal stepsize (d(beta) = beta * sqrt(ulp))
            const double sqrtULP = std::sqrt(std::numeric_limits<double>::epsilon()); 

            for(int k = 0, done = 0; k < kmax && done <= ndone; k++)
            {
                // recalculate the jacobian and the curvature matrix, 
                // since the value for beta has changed
                if (updateJacobian)
                {
                    PTy step = sqrtULP * ((beta.array() != 0).select(beta, 1).matrix());

                    JCalc::JacobianInPlace(jac, func, beta, step);
                    alpha = jac.transpose()*jac;
                    diag = alpha.diagonal().asDiagonal(); 
                    updateJacobian = false;
                }

                // calculate step
                PTy delta = (alpha+lambda*diag).colPivHouseholderQr().solve(jac.transpose()*chi);
                
                FTy chiUp = yEval - std::invoke(func, beta + delta);
                double chi2Up = chiUp.squaredNorm();

                // check if step was small
                if (std::abs(chi2 - chi2Up) < rtol*std::max(chi2, chi2Up))
                    done++;

                if (chi2Up > chi2)
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
                    chi = std::move(chiUp);
                    chi2 = chi2Up;
                    updateJacobian = true;
                }

                // if (delta.norm() <= eps2*beta.norm()+eps2)
                //     break;
            }

            // estimate the covariance matrix
            ATy cov = alpha.inverse(); 

            return std::make_tuple(beta, cov, chi2);
        }
    };
    
}

#endif
