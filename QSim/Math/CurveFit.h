// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_CurveFit_H_
#define QSim_Math_CurveFit_H_

#include <vector>
#include <utility>
#include <memory>

#include "../Util/Progress.h"
#include <Eigen/Dense>

namespace QSim
{

    namespace Internal
    {
		// invoke a function by expanding a array to a parameter pack
		template<int argCnt>
		struct TFunctionInvoker
		{
			template<typename sequence>
			struct TFunctionInvoker_impl;
			template<int... indices>
			struct TFunctionInvoker_impl<std::integer_sequence<int, indices...>>
			{
                template<int idx>
                static double GetElement(double* arrs) { return arrs[idx]; }

				template <typename Func>
				static double Invoke(Func func, double x, double* args)
				{
					return func(x, GetElement<indices>(args)...);
				}
			};

			template <typename Func>
			static double Invoke(Func func, double x, double* args)
			{
				return TFunctionInvoker_impl<decltype(std::make_index_sequence<argCnt>{}) > ::invoke(func, x, args);
			}
		};


        template<typename Func, int N>
        void CurveFitCalcJacobian(
            Func func, const Eigen::Ref<const Eigen::VectorXd>& x,
            Eigen::Matrix<double, Eigen::Dynamic, N>& J,
            const Eigen::Matrix<double, N, 1>& beta, 
            const Eigen::Matrix<double, N, 1>& step)
        {  
            int n = x.size();
            const double* beg = x.data();
            const double* end = beg + n;
            Eigen::VectorXd ly(n);
            Eigen::VectorXd ry(n);

            auto betaCpy = beta;
            for (int i = 0; i < N; i++)
            {
                betaCpy[i] = beta[i] - step[i];
                auto lfunc = [=](double t) { return func(t, betaCpy); };
                betaCpy[i] = beta[i] + step[i];
                auto rfunc = [=](double t) { return func(t, betaCpy); };

                for (int j = 0; j < n; j++)
                    ly[j] = func(beg[j], betaCpy);
                for (int j = 0; j < n; j++)
                    ry[j] = func(beg[j], betaCpy);

                for (int j = 0; j < n; j++)
                    J(j, i) = (ry[j] - ly[j]) / (2 * step[i]);                
            }
        }
    }

    template<typename Func, int N>
    Eigen::Matrix<double, N, 1> CurveFitV(
        Func func, const Eigen::Ref<const Eigen::VectorXd>& x, 
        const Eigen::Ref<const Eigen::VectorXd>& y, 
        const Eigen::Matrix<double, N, 1>& beta, 
        const Eigen::Matrix<double, N, 1>& step,
        const Eigen::Matrix<double, N, 1>& order)
    {
        int n = x.size();
        assert(x.size() == y.size());

        auto IdN = Eigen::Matrix<double, N, N>::Identity(N, N);
        int kmax = 250;

        double v = 2;
        double tau = 1e-6 * 1e-2;
        double eps1 = 1e-12 * 1e-3;
        double eps2 = 1e-12 * 1e-3;

        double eps1_sq = eps1*eps1;
        double eps2_sq = eps2*eps2;

        Eigen::Matrix<double, N, 1> currBeta = beta.array() / order.array();
        Eigen::Matrix<double, N, 1> newBeta = currBeta;

        // define function with "normalized parameters"
        auto nfunc = [&](double t, const Eigen::Matrix<double, N, 1>& params)
        { 
            return func(t, (params.array() * order.array()).matrix()); 
        };

        // calculate f vector
        Eigen::VectorXd f(n);
        for (int i = 0; i < n; i++)
            f[i] = nfunc(x[i], currBeta) - y[i];
        double currFerr = f.squaredNorm();

        Eigen::VectorXd newF = f;
        double newFerr = currFerr;

        // calculate jacobian
        Eigen::Matrix<double, Eigen::Dynamic, N> J(n, N);
        Internal::CurveFitCalcJacobian(nfunc, x, J, currBeta, step); 
        
        Eigen::Matrix<double, N, N> A = J.transpose() * J;
        Eigen::Matrix<double, N, 1> g = J.transpose() * f;    

        double l = tau * A.diagonal().maxCoeff();
        Eigen::Matrix<double, N, 1> delta(N);
        for (size_t k = 0; (g.squaredNorm() > eps1_sq) && k < kmax; k++)
        {
            // solve (A + l * I) * delta = -g
            delta = (A + l * IdN).colPivHouseholderQr().solve(-g);
            newBeta = currBeta + delta;

            // check termination condition
			if (delta.norm() <= eps2*(currBeta.norm() + eps2))
                break;

            // calculate new f vector
            for (int i = 0; i < n; i++)
                newF[i] = nfunc(x[i], newBeta) - y[i];
            newFerr = newF.squaredNorm();

            double dF = currFerr - newFerr;
            double dL = delta.dot(l*delta - g);
            double rho = dF / dL;
            
            if (rho > 0)
            {
                // update variables
                f = newF;
                newFerr = currFerr;
                currBeta = newBeta;

                Internal::CurveFitCalcJacobian(nfunc, x, J, currBeta, step); 
                A = J.transpose() * J;
                g = J.transpose() * f;

                v = 2;
                l *= std::max(1.0 / 3.0, 1 - std::pow(2*rho - 1, 3));
            }
            else
            {
                l *= v;
                v *= 2;
            }   
        }

        // reverse the normalization
        return currBeta.array() * order.array();
    }


    namespace Internal
    {
        void CurveFitAssignFromPack(double* dest) {}
        template<typename Arg, typename... Args>
        void CurveFitAssignFromPack(double* dest, Arg arg, Args... args)
        {
            *dest = arg;
            CurveFitAssignFromPack(dest + 1, std::forward<Args>(args)...);
        }

        void CurveFitAssignToPack(const double* dest) {}
        template<typename Arg, typename... Args>
        void CurveFitAssignToPack(const double* dest, Arg& arg, Args&... args)
        {
            arg = *dest;
            CurveFitAssignToPack(dest + 1, args...);
        }

        template<int n>
        struct CurveFitInvocationHelper
        {
            template<typename Func, typename... Args>
            static auto Invoke(Func func, const double* pArgs, Args... args)
            {
                return CurveFitInvocationHelper<n - 1>::Invoke(
                    func, pArgs + 1, std::forward<Args>(args)..., *pArgs);
            }
        };
        template<>
        struct CurveFitInvocationHelper<0>
        {
            template<typename Func, typename... Args>
            static auto Invoke(Func func, const double* pArgs, Args... args)
            {
                return func(std::forward<Args>(args)...);
            }
        };
    }

    template<typename Func, typename... Args>
    void CurveFit(
        Func func, 
        const Eigen::Ref<const Eigen::VectorXd>& x, 
        const Eigen::Ref<const Eigen::VectorXd>& y,
        Args&... args)
    {
        constexpr int N = sizeof...(Args);
        using VType = Eigen::Matrix<double, N, 1>;
        
        auto funcV = [&](double t, const VType& params)
        {
            return Internal::CurveFitInvocationHelper<sizeof...(Args)>::Invoke(
                func, params.data(), t);
        };

        VType beta(N);
        Internal::CurveFitAssignFromPack(&beta[0], std::forward<Args>(args)...);

        VType order(N);
        order = beta.cwiseAbs();

        VType steps(N);
        steps.setConstant(N, 1e-3);

        beta = CurveFitV(funcV, x, y, beta, steps, order);
        Internal::CurveFitAssignToPack(beta.data(), args...);
    }

}

#endif
