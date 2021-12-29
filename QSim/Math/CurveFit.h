// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_CurveFit_H_
#define QSim_Math_CurveFit_H_

#include "Matrix.h"
#include "../Executor/Executor.h"

#include <vector>
#include <utility>

namespace QSim
{

    namespace Internal
    {
		// invoke a function by expanding a array to a parameter pack
		template<std::size_t argCnt>
		struct TFunctionInvoker
		{
			template<typename sequence>
			struct TFunctionInvoker_impl;
			template<std::size_t... indices>
			struct TFunctionInvoker_impl<std::integer_sequence<std::size_t, indices...>>
			{
                template<std::size_t idx>
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


        template<typename Func, std::size_t N, typename VT, typename Ex>
        void CurveFitCalcJacobian(
            Func func, const TVector<VT>& x,
            THybridMatrix<double, N, false>& J,
            const TStaticColVector<double, N>& beta, 
            const TStaticColVector<double, N>& step,
            TExecutor<Ex>& executor)
        {  
            std::size_t n = (~x).Size();
            const double* beg = (~x).Data();
            const double* end = beg + n;
            TDynamicColVector<double> ly(n);
            TDynamicColVector<double> ry(n);

            auto betaCpy = (~beta);
            for (std::size_t i = 0; i < N; i++)
            {
                betaCpy[i] = (~beta)[i] - (~step)[i];
                auto lfunc = [=](double t) { return func(t, betaCpy); };
                betaCpy[i] = (~beta)[i] + (~step)[i];
                auto rfunc = [=](double t) { return func(t, betaCpy); };

                (~executor).MapNonBlocking(lfunc, ly, beg, end);
                (~executor).MapNonBlocking(rfunc, ry, beg, end);
                (~executor).WaitUntilFinnished();

                for (std::size_t j = 0; j < n; j++)
                    J(j, i) = (ry[j] - ly[j]) / (2 * (~step)[i]);                
            }
        }
    }

    template<typename Func, std::size_t N, typename VT, typename Ex>
    TStaticColVector<double, N> CurveFitV(
        Func func, const TVector<VT>& x, const TVector<VT>& y, 
        const TStaticColVector<double, N>& beta, 
        const TStaticColVector<double, N>& step,
        TExecutor<Ex>& executor)
    {
        std::size_t n = (~x).Size();
        assert((~x).Size() == (~y).Size());

        auto IdN = CreateIdentityStatic<double, N>();
        std::size_t kmax = 250;

        double v = 2;
        double tau = 1e-6;
        double eps1 = 1e-12;
        double eps2 = 1e-12;

        double eps1_sq = eps1*eps1;
        double eps2_sq = eps2*eps2;

        auto currBeta = beta;
        auto newBeta = currBeta;

        // calculate f vector
        TDynamicColVector<double> f(n);
        for (std::size_t i = 0; i < n; i++)
            (~executor).AddTask([&, i](){ f[i] = func((~x)[i], ~currBeta) - (~y)[i]; });
        (~executor).WaitUntilFinnished();
        double currFerr = VectorLen2(f);

        TDynamicColVector<double> newF = f;
        double newFerr = currFerr;

        // calculate jacobian
        THybridMatrix<double, N, false> J(n, N);
        Internal::CurveFitCalcJacobian(func, x, J, ~currBeta, ~step, ~executor); 
        auto JT = J.Transpose();
        
        TStaticMatrix<double, N, N> A = JT * J;
        TStaticColVector<double, N> g = JT * f;    

        double l = tau * MatrixDiagMax(A);
        for (size_t k = 0; (VectorLen2(g) > eps1_sq) && k < kmax; k++)
        {
            // solve (A + l * I) * delta = -g
            auto delta = LinearSolve(A + l * IdN, -g);
            newBeta = currBeta + delta;

            // check termination condition
			if (VectorLen(delta) <= eps2*(VectorLen(beta) + eps2))
                break;

            // calculate new f vector
            for (std::size_t i = 0; i < n; i++)
                (~executor).AddTask([&, i](){ newF[i] = func((~x)[i], ~newBeta) - (~y)[i]; });
            (~executor).WaitUntilFinnished();
            newFerr = VectorLen2(newF);

            double dF = currFerr - newFerr;
            double dL = VectorDot(delta, l*delta - g);
            double rho = dF / dL;
            
            if (rho > 0)
            {
                // update variables
                f = newF;
                newFerr = currFerr;
                currBeta = newBeta;

                Internal::CurveFitCalcJacobian(func, x, J, ~currBeta, ~step, ~executor); 
                auto JT = J.Transpose();
                A = JT * J;
                g = JT * f;

                v = 2;
                l *= std::max(1.0 / 3.0, 1 - std::pow(2*rho - 1, 3));
            }
            else
            {
                l *= v;
                v *= 2;
            }   
        }

        return currBeta;
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

        template<std::size_t n>
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

    template<typename Func, typename VT, typename Ex, typename... Args>
    void CurveFit(
        TExecutor<Ex>& executor, Func func, 
        const TVector<VT>& x, const TVector<VT>& y, 
        Args&... args)
    {
        constexpr std::size_t N = sizeof...(Args);
        using VType = TStaticColVector<double, N>;
        
        auto funcV = [&](double t, const VType& params)
        {
            return Internal::CurveFitInvocationHelper<sizeof...(Args)>::Invoke(
                func, params.Data(), t);
        };

        VType beta(N);
        Internal::CurveFitAssignFromPack(&beta[0], std::forward<Args>(args)...);

        VType steps(N);
        MatrixAbs(steps, beta / 1e3);

        beta = CurveFitV(funcV, x, y, beta, steps, ~executor);
        Internal::CurveFitAssignToPack(beta.Data(), args...);
    }

}

#endif
