// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Numerov_H
#define QSim_Math_Numerov_H

#include <functional>
#include "MatrixTraits.h"

namespace QSim
{

    class Numerov
    {
    public:

        template<typename KFunc, typename XTyA, typename XTyB, typename YTyA, typename YTyB, typename DivFunc>
        static auto Integrate(KFunc&& kfunc, XTyA xstart, XTyB xstop, 
            std::size_t n, YTyA init1, YTyB init2, DivFunc divergenceHandler)
        {
            using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
            using YTy = TMatrixAddResultFP_t<
                TMatrixMulResultFP_t<YTyA, XTy>, TMatrixMulResultFP_t<YTyB, XTy>, 
                TMatrixMulResultFP_t<std::invoke_result_t<KFunc, XTy>, XTy>>;

            Eigen::Matrix<XTy, Eigen::Dynamic, 1> xs(n);
            Eigen::Matrix<YTy, Eigen::Dynamic, 1> ys(n);

            XTy dx = (xstop-xstart) / (n - 1);
            XTy c1 = dx*dx / 12.0;
            XTy c5 = 5*c1;

            YTy kprev2 = kfunc(xstart - 2*dx);
            YTy kprev1 = kfunc(xstart - dx);

            YTy yprev2 = init2;
            YTy yprev1 = init1;

            // generate x axis
            for (std::size_t i = 0; i < n; i++)
                xs[i] = xstart + i*dx;

            // start the evaluation with init1 and init2
            YTy kcurr; YTy ycurr;
            for (int i = 0; i < n; i++)
            {
                // execute numerov step
                kcurr = kfunc(xs[i]);
                ycurr = (2*(1-c5*kprev1)*yprev1 - (1+c1*kprev2)*yprev2) / (1+c1*kcurr);

                // update trailing variables
                yprev2 = yprev1;
                yprev1 = ycurr;
                kprev2 = kprev1;
                kprev1 = kcurr;

                // store data
                ys[i] = ycurr;

                if(!divergenceHandler(i, xs, ys))
                    break;
            }
            
            return std::make_pair(xs, ys);
        }

        template<typename KFunc, typename XTyA, typename XTyB, typename YTyA, typename YTyB>
        static auto Integrate(KFunc&& kfunc, XTyA xstart, XTyB xstop, std::size_t n, YTyA init1, YTyB init2)
        {
            return Integrate(kfunc, xstart, xstop, n, init1, init2, 
                [](auto i, auto& xs, auto& ys){ return true; });
        }
    };


}

#endif
