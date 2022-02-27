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
        template<typename KFunc, typename XTyA, typename XTyB, typename YTyA, typename YTyB>
        static auto Integrate(KFunc&& kfunc, XTyA xstart, XTyB xstop, std::size_t n, YTyA init1, YTyB init2)
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

            // start the evaluation with init1 and init2
            

            YTy kcurr; YTy ycurr;
            for (std::size_t i = 0; i < n; i++)
            {
                // execute numerov step
                kcurr = kfunc(xstart + i*dx);
                ycurr = (2*(1-c5*kprev1)*yprev1 - (1+c1*kprev2)*yprev2) / (1+c1*kcurr);

                // update trailing variables
                yprev2 = yprev1;
                yprev1 = ycurr;
                kprev2 = kprev1;
                kprev1 = kcurr;

                // store data
                xs[i] = xstart + i*dx;
                ys[i] = ycurr;
            }
            
            return std::make_pair(xs, ys);
        }
    };


}

#endif
