// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Differentiation_H
#define QSim_Math_Differentiation_H

#include <cstdint>
#include <Eigen/Dense>

#include "../Util/CRTP.h"

namespace QSim
{
    // Midpoint integrator
    template<typename IT, typename XTy>
    class TDifferentiator : public TCRTP<IT>
    {
    public:
        using XNormTy = XTy;

        template<typename Func>
        auto Differentiate(Func& func, XTy x, XTy dx)
        {
            return (~(*this)).DifferentiateEx(func, x, dx, dx);
        }
    };

    // Midpoint integrator
    template<typename IT, typename Ty, int N, int M>
    class TDifferentiator<IT, Eigen::Matrix<Ty, N, M>> : public TCRTP<IT>
    {
    public:
        using XNormTy = Ty;

        template<typename Func, typename XTy, typename DXTy>
        auto Differentiate(Func& func, const Eigen::MatrixBase<XTy>& x, const Eigen::MatrixBase<DXTy>& dx)
        {
            return (~(*this)).DifferentiateEx(func, x, dx, dx.norm());
        }
    };

    // Midpoint integrator
    template<typename XTy>
    class TDiff1O2 : public TDifferentiator<TDiff1O2<XTy>, XTy>
    {
        using Parent = TDifferentiator<TDiff1O2<XTy>, XTy>;
    public:
        template<typename Func>
        auto DifferentiateEx(Func& func, XTy x, XTy dx, typename Parent::XNormTy dxLen)
        {
            return (func(x + dx) - func(x - dx)) / (2*dxLen);
        }
    };


}

#endif
