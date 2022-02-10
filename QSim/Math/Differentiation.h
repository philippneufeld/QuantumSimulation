// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Differentiation_H
#define QSim_Math_Differentiation_H

#include <cstdint>
#include <Eigen/Dense>

#include "../Util/CRTP.h"

namespace QSim
{
    // Differentiator CRTP base class
    // CRTP interface for different differentiator implementation that handles
    // single variable and multi-variable differentiation. Each implementation
    // is supposed to implement a "DifferentiateFast" method, which takes as an 
    // extra argument the length of the step (which is just the step in case of 
    // the single variable case and the vector norm in the multivariable case)
    template<typename IT, typename XTy>
    class TDifferentiator : public TCRTP<IT>
    {
    public:
        using XNormTy = XTy;

        template<typename Func>
        auto Differentiate(Func& func, XTy x, XTy dx)
        {
            return (~(*this)).DifferentiateFast(func, x, dx, dx);
        }
    };

    // Multivariable differentiator template specialization
    template<typename IT, typename Ty, int N, int M>
    class TDifferentiator<IT, Eigen::Matrix<Ty, N, M>> : public TCRTP<IT>
    {
    public:
        using XNormTy = Ty;

        template<typename Func, typename XTy, typename DXTy>
        auto Differentiate(Func& func, const Eigen::MatrixBase<XTy>& x, const Eigen::MatrixBase<DXTy>& dx)
        {
            return (~(*this)).DifferentiateFast(func, x, dx, dx.norm());
        }
    };


    //
    // Differentiator implementations
    //

    // Central difference differentiator
    // Differentiator that calculates the first derivative with a 
    // precision up to second order in dx
    template<typename XTy>
    class TDiff1O2 : public TDifferentiator<TDiff1O2<XTy>, XTy>
    {
        using Parent = TDifferentiator<TDiff1O2<XTy>, XTy>;
    public:
        template<typename Func>
        auto DifferentiateFast(Func& func, XTy x, XTy dx, typename Parent::XNormTy dxLen)
        {
            return (func(x + dx) - func(x - dx)) / (2*dxLen);
        }
    };

}

#endif
