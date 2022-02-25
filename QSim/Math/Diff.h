// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Differentiation_H
#define QSim_Math_Differentiation_H

#include <functional>
#include "MatrixTraits.h"

namespace QSim
{
    // Differentiation policy
    // Implements at least the following functions:
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
                std::forward<DXTy>(dx), TMatrixNorm<XTy>::Get(dx));
        }

        template<typename Func, typename XTy, typename DXTy>
        static TMatrixEvalType_t<std::invoke_result_t<Func, XTy>>
            Differentiate(Func&& func, XTy&& x, DXTy&& dx, TMatrixNorm_t<XTy> dxLen)
        {
            return ((std::invoke(func, x + dx) - std::invoke(func, x - dx)) / (2*dxLen));
        }
    };

    // Helper class for usecases where inheritance from the policy is not desired
    template<typename DiffPolicy=Diff1O2Policy>
    class TDifferentiator : public DiffPolicy {};

}

#endif
