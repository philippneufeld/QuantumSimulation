// Philipp Neufeld, 2021

#ifndef QSim_Math_Quadrature_H
#define QSim_Math_Quadrature_H

#include <cstdint>

namespace QSim
{

    // Rectangle integrator
    class QuadratureMidpoint
    {
    public:
        template<typename XTy, typename Func>
        auto Integrate(const Func& func, XTy x0, XTy x1, std::size_t cntFev)
        {
            using YTy = decltype(func(std::declval<XTy>()));
            if (cntFev == 0)
                return YTy{};
            
            XTy dx = (x1 - x0) / cntFev;
            XTy x = x0 + dx/2;

            // use first iteration for initialization of the result (no addition needed)
            YTy result = func(x);
            x += dx;

            // execute integration
            for (std::size_t i = 1; i < cntFev; i++, x+=dx)
                result += func(x);
            
            return result * dx;
        }
    };

}

#endif
