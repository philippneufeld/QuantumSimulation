// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_HydrogenicSystem_H_
#define QSim_Rydberg_HydrogenicSystem_H_

#include <cstdint>
#include <utility>
#include <cmath>

namespace QSim
{

    class HydrogenicSystem
    {
    public:

        // y''(x) = k(x)*y(x)
        template<typename Func>
        static auto NumerovWF(Func&& potential, double rmin, 
            double rmax, double step, double init1, double init2)
        {
            assert(rmax > rmin);
            step = std::abs(step);
            std::size_t n = static_cast<std::size_t>((rmax - rmin) / step) + 1;
            
            // variable transormation:
            // x = sqrt(r)
            // y(x) = x^{3/2} * R(x^2)

            auto kfunc = [&potential](double x)
            {
                return potential(x);
                double r = x*x;
                return 3.0 / (4.0 * r) + 4.0 * r * potential(r);
            };

            double xmin = rmin;
            double xmax = rmax;
            // double h = (xmax - xmin) / (n - 1);

            // double xmin = std::sqrt(rmin);
            // double xmax = std::sqrt(rmax);
            double h = (xmax - xmin) / (n - 1);
            
            Eigen::VectorXd rs(n);
            Eigen::VectorXd wfs(n);

            double yprev2 = init2;
            double yprev1 = init1;
            double kprev2 = kfunc(xmax + 2*h);
            double kprev1 = kfunc(xmax + h);

            double h12Div12 = h*h / 12.0;
            for (int i = 0; i < n; i++)
            {
                double x = xmax - i*h;
                double k = kfunc(x);
                double y = (2*(1-5*h12Div12*kprev1)*yprev1 - (1+h12Div12*kprev2)*yprev2) / (1+h12Div12*k);

                // double y = ((2 - potential(x) * h * h) * yprev1 - (1 - h / x) * yprev2) / (1 + h / x);
                
                double sqrtx = std::sqrt(x);
                rs[(n-1)-i] = x; //x*x;
                wfs[(n-1)-i] = x*x*y; //*sqrtx*sqrtx*sqrtx;
                
                // shift previous variables
                kprev2 = kprev1;
                yprev2 = yprev1;
                kprev1 = k;
                yprev1 = y;
            }
            
            return std::make_pair(rs, wfs);
        }
    };
    
}

#endif
