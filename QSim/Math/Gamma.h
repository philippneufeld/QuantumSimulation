// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Gamma_H
#define QSim_Math_Gamma_H

#include <cmath>
#include <limits>
#include <array>

namespace QSim
{
    // Gamma function variants calculator
    // See "Numerical recipes" book for reference
    class GammaFunction
    {
    private:
        template<typename Ty>
        constexpr static bool IsValidExp_v = 
            std::is_integral_v<Ty> || std::is_floating_point_v<Ty>;
        template<typename Ty1, typename Ty2=Ty1> 
        using EnableFunc_t = std::enable_if_t<IsValidExp_v<Ty1> && IsValidExp_v<Ty2>>;
        
        template<typename Ty1, typename Ty2=Ty1> 
        using EvalTy_t = std::conditional_t<
            std::is_integral_v<std::common_type_t<Ty1, Ty2>>, double, 
            std::decay_t<std::common_type_t<Ty1, Ty2>>>;

    public:
        // Normal gamma function
        // Fallback to the C++11 standard library implementation
        template<typename Ty, typename=EnableFunc_t<Ty>>
        static auto Gamma(Ty x) { return std::tgamma(x); }

        // Natural logarithm of the gamma function
        // Fallback to the C++11 standard library implementation
        template<typename Ty, typename=EnableFunc_t<Ty>>
        static auto GammaLn(Ty x) { return std::lgamma(x); }

        // Calculate the incomplete gamma function
        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto GammaP(Ty1 a, Ty2 x) 
        {
            using Ty = EvalTy_t<Ty1, Ty2>;
            if (x < 0.0 || a <= 0.0) return std::numeric_limits<Ty>::quiet_NaN();
            if (x == 0.0) return static_cast<Ty>(0.0);
            if (x > 100) return GammaPGL(a, x);
            if (x < a + 1.0) return GammaPSer(a, x);
            else return 1.0-GammaQSer(a, x);
        }

        // Calculate the composite incomplete gamma function
        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto GammaQ(Ty1 a, Ty2 x) 
        {
            return 1.0 - GammaP(a, x);
        }

        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto InvGammaP(Ty1 a, Ty2 p) 
        {
            using Ty = EvalTy_t<Ty1, Ty2>;
            constexpr Ty eps = 1e-8;
            if (a <= 0.0) return std::numeric_limits<Ty>::quiet_NaN();
            if (p >= 1) return std::max<Ty>(100, a+100*std::sqrt(a));
            if (p <= 0.0) return 0.0;

            Ty gln = GammaLn(a);

            // make an initial guess for x
            Ty x = 0;
            Ty a1 = a - 1.0;
            Ty lna1 = std::log(a1);
            Ty afac = std::exp(a1*(lna1 - 1) - gln);

            if (a <= 1)
            {
                Ty t = 1.0 - a*(0.532+0.12*a);
                if (p < t) x = std::pow(p/t, 1.0/a);
                else x = 1-std::log(1-(p-t)/(1-t));
            }
            if (a > 1)
            {
                Ty pp = (p < 0.5)? p : 1 - p;
                Ty t = std::sqrt(-2*std::log(pp));
                x = (2.30753+t*0.27061)/(1 + t*(0.99229 + t*0.04481)) - t;
                x = (p >= 0.5) ? x : -x;
                x = std::max<Ty>(1e-3, a*std::pow(1 - 1/(9*a) - x/(3*std::sqrt(a)), 3));
            }

            // 
            for (int i = 0; i < 12; i++)
            {
                if (x <= 0) return 0.0;

                Ty t = 0;
                if (a > 1) t = afac*std::exp(-(x-a1) + a1*(std::log(x) - lna1));
                else t = std::exp(-x + a1*std::log(x) - gln);
                Ty u = (GammaP(a, x) - p) / t;

                x -= (t = u/(1 - 0.5*std::min<Ty>(1, u*(a1/x - 1))));
                if (x <= 0) x = 0.5*(x + t);
                if (std::abs(t) < eps*x ) break;
            }
            
            return x;
        }


    private:
        // Calculate the incomplete gamma function by a series expansion
        // valid if x < a + 1
        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto GammaPSer(Ty1 a, Ty2 x) 
        {
            // P(a, b) = gamma(a, x)/Gamma(a)
            // gamma(a, x) = e^(-x) x^a sum_n (Gamma(a)/Gamma(a+1+n) * x^n)
            using Ty = EvalTy_t<Ty1, Ty2>;
            constexpr Ty eps = std::numeric_limits<Ty>::epsilon();
            Ty term = 1.0 / a;
            Ty sum = term;
            for (Ty ap=a+1; std::abs(term) > std::abs(sum)*eps; ap+=1)
            {
                term *= x/ap;
                sum += term;
            }
            return sum*std::exp(-x+std::log(x)*a-GammaLn(a));
        }

        // Calculate the composite incomplete gamma function by a series expansion
        // valid if x > a + 1
        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto GammaQSer(Ty1 a, Ty2 x) 
        {
            // 1/(x+1-a - 1(1-a)/(x+3-a - 2*(2-a)/(x+5-a - ...)))
            // calculation using modified Lentz's rule
            using Ty = EvalTy_t<Ty1, Ty2>;
            constexpr Ty eps = std::numeric_limits<Ty>::epsilon();
            constexpr Ty fpmin = std::numeric_limits<Ty>::min() / eps;
            Ty b = x + 1 - a;
            Ty c = 1.0 / fpmin;
            Ty d = 1 / b;
            Ty res = d;
            for (int n=1;;n++)
            {
                Ty an = -n*(n-a);
                b += 2;
                c = (std::abs(c) > fpmin ? b + an/c : b + an/fpmin);
                d = an*d + b;
                d = (std::abs(d) > fpmin ? 1.0 / d : 1.0 / fpmin);
                
                Ty cd = c*d;
                res *= cd;
                if (std::abs(cd - 1.0) <= eps)
                    break;
            }
            return res*std::exp(-x+std::log(x)*a-GammaLn(a));
        }

        // Calculate the composite incomplete gamma function by Gauss-Legendre quadrature
        // valid if x > a + 1
        template<typename Ty1, typename Ty2, typename=EnableFunc_t<Ty1, Ty2>>
        static auto GammaPGL(Ty1 a, Ty2 x) 
        {
            using Ty = EvalTy_t<Ty1, Ty2>;
            
            constexpr static int ngau = 18;
            constexpr static std::array<Ty, ngau> GLy_v = {
                0.0021695375159141994, 0.011413521097787704, 0.027972308950302116,
                0.051727015600492421, 0.082502225484340941, 0.12007019910960293,
                0.16415283300752470, 0.21442376986779355, 0.27051082840644336, 
                0.33199876341447887, 0.39843234186401943, 0.46931971407375483, 
                0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 
                0.78649910768313447, 0.87126389619061517, 0.95698180152629142
            };
            constexpr static std::array<Ty, ngau> GLw_v = {
                0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 
                0.027298621498568734, 0.034213810770299537, 0.040875750923643261,
                0.047235083490265582, 0.053244713977759692, 0.058860144245324798,
                0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
                0.076598410645870640, 0.079687828912071670, 0.082187266704339706,
                0.084078218979661945, 0.085346685739338721, 0.085983275670394821
            };

            Ty xu = 0;
            Ty a1 = a-1, lna1 = std::log(a1), sqrta1 = std::sqrt(a1);
            if (x > a1) xu = std::max<Ty>(a1 + 11.5*sqrta1, x + 6*sqrta1);
            else xu = std::max<Ty>(0, std::min<Ty>(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
            Ty sum = 0;
            for (int j=0;j<ngau;j++) 
            {
                Ty t = x + (xu-x)*GLy_v[j];
                sum += GLw_v[j]*std::exp(-(t-a1)+a1*(std::log(t)-lna1));
            }
            Ty ans = sum*(xu-x)*std::exp(a1*(lna1-1.0)-GammaLn(a));
            return ans>0.0 ? 1.0-ans : -ans;
        }

    };

}

#endif 
