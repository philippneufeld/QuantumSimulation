// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_QuadAdapt_H
#define QSim_Math_QuadAdapt_H

#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>

#include "MatrixTraits.h"
#include "Quad.h"

namespace QSim
{
    // Quadrature policy (see Math/Quad.h)
    // Adaptive step size integrator
    class QuadAdaptivePolicy
    {
        template<typename Func, typename XTyA, typename XTyB>
        using YTy_t = TMatrixMulResultFP_t<
            std::invoke_result_t<Func, TMatrixAddResultFP_t<XTyA, XTyB>>, 
            TMatrixAddResultFP_t<XTyA, XTyB>>;

    protected:
        ~QuadAdaptivePolicy() = default;

    public:
        QuadAdaptivePolicy()
            : m_rtol(std::numeric_limits<double>::epsilon()), 
            m_atol(std::numeric_limits<double>::min()), m_depth(10) {}

        void SetIntegrationRTol(double rtol) { m_rtol = rtol; }
        void SetIntegrationATol(double atol) { m_atol = atol; }
        void SetIntegrationDepth(double depth) { m_depth = depth; }

        template<typename Func, typename XTyA, typename XTyB>
        auto Integrate(Func&& func, XTyA&& a, XTyB&& b, std::size_t n) const
        {
            return Integrate(func, a, b, n, m_rtol, m_atol, m_depth);
        }

        template<typename Func, typename XTyA, typename XTyB>
        static auto IntegrateFevs(
            Func&& func, XTyA&& a, XTyB&& b, std::size_t n, 
            double rtol, double atol, std::size_t depth)
        {
            // define types
            using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
            using YTy = YTy_t<Func, XTyA, XTyB>;

            // adjust n to fit the simpson method used
            if (n < 5)
                n = 5;
            n -= (n - 1) % 4;

            // separate range into sections
            std::size_t sections = (n - 1) / 4;
            XTy dx = static_cast<XTy>(b - a) / sections;
            auto dxlen = TMatrixNorm<XTy>::Get(dx);
            atol /= sections;

            // invoke functions
            std::vector<YTy> fs;
            fs.reserve(n);
            for (size_t i = 0; i < n; i++)
                fs.emplace_back(std::invoke(func, a+i*0.25*dx));
            YTy fsSum = std::accumulate(fs.begin() + 1, fs.end(), fs[0]);
            std::size_t fevs = n;

            // estimate using the trapezoid rule (divided by section count)
            YTy IestSec = (fsSum - 0.5*(fs[0]-fs[n-1])) * (0.25* dxlen / sections);

            // first iteration
            auto [res, sub_fevs] = IntegrateHelper(func, a, dx, dxlen, 
                1, fs[0], fs[1], fs[2], fs[3], fs[4], IestSec, rtol, atol, depth);
            YTy result = std::move(res);
            fevs += sub_fevs;

            // rest of the iterations
            for (size_t i = 1; i < sections; i++)
            {
                auto [res, sub_fevs] = IntegrateHelper(func, a+i*dx, dx, dxlen, 1, 
                    fs[4*i], fs[4*i+1], fs[4*i+2], fs[4*i+3], fs[4*i+4], IestSec, rtol, atol, depth);
                result += res;
                fevs += sub_fevs;
            }
            
            return std::make_pair(result, fevs);
        }

        template<typename Func, typename XTyA, typename XTyB>
        static auto Integrate(
            Func&& func, XTyA&& a, XTyB&& b, std::size_t n, 
            double rtol, double atol, std::size_t depth)
        {
            return IntegrateFevs(func, a, b, n, rtol, atol, depth).first;
        }

    private:
        template<typename Func, typename XTyA, typename DXTy>
        static std::pair<YTy_t<Func, XTyA, DXTy>, std::size_t> IntegrateHelper(
            Func& func, XTyA&& a, DXTy&& dx, TMatrixNorm_t<DXTy> dxlen, std::size_t d,
            const YTy_t<Func, XTyA, DXTy>& f0, const YTy_t<Func, XTyA, DXTy>& f1,
            const YTy_t<Func, XTyA, DXTy>& f2, const YTy_t<Func, XTyA, DXTy>& f3,
            const YTy_t<Func, XTyA, DXTy>& f4, const YTy_t<Func, XTyA, DXTy>& Iest,
            double rtol, double atol, std::size_t depth)
        {
            // define types
            using XTy = TMatrixAddResultFP_t<XTyA, DXTy>;
            using YTy = YTy_t<Func, XTyA, DXTy>;
            
            std::size_t fevs = 0;

            // calculate area via one and two Simpson's rules
            YTy I1 = (dxlen / 6) * (f0 + 4*f2 + f4);
            YTy I2 = (dxlen / 12) * (f0 + 4*(f1 + f3) + 2*f2 + f4);
            
            // add error to I2 (I is equivalent to Boole's rule)
            YTy Ierr = (I2 - I1) / 15;
            YTy I = I2 + Ierr;

            // check if accuracy is sufficient (otherwise recurse)
            YTy IerrAbs = TMatrixCwiseAbs<YTy>::Get(Ierr);
            YTy atolY = TMatrixOnesLike<YTy>::Get(Ierr) * atol;
            YTy IestAdj = TMatrixCwiseAbsMax<YTy>::Get(Iest, I);
            if (d < depth && TMatrixAnyCwiseLess<YTy>::Get(rtol*IestAdj, IerrAbs-atolY))
            {
                // calculate half of integration properties
                XTy dx2 = dx / 2;
                YTy Iest2 = IestAdj / 2;
                auto dxlen2 = dxlen / 2;

                // calculate intermediate points
                YTy f01 = std::invoke(func, a + 0.25*dx2);
                YTy f12 = std::invoke(func, a + 0.75*dx2);
                YTy f23 = std::invoke(func, a + 1.25*dx2);
                YTy f34 = std::invoke(func, a + 1.75*dx2);

                // NOTE:
                // static_cast<XTy>(a+dx2) in IntegrateHelper is necessary in order to
                // prevent the expression template to be passed directly and thus
                // exceeding the template recursion limit
                auto [res1, sub_fevs1] = IntegrateHelper(func, a, dx2, 
                    dxlen2, d+1, f0, f01, f1, f12, f2, Iest2, rtol, atol / 2, depth);
                auto [res2, sub_fevs2] = IntegrateHelper(func, static_cast<XTy>(a+dx2), 
                    dx2, dxlen2, d+1, f2, f23, f3, f34, f4, Iest2, rtol, atol / 2, depth);
                
                I = res1 + res2;
                fevs += 4 + sub_fevs1 + sub_fevs2;
            }

            return std::make_pair(I, fevs);
        }

    private:
        double m_rtol;
        double m_atol;
        std::size_t m_depth;
    };


    // Quadrature policy (see Math/Quad.h)
    // Adaptive step size integrator with fixed number of function evaluations
    class QuadFixedAdaptivePolicy
    {
        template<typename Func, typename XTyA, typename XTyB>
        using YTy_t = TMatrixMulResultFP_t<
            std::invoke_result_t<Func, TMatrixAddResultFP_t<XTyA, XTyB>>, 
            TMatrixAddResultFP_t<XTyA, XTyB>>;

        template<typename XTy, typename YTy>
        struct Section
        {
            YTy m_I;
            YTy m_Ierr;
            XTy m_a;
            std::array<YTy, 5> m_fs;
            int m_depth;

            Section(const TMatrixNorm_t<XTy>& dxlen, const YTy& f0, const YTy& f1, 
                const YTy& f2, const YTy& f3, const YTy& f4, const XTy& a, int depth)
                : m_fs{f0, f1, f2, f3, f4}, m_a(a), m_depth(depth)
            {
                YTy I1 = (dxlen / 6) * (m_fs[0] + 4*m_fs[2] + m_fs[4]);
                YTy I2 = (dxlen / 12) * (m_fs[0] + 4*(m_fs[1] + m_fs[3]) + 2*m_fs[2] + m_fs[4]);
                YTy Ierr = (I2 - I1) / 15;
                m_Ierr = TMatrixCwiseAbs<YTy>::Get(Ierr);
                m_I = I2 + Ierr;
            }

            bool operator<(const Section& rhs) const { return TMatrixAnyCwiseLess<YTy>::Get(m_Ierr, rhs.m_Ierr); }
        };

    protected:
        ~QuadFixedAdaptivePolicy() = default;

    public:
        QuadFixedAdaptivePolicy() = default;

        template<typename Func, typename XTyA, typename XTyB>
        static auto Integrate(Func&& func, XTyA&& a, XTyB&& b, std::size_t n)
        {
            return IntegrateFevs(func, a, b, n).first;
        }

        template<typename Func, typename XTyA, typename XTyB>
        static auto IntegrateFevs(Func&& func, XTyA&& a, XTyB&& b, std::size_t n)
        {
            // define types
            using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
            using YTy = YTy_t<Func, XTyA, XTyB>;

            // reduce count of points to be used in the coarse run
            std::size_t n0 = static_cast<std::size_t>(n / 2.5);

            // adjust n0 to fit the simpson method used
            if (n0 < 5)
                n0 = 5;
            n0 -= (n0 - 1) % 4;

            // separate range into sections
            std::size_t secCnt = (n0 - 1) / 4;
            XTy dx = static_cast<XTy>(b - a) / secCnt;
            auto dxlen = TMatrixNorm<XTy>::Get(dx);
            
            std::vector<Section<XTy, YTy>> queue_container;
            queue_container.reserve(static_cast<std::size_t>(std::ceil(n / 4.0)));
            std::priority_queue<Section<XTy, YTy>> section_queue(
                std::less<Section<XTy, YTy>>(), std::move(queue_container));

            YTy f0 = std::invoke(func, a);
            for (size_t i = 0; i < secCnt; i++)
            {
                YTy f1 = std::invoke(func, a + (i+0.25)*dx);
                YTy f2 = std::invoke(func, a + (i+0.5)*dx);
                YTy f3 = std::invoke(func, a + (i+0.75)*dx);
                YTy f4 = std::invoke(func, a + (i+1)*dx);
                section_queue.emplace(dxlen, f0, f1, f2, f3, f4, a + i*dx, 0);
                f0 = f4;
            }

            std::size_t fevs = n0;
            while (fevs < n)
            {
                auto sec = section_queue.top();
                section_queue.pop();

                double den = std::pow(2.0, sec.m_depth + 1);
                XTy dx2 = dx / den;
                auto dxlen2 = dxlen / den;

                // calculate intermediate points
                YTy f01 = std::invoke(func, sec.m_a + 0.25*dx2);
                YTy f12 = std::invoke(func, sec.m_a + 0.75*dx2);
                YTy f23 = std::invoke(func, sec.m_a + 1.25*dx2);
                YTy f34 = std::invoke(func, sec.m_a + 1.75*dx2);
                fevs += 4;

                section_queue.emplace(dxlen2, sec.m_fs[0], f01, sec.m_fs[1], f12, sec.m_fs[2], sec.m_a, sec.m_depth + 1);
                section_queue.emplace(dxlen2, sec.m_fs[2], f23, sec.m_fs[3], f34, sec.m_fs[4], sec.m_a+dx2, sec.m_depth + 1);
            }

            YTy res = section_queue.top().m_I;
            section_queue.pop();
            for (; !section_queue.empty(); section_queue.pop())
                res += section_queue.top().m_I;

            return std::make_pair(res, fevs);
        }
    };

}

#endif
