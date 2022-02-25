// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/HydrogenicSystem.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

#include <QSim/Math/Quad.h>

using namespace QSim;
using namespace Eigen;

class Numerov
{
public:
    template<typename KFunc, typename XTyA, typename XTyB, typename YTyA, typename YTyB>
    static void IntegrateInPlace(KFunc&& kfunc, XTyA xmin, XTyB xmax, std::size_t n, YTyA y0, YTyB y1)
    {
        using XTy = TMatrixAddResultFP_t<XTyA, XTyB>;
        //using YTy = TMatrixAddResultFP_t<



    }


};



double potential(double r, unsigned int n, unsigned int l)
{ 
    double E = -RydbergEnergy_v / (n*n);
    constexpr double k1 = 2*ElectronMass_v / ConstexprPow(ReducedPlanckConstant_v, 2);
    constexpr double k2 = ConstexprPow(ElementaryCharge_v, 2) / (4* Pi_v* VacuumPermittivity_v);

    return k1*(E + k2 / (r)) - l*(l+1) / (r*r);
}

double V(double r, int l)
{
    return l * (l + 1.0) / (r * r) - 2 / r;
}

std::pair<VectorXd, VectorXd> SimulateWavefunction(double dx, double xmax, unsigned int nE, unsigned int l)
{
	size_t n = (size_t)std::ceil(xmax / dx);

    VectorXd xs = VectorXd::LinSpaced(n, xmax, xmax/n + xmax*2/n);
    VectorXd psis = VectorXd::Zero(n);
    psis[0] = 0;
    psis[1] = 0.01;

    dx = xs[1]-xs[0];
    double c = dx*dx/12;

    // auto kfunc = [&](double x){return (E-V(x, l)); };
    auto kfunc = [&](double x){ return potential(x, nE, l); };

    for (int i = 2; i < xs.size(); i++)
    {
        double den = (1+c*kfunc(xs[i]));
        double s1 = 2*(1-5*c*kfunc(xs[i-1]))*psis[i-1];
        double s2 = (1+c*kfunc(xs[i-2]))*psis[i-2];
        psis[i] = (s1 - s2) / den;

        // psis[i] = ((2 - kfunc(xs[i]) * dx * dx) * psis[i-1] - (1 - dx / xs[i]) * psis[i-2]) / (1 + dx / xs[i]);
        // psis[i] = (2 - kfunc(xs[i-1]) * dx * dx) * psis[i-1] - psis[i-2];
    }
    
    psis = (4*Pi_v*xs.array()*xs.array()*psis.array()*psis.array()).matrix().eval();
    psis = psis / std::sqrt(psis.sum());

	return std::make_pair(xs, psis);
}



int main(int argc, const char* argv[])
{

    using type1 = typename TMatrixMulResultFP<short, VectorXd>::type;
    
    type1 t1;

    HydrogenicSystem hyd;

    // auto [r, y] = hyd.NumerovWF(potential, 0, 2, 1e-3, 1, 0);
    

    VectorXd rs = VectorXd::LinSpaced(1000, 0.1*BohrRadius_v, 2.5*BohrRadius_v);
    VectorXd ys = rs.unaryExpr([](double x){return potential(x, 1, 0);});

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    // ax.Plot(rs.data(), ys.data(), rs.size());
    
    for (int n=4; n<=15; n++)
    {
        for (int l=0; l < 1; l++)
        {
            auto [r, psi] = SimulateWavefunction(1e-2*BohrRadius_v, 3*n*n*BohrRadius_v, n, l);
            VectorXd x = r.cwiseSqrt();
            ax.Plot(x.data(), psi.data(), r.size());
        }
    }

    matplotlib.RunGUILoop();
#endif

    return 0;
}
