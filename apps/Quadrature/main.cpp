// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Eigen>

// #ifdef QSIM_PYTHON3
#include <functional>
#include <QSim/Constants.h>
#include <QSim/Math/Quad.h>
#include <QSim/Math/QuadAdapt.h>
#include <QSim/Python/Plotting.h>

#include <QSim/Util/DataFile.h>
#include <QSim/Util/PathUtil.h>

using namespace QSim;
using namespace Eigen;

template<typename Quad>
void GetIntegratorErrors(DataFileGroup group, std::function<double(double)> func, 
    double a, double b, double exact, int n0, int n1)
{
    VectorXi ns = VectorXi::Zero(n1 - n0 + 1);
    VectorXd res = VectorXd::Zero(ns.size());

    for (int n=n0, i=0; n<=n1; n++, i++)
    {
        ns[i] = n;
        res[i] = TQuadrature<Quad>::Integrate(func, a, b, n);
    };

    group.SetAttribute("a", a);
    group.SetAttribute("b", b);
    group.SetAttribute("exact", exact);
    group.CreateDataset("ns", ns);
    group.CreateDataset("res", res);
}

template<typename Quad>
void GetIntegratorSteps(DataFileGroup group, std::function<double(double)> func, 
    double a, double b, double exact, int n)
{
    std::vector<double> xs;
    std::vector<double> ys;

    auto func_dec = [&](double x)
    {
        xs.push_back(x);
        ys.push_back(func(x));
        return ys.back();
    };

    double res = TQuadrature<Quad>::Integrate(func_dec, a, b, n);
    
    group.SetAttribute("a", a);
    group.SetAttribute("b", b);
    group.SetAttribute("n", n);
    group.SetAttribute("exact", exact);
    group.SetAttribute("res", res);
    group.CreateDataset("xs", Map<VectorXd>(xs.data(), xs.size()).eval());
    group.CreateDataset("ys", Map<VectorXd>(ys.data(), ys.size()).eval());
}


int main()
{
    DataFile file;
    file.Open(GetHomeDirPath() + "/git/Masterarbeit/Data/Quad.h5", DataFile_TRUNCATE);
    auto root = file.OpenRootGroup();
    
    std::function<double(double)> func = [](double x){ return 1 / (Pi_v * (x*x + 1)); };
    auto [a, b, nmin, nmax] = std::make_tuple(-100.0, 100.0, 1, 500);
    double exact = TQuadrature<QuadMidpointPolicy>::Integrate(func, a, b, 10000000);

    auto errs = root.CreateSubgroup("Errors");
    GetIntegratorErrors<QuadSimpsonPolicy>(errs.CreateSubgroup("Simpson"), func, a, b, exact, nmin, nmax);
    GetIntegratorErrors<QuadFixedAdaptivePolicy>(errs.CreateSubgroup("Adaptive"), func, a, b, exact, nmin, nmax);

    int n = 55;
    auto steps = root.CreateSubgroup("Steps");
    GetIntegratorSteps<QuadSimpsonPolicy>(steps.CreateSubgroup("Simpson"), func, a, b, exact, n);
    GetIntegratorSteps<QuadFixedAdaptivePolicy>(steps.CreateSubgroup("Adaptive"), func, a, b, exact, n);

    VectorXd xs = VectorXd::LinSpaced(10000, a, b);
    VectorXd ys = xs.unaryExpr(func);
    root.CreateDataset("xs", xs);
    root.CreateDataset("ys", ys);

    return 0;
}

