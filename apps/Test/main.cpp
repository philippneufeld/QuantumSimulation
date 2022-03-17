// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergSystem.h>
#include <QSim/Rydberg/StarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

class Hydrogen : public RydbergSystem
{
public:

    Hydrogen() : RydbergSystem(AtomicMassUnit_v) {}

    virtual double GetQuantumDefect(int n, int l, double j) const override
    {
        return 0.0;
    }

};

class Rubidium : public RydbergSystem
{
public:

    Rubidium() : RydbergSystem(84.9117897379*AtomicMassUnit_v) {}

    
    virtual double GetQuantumDefect(int n, int l, double j) const override
    {
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.67.052502 (l=0...2)
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.74.054502 (l=3)
        // see https://journals.aps.org/pra/pdf/10.1103/PhysRevA.102.062817 (l=4)

        using DefectCoeffs = std::array<double, 2>;
        constexpr DefectCoeffs coeffs1[] = { // coeffs for j = l - s
            {3.1311804, 0.1784}, // s
            {2.6548849, 0.2900}, // p
            {1.34809171, -0.60286}, // d
            {0.0165192, -0.085}, // f
            {0.0039990, -0.0202}, // g
        };
        constexpr DefectCoeffs coeffs2[] = { // coeffs for j = l + s
            {3.1311804, 0.1784}, // s
            {2.6416737, 0.2950}, // p
            {1.34646572, -0.59600}, // d
            {0.0165437, -0.086}, // f
            {0.0039990, -0.0202}, // g
        };

        // select right table
        double s = 0.5;
        int idx = static_cast<int>(std::round(((j - l) / s + 1) * 0.5));

        const DefectCoeffs* pDefectTable = nullptr;
        if (idx == 0) pDefectTable = coeffs1;
        else if(idx == 1) pDefectTable = coeffs2;
        else return 0.0;

        int ldef = std::min(l, 4);
        double defect = GetQuantumDefectFromCoeffs(n, pDefectTable[ldef].data(), 2);

        if (l > 4)
            defect = GetExtrapolatedQuantumDefect(defect, ldef, l);
        
        return defect;
    }

};

int main(int argc, const char* argv[])
{
    Rubidium atom;
    StarkMap starkMap(atom, 28, 0, 0.5, 0.5, 23, 31, 20);

    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

    VectorXd eField = VectorXd::LinSpaced(600, 0.0, 60000.0);
    MatrixXd energies(eField.size(), cnt);

    ThreadPool pool; 
    ProgressBar progress(eField.size());
    
    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            energies.row(i) = starkMap.GetEnergies(eField[i]) / (PlanckConstant_v*SpeedOfLight_v) / 100;
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    for (int i=0; i < cnt; i++)
        ax.Plot((eField / 100).eval().data(), energies.col(i).eval().data(), eField.size(), "", ".C0");

    ax.SetYLimits(-188.0, -167.5);
    // ax.SetYLimits(-150.0, -130.0);
    ax.SetXLimits(0, 600);
    ax.SetXLabel("Electric Field (V/cm)");
    ax.SetYLabel("Energy $\\frac{E}{hc}$ (cm$^{-1}$)");    

    /*double time = 0;
    int N = 600;
    for (int i=0; i < N; i++)
    {
        double el = i * 60000 / N;
        auto ts = std::chrono::high_resolution_clock::now();
        VectorXd energies = starkMap.GetEnergies(el) / (PlanckConstant_v*SpeedOfLight_v) / 100;
        time += (std::chrono::high_resolution_clock::now() - ts).count() / 1e6;
        VectorXd eField = VectorXd::Ones(cnt) * el / 100;
        
        ax.Plot(eField.data(), energies.data(), cnt, "", ".C0");
        // ax.SetYLimits(-150.0, -130.0);
        ax.SetYLimits(-188.0, -167.5);
        ax.SetXLimits(0, 600);
        ax.SetXLabel("Electric Field (V/cm)");
        ax.SetYLabel("Energy $\\frac{E}{hc}$ (cm$^{-1}$)");
    }
    
    std::cout << time / N << "ms per diagonalization step" << std::endl;*/

    matplotlib.RunGUILoop();
#endif


    return 0;
}
