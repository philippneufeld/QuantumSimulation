// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergAtom.h>
#include <QSim/Rydberg/StarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    Rubidium atom;
    AtomStarkMap starkMap(atom, 28, 0, 0.5, 0.5, 23, 31, 20);

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

    matplotlib.RunGUILoop();
#endif


    return 0;
}
