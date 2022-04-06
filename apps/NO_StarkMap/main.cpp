// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergDiatomic.h>
#include <QSim/Rydberg/DiatomicStarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    NitricOxide molecule;
    
    int nmin = 25;
    int nmax = 250;
    int Rmax = 4;

    int cnt = 0;
    for (int i = nmin; i <= nmax; i++) cnt += i;

    MatrixXd rydbergSeries(cnt, Rmax + 1);
    for (int R=0; R<=Rmax; R++)
    {
        int idx = 0;
        for (int n=nmin; n<=nmax; n++)
        {
            for (int l=0; l<n; l++)
            {
                RydbergDiatomicState_t state = std::make_tuple(n, l, R, 0, 0);
                rydbergSeries(idx++, R) = molecule.GetEnergy(state) / (100*PlanckConstant_v*SpeedOfLight_v); // cm^-1
            }
        }
    }

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;

    //
    // Level plot
    //
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();
    for (int i=0; i < rydbergSeries.cols(); i++)
    {
        VectorXd energies = rydbergSeries.col(i);
        VectorXd Ns = VectorXd::Ones(energies.size()) * i;
        ax.Plot(Ns.data(), energies.data(), Ns.size(), "", "C0.");
    }
    ax.SetYLimits(-60.0, 40.0);
    ax.SetXLabel("Rotational quantum number $N^+$");
    ax.SetYLabel("$E/hc$ (cm${}^{-1}$)");
        
    matplotlib.RunGUILoop();
#endif


    return 0;
}
