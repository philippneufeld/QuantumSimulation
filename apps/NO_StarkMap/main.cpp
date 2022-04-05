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
    int Nmax = 4;

    int cnt = 0;
    for (int i = nmin; i <= nmax; i++) cnt += i;

    MatrixXd rydbergSeries(cnt, Nmax + 1);
    for (int N=0; N<=Nmax; N++)
    {
        int idx = 0;
        for (int n=nmin; n<=nmax; n++)
        {
            for (int l=0; l<n; l++)
            {
                RydbergDiatomicState_t state = std::make_tuple(n, l, 0, N, 0);
                rydbergSeries(idx++, N) = molecule.GetEnergy(state) / (PlanckConstant_v*SpeedOfLight_v);
            }
        }
    }

#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    for (int i=0; i < rydbergSeries.cols(); i++)
    {
        VectorXd energies = rydbergSeries.col(i);
        VectorXd Ns = VectorXd::Ones(energies.size()) * i;
        ax.Plot(Ns.data(), energies.data(), Ns.size(), "", "C0.");
    }
        
    matplotlib.RunGUILoop();
#endif


    return 0;
}
