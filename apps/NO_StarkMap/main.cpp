// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergDiatomic.h>
#include <QSim/Rydberg/DiatomicStarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

#include <QSim/Util/DataFile.h>
#include <QSim/Util/PathUtil.h>

#ifdef QSIM_PYTHON3
#include <QSim/Python/Plotting.h>
#endif

using namespace QSim;
using namespace Eigen;

void plotNOEnergyDiagram()
{
    NitricOxide molecule;

    int nmin = 22;
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
    ax.SetYLimits(-150.0, 50.0);
    ax.SetXLabel("Rotational quantum number $N^+$");
    ax.SetYLabel("$E/hc$ (cm${}^{-1}$)");
        
    matplotlib.RunGUILoop();
#endif
}

void plotNOStarkMap()
{

    constexpr double dE = 4 * EnergyInverseCm_v;

    NitricOxide molecule;
    DiatomicStarkMap starkMap(molecule, RydbergDiatomicState_t(38, 2, 2, 0, 0), 30, 60, 3, dE);
    
    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

    VectorXd eField = VectorXd::LinSpaced(600, 0.0, 25.0); // V cm^-1
    MatrixXd energies(eField.size(), cnt); // cm^-1

    ThreadPool pool; 
    ProgressBar progress(eField.size());
    
    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            energies.row(i) = starkMap.GetEnergies(eField[i] * 100) / EnergyInverseCm_v;
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

    
    // generate filename
    std::string path = GetHomeDirSubfolderPath("remote_home");
    path += "/Masterarbeit/06_StarkMap/03_NO/";
    path += GenerateFilename("NOStarkMap") + ".h5";
     
    DataFile file;
    file.Open(path, DataFile_MUST_NOT_EXIST);
    auto root = file.OpenRootGroup();
    root.CreateAttribute("Nmin", { 1 });
    auto energyDS = root.CreateDataset("Energies", { (std::size_t)energies.rows(), (std::size_t)energies.cols() });
    energyDS.StoreMatrix(energies);
    file.Close();


#ifdef QSIM_PYTHON3
    PythonMatplotlib matplotlib;
    auto fig = matplotlib.CreateFigure();
    auto ax = fig.AddSubplot();

    for (int i=0; i < cnt; i++)
        ax.Plot(eField.data(), energies.col(i).eval().data(), eField.size(), "", "C0.");

    ax.SetYLimits(-66.5, -61.5);
    ax.SetXLimits(0, eField.maxCoeff());
    ax.SetXLabel("Electric Field (V/cm)");
    ax.SetYLabel("Energy $\\frac{E}{hc}$ (cm$^{-1}$)");    

    matplotlib.RunGUILoop();
#endif
}

int main(int argc, const char* argv[])
{
    
    plotNOStarkMap();

    return 0;
}
