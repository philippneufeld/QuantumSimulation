// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Constants.h>
#include <QSim/Rydberg/RydbergAtom.h>
#include <QSim/Rydberg/AtomStarkMap.h>

#include <QSim/Execution/ThreadPool.h>
#include <QSim/Util/ProgressBar.h>

using namespace QSim;
using namespace Eigen;

int main(int argc, const char* argv[])
{
    Rubidium atom;
    AtomStarkMap starkMap(atom, RydbergAtomState_t(28, 0, 0.5, 0.5), 23, 31, 20);

    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

    VectorXd eField = VectorXd::LinSpaced(600, 0.0, 60000.0);
    MatrixXd energies(eField.size(), cnt);
    MatrixXd overlaps(eField.size(), cnt);

    ThreadPool pool; 
    ProgressBar progress(eField.size());
    
    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            auto [ens, ovs] = starkMap.GetEnergies(eField[i]);
            energies.row(i) = ens / (PlanckConstant_v*SpeedOfLight_v) / 100;
            overlaps.row(i) = ovs;
            progress.IncrementCount();
        });
    }
    progress.WaitUntilFinished();

    // TODO: Save data

    return 0;
}
