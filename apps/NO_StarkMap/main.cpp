// Philipp Neufeld, 2021-2022

#include <iostream>
#include <array>
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

std::array<double, 5> stateToArray(const RydbergDiatomicState_t& state)
{
    auto createArray = [](auto&&... x) { return std::array<double, 5>{static_cast<double>(x)...}; };
    return std::apply(createArray, state);
}

int main(int argc, const char* argv[])
{
    NitricOxide molecule;
    RydbergDiatomicState_t state(38, 2, 2, 0, 0);
    constexpr double dE = 2 * EnergyInverseCm_v;
    DiatomicStarkMap starkMap(molecule, state, 30, 60, 3, dE);
    
    int cnt = starkMap.GetBasis().size();
    std::cout << "Basis size: " << cnt << std::endl;

    VectorXd eField = VectorXd::LinSpaced(600, 0.0, 25.0); // V cm^-1
    MatrixXd energiesData(eField.size(), cnt); // cm^-1

    // generate filename
    std::string path = GetHomeDirSubfolderPath("remote_home");
    path += "/Masterarbeit/06_StarkMap/03_NO/";
    path += GenerateFilename("NOStarkMap") + ".h5";
     
    DataFile file;
    file.Open(path, DataFile_MUST_NOT_EXIST);
    auto root = file.OpenRootGroup();
    
    // write stark map parameters
    root.CreateAttribute("State", { 5 });
    root.StoreAttribute("State", stateToArray(state).data());
    root.CreateAttribute("Energy_Range", { 1 });
    root.StoreAttribute("Energy_Range", &dE);

    std::mutex mutex;
    ThreadPool pool; 
    ProgressBar progress(eField.size());

    int groupNameLen = std::to_string(eField.size()).size();
    for (int i=0; i<eField.size(); i++)
    {
        pool.Submit([&,i=i](){
            auto [energies, states] = starkMap.GetEnergiesAndStates(eField[i] * 100);
            energies /= EnergyInverseCm_v;
            energiesData.row(i) = energies;

            // synchronize data file access
            std::unique_lock<std::mutex> lock(mutex);
            
            // create data group
            std::string groupName = std::to_string(i);
            groupName.insert(0, groupNameLen - groupName.size(), '0');
            auto group = root.CreateSubgroup(groupName);

            // define datasets
            group.CreateAttribute("Electric_Field", { 1 });
            auto energyStorage = group.CreateDataset("Energies", { static_cast<std::size_t>(energies.size()) });
            auto stateStorage = group.CreateDataset("States", { static_cast<std::size_t>(states.rows()), 
                static_cast<std::size_t>(states.cols()) });

            // store data
            group.StoreAttribute("Electric_Field", &eField[i]);
            energyStorage.Store(energies.data());
            stateStorage.StoreMatrix(states);
            
            progress.IncrementCount();
        });
    }

    progress.WaitUntilFinished();
    file.Close();

    return 0;
}
