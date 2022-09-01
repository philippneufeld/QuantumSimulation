// Philipp Neufeld, 2021-2022

#include <Eigen/Dense>

#include <QSim/Rydberg/RydbergDiatomic.h>
#include <QSim/Util/DataFile.h>
#include <QSim/Util/PathUtil.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace QSim;
using namespace Eigen;


int main(int argc, const char* argv[])
{
    NitricOxide rydberg;

    // 
    // Rydberg wavefunction
    //

    DataFile wfFile;
    wfFile.Open(GetHomeDirPath() + "/RydbergWF.h5", DataFile_TRUNCATE);
    auto root = wfFile.OpenRootGroup();

    std::pair<int, int> nl[] = {{25, 0}, {30, 0}, {35, 0}};
    for (auto [n, l]: nl)
    {
        auto state = RydbergDiatomicState_t(n, l, 0, l, 0);
        auto [r, R] = rydberg.GetRadialWF(state, 100);

        std::string name = "n=" + std::to_string(n) + " l=" + std::to_string(l);
        auto grp = root.CreateSubgroup(name);
        grp.CreateDataset("r", r);
        grp.CreateDataset("R", R);
    }

    wfFile.Close();

    //
    // Tikz graphics
    //

    std::vector<std::vector<double>> energies;
    energies.resize(5);

    for (int r=0; r<energies.size(); r++)
    {
        for (int n = 10; n<100; n++)
        {
            for (int l=0; l<std::min(n,5); l++)
            {
                double energy = rydberg.GetEnergy(RydbergDiatomicState_t(n, l, r, l, 0));
                energies[r].push_back(energy);
            }
        }
    }


    std::ofstream file;

    file.open(GetHomeDirPath() + "/git/Masterarbeit/Plots/RydbergLevels.tex", std::ios::out | std::ios::trunc);
    file << "\\begin{tikzpicture}" << std::endl;

    double height = 4.5;
    double lineWidth = 1;
    double lineSpacing = 0.1;

    double energy0 = -150.75 * EnergyInverseCm_v;
    double energyMax = 50 * EnergyInverseCm_v;
    double Nlabely = 0.275*height;

    file << "\\draw[-stealth, line width=1] (" << -lineWidth << "," << -1.05*height << ") -- (" 
         << -lineWidth << "," << 1.75*Nlabely << ");" << std::endl;
    file << "\\draw[line width=1] (" << -1.075*lineWidth << "," << -1.06*height << ") -- (" 
         << -0.925*lineWidth << "," << -1.04*height << ");" << std::endl;
    file << "\\draw[line width=1] (" << -1.075*lineWidth << "," << -1.08*height << ") -- (" 
         << -0.925*lineWidth << "," << -1.06*height << ");" << std::endl;


    file << "\\node[above, rotate=90] at (" << -1.75*lineWidth << "," << (-height+Nlabely)/2 << ") {" 
         << "Energy [$\\si{cm^{-1}}$]" << "};" << std::endl;

    for (double tick=30;; tick -= 30)
    {
        double tickEnergy = tick * EnergyInverseCm_v;
        if (tickEnergy < energy0) break;
        double y = height * tickEnergy / std::abs(energy0);
        file << "\\node[left] at (" << -lineWidth << "," << y << ") {" << tick << "};" << std::endl;
        file << "\\draw[line width=1] (" << -1.075*lineWidth << "," << y << ") -- (" 
         << -0.925*lineWidth << "," << y << ");" << std::endl;
    }

    file << "\\node[above] at (" << -lineWidth/2 << "," << Nlabely << ") {$N^+$};" << std::endl;

    for (int nu=0; nu<=0; nu++)
    {
        for (int R=0; R<=4; R++)
        {
            double x = (nu*5 + R)*(lineWidth+lineSpacing);
            file << "\\node[above] at (" << x + lineWidth/2 << "," << Nlabely << ")";
            file << " {" << R << "};" << std::endl;

            for (int n = 1; n<200; n++)
            {
                // for (int l=0; l<std::min(5, n); l++)
                for (int l=n-1; l<n; l++)
                {
                    double energy = rydberg.GetEnergy(RydbergDiatomicState_t(n, l, R, l, 0)) + nu*2376*EnergyInverseCm_v;
                    if (energy < energy0 || energy > energyMax) continue;

                    double y = height * energy / std::abs(energy0);

                    file << "\t\\draw[line width=1.0, color=cmap0] (";
                    file << x << "," << y << ") -- (";
                    file << x + lineWidth << "," << y << ");" << std::endl;
                }
            }
        }
    }
    
    file << "\\end{tikzpicture}" << std::endl;
    file.close();
    
}
