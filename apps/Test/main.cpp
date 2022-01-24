// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/NLevel/NLevelSystem2.h>
#include <QSim/NLevel/NLevelSystemQM2.h>


#include <QSim/NLevel/NLevelSystemQM.h>

using namespace QSim;

int main(int argc, const char* argv[])
{
    constexpr double dip = 4.227 * ElementaryCharge_v * BohrRadius_v;
    constexpr double intProbe = GetIntensityFromRabiFrequency(dip, 3.5e6);
    

    TNLevelSystemQM<2> system1;
    system1.SetLevel(0, 0.0);
    system1.SetLevel(1, SpeedOfLight_v / 780.241e-9);
    system1.SetDecay(1, 0, 6.065e6);
    system1.SetDipoleElement(0, 1, dip);
    system1.AddLaser("test", 0, 1, intProbe, false);

    TNLevelSystemQM2<DynamicDim_v> system2(2);
    system2.SetLevel(0, 0.0);
    system2.SetLevel(1, SpeedOfLight_v / 780.241e-9);
    system2.SetDecay(1, 0, 6.065e6);
    system2.SetDipoleElement(0, 1, dip);
    system2.AddLaser(0, 1, intProbe, false);

    auto ss1 = system1.GetDensityMatrixSS(CreateZeros<double>(1), 0.0);
    for (std::size_t i = 0; i < ss1.Rows(); i++)
    {
        for (std::size_t j = 0; j < ss1.Cols(); j++)
            std::cout << ss1(i, j) << " ";
        std::cout << std::endl;
    }

    std::cout << system2.GetDensityMatrixSS(Eigen::VectorXd::Zero(1), 0.0) << std::endl;

    return 0;
}
