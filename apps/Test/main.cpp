// Philipp Neufeld, 2021

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>

#include <QSim/Python/Plotting.h>
#include <QSim/Util/Argparse.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>

int main(int argc, const char* argv[])
{
    // parse command line arguents
    QSim::ArgumentParser parser;
    parser.AddOptionDefault("f,file", "Output filename.", "./data.txt");
    parser.AddOption("h,help", "Print this help string.");
    parser.AddOption("noplot", "Don't plot the results.");
    parser.AddOption("nocalc", "Just load the data from the specified file and skip the calculation");
    auto cmdArgs = parser.Parse(argc, argv);

    if (cmdArgs.IsError())
    {
        std::cout << cmdArgs.GetError() << std::endl;
        return -1;
    }
    else if (cmdArgs.IsOptionPresent("help"))
    {
        std::cout << parser.GetHelpString() << std::endl;
        return 0;
    }

    std::vector<double> x_axis;
    std::vector<double> y_axis;

    if (!cmdArgs.IsOptionPresent("nocalc"))
    {
        QSim::ThreadPool pool;

        // define parameters
        constexpr double lvl5S = 0;
        constexpr double lvl5P = lvl5S + QSim::SpeedOfLight_v / 780e-9;
        constexpr double lvl53D = lvl5P + QSim::SpeedOfLight_v / 480e-9;
        constexpr double lvl54P = lvl53D - 14.3e9;

        constexpr double dip5S5P = 5.178 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dip5P53D = 2.394e-2 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double dip54P53D = 3611 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;

        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip5S5P, 1.5e6);
        constexpr double intCoupling = QSim::GetIntensityFromRabiFrequency(dip5P53D, 0.25e6);
        constexpr double intMicrowave = QSim::GetIntensityFromRabiFrequency(dip54P53D, 10e6);

        // Create system
        QSim::TNLevelSystemQM<4> system({"5S", "5P", "54P", "53D"}, {lvl5S, lvl5P, lvl54P, lvl53D});
        system.SetMass(1.44316060e-25);
        system.SetDipoleElementByName("5S", "5P", dip5S5P);
        system.SetDipoleElementByName("5P", "53D", dip5P53D);
        system.SetDipoleElementByName("54P", "53D", dip54P53D);
        system.AddLaserByName("Probe", "5S", "5P", intProbe, false);
        system.AddLaserByName("Coupling", "5P", "53D", intCoupling, false);
        system.AddLaserByName("Microwaves", "54P", "53D", intMicrowave, false);
        system.SetDecayByName("5P", "5S", 6.065e6);
        
        // https://atomcalc.jqc.org.uk
        system.SetDecayByName("53D", "5P", 6.053e2);
        system.SetDecayByName("54P", "53D", 7.964e1);
        system.SetDecayByName("54P", "5S", 3.583e1);
        
        system.SetDopplerIntegrationSteps(35000);

        auto laserDetunings = QSim::CreateLinspaceRow(-3e7, 3e7, 201);
        QSim::TDynamicMatrix<double> detunings(system.GetLaserCount(), laserDetunings.Size());
        QSim::SetRow(detunings, laserDetunings, 0);
        
        auto start_ts = std::chrono::high_resolution_clock::now();

        auto func = [&](auto dets)
        { 
            auto rho = system.GetDensityMatrixSS(dets);
            return rho.GetAbsCoeff("5S", "5P");
        }; 
        auto absCoeffs = pool.Map(func, 
            QSim::GetColIteratorBegin(detunings), 
            QSim::GetColIteratorEnd(detunings));

        std::cout << "Calculation took " << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 << "s" << std::endl;

        // Write to file
        if (!cmdArgs.GetOptionStringValue("file").empty())
        {
            std::ofstream file;
            file.open(cmdArgs.GetOptionStringValue("file"), std::ios::out);
            for (std::size_t i = 0; i < laserDetunings.Size(); i++)
                file << laserDetunings[i] << " " << absCoeffs[i] << std::endl;
            file.close();
        }

        x_axis.assign(laserDetunings.Data(), laserDetunings.Data() + laserDetunings.Size());
        y_axis = absCoeffs;
    }
    else
    {
        // Load from file
        if (!cmdArgs.GetOptionStringValue("file").empty())
        {
            std::ifstream file;
            file.open(cmdArgs.GetOptionStringValue("file"), std::ios::in);
            while (!file.eof())
            {
                double det = 0;
                double absCoeff = 0;
                file >> det;
                file >> absCoeff;

                if (!file.fail())
                {
                    x_axis.push_back(det);
                    y_axis.push_back(absCoeff);
                }
                else
                {
                    break;
                }
            }
            file.close();
        }
    }

    // Plot
    if (!cmdArgs.IsOptionPresent("noplot"))
    {
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.MakeFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
    }

    return 0;
}
