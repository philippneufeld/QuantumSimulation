// Philipp Neufeld, 2021

#include <iostream>
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
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10.0e6);
        constexpr double decay = 6.065e6;
        constexpr double mass = 1.44316060e-25;

        // Create system
        QSim::TNLevelSystemQM<3> system({"S1_2_F1", "S1_2_F2", "P3_2"}, {-4.271e9, 2.563e9, QSim::SpeedOfLight_v / 780.241e-9});
        system.SetDipoleElementByName("S1_2_F1", "P3_2", dip);
        system.SetDipoleElementByName("S1_2_F2", "P3_2", dip);
        system.AddLaserByName("Probe", "S1_2_F1", "P3_2", intProbe, false);
        system.AddLaserByName("Pump", "S1_2_F2", "P3_2", intPump, false);
        system.SetDecayByName("P3_2", "S1_2_F1", 3.0/8.0 * decay);
        system.SetDecayByName("P3_2", "S1_2_F2", 5.0/8.0 * decay);
        system.SetMass(mass);

        // Generate detuning axis
        constexpr static std::size_t cnt = 501;
        QSim::TStaticMatrix<double, 2, cnt> detunings;
        auto probeDetunings = QSim::CreateLinspaceRow(-100.0e6, 100.0e6, cnt);
        QSim::SetRow(detunings, probeDetunings, 0);

        auto start_ts = std::chrono::high_resolution_clock::now();

        auto func = [&](auto dets)
        { 
            auto rho = system.GetDensityMatrixSS(dets);
            return rho.GetAbsCoeff("S1_2_F1", "P3_2");
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
            for (std::size_t i = 0; i < probeDetunings.Size(); i++)
                file << probeDetunings[i] << " " << absCoeffs[i] << std::endl;
            file.close();
        }

        x_axis.assign(probeDetunings.Data(), probeDetunings.Data() + probeDetunings.Size());
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

