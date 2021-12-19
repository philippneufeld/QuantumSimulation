// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/Python/Plotting.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/NLevel/Doppler.h>
#include <QSim/Util/ThreadPool.h>

class CRb87SASApp : public QSim::CalcApp
{
    constexpr static double decay = 6.065e6;
public:
    CRb87SASApp()
    {
        // calculate parameters
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10e6);       

        // Create system
        std::array<std::string, 2> lvlNames = {"S1_2", "P3_2"};
        std::array<double, 2> levels = {0, QSim::SpeedOfLight_v / 780.241e-9};
        m_system = QSim::TNLevelSystemSC<2>(lvlNames, levels);
        m_system.SetDipoleElementByName("S1_2", "P3_2", dip);
        m_system.AddLaserByName("Probe", "S1_2", "P3_2", intProbe, false);
        m_system.AddLaserByName("Pump", "S1_2", "P3_2", intPump, true);
        m_system.SetDecayByName("P3_2", "S1_2", decay);
        
        m_doppler.SetMass(1.44316060e-25);
    }

    virtual void DoCalculation() override
    {
        QSim::ThreadPool pool;

        // dt << Rabi^-1, Doppler^-1, detuning^-1
        double dt = 1e-10;
        constexpr double tint = 15 / decay;

        auto laserDetunings = QSim::CreateLinspaceRow(-1e9, 1e9, 1001);
        QSim::TDynamicMatrix<double> detunings(2, laserDetunings.Size());
        detunings.SetRow(laserDetunings, m_system.GetLaserIdxByName("Probe"));
        detunings.SetRow(laserDetunings, m_system.GetLaserIdxByName("Pump"));

        auto start_ts = std::chrono::high_resolution_clock::now();

        auto rho0 = m_system.CreateGroundState();
        auto func = [&](auto dets)
        {
            return m_doppler.Integrate([&](double vel)
            { 
                auto rho = m_system.GetDensityMatrixAv(
                    dets, rho0, vel, 0.0, tint, 0.25*tint, dt);
                return rho.GetPopulation("P3_2");
            });
        };

        QSim::TDynamicRowVector<double> populations = pool.Map(
            func, detunings.GetColIterBegin(), detunings.GetColIterEnd());

        this->StoreMatrix("Detunings", detunings);
        this->StoreMatrix("Population P3/2", populations);
    }

    virtual void Plot() override
    {
        auto detunings = this->LoadMatrix("Detunings");
        auto x_axis = detunings.GetRow(m_system.GetLaserIdxByName("Probe"));
        auto y_axis = this->LoadMatrix("Population P3/2");
        
        QSim::PythonMatplotlib matplotlib;
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.Data(), y_axis.Data(), x_axis.Size());
        matplotlib.RunGUILoop();
    }

private:
    QSim::TNLevelSystemSC<2> m_system;
    QSim::DopplerIntegrator m_doppler;
};

int main(int argc, const char* argv[])
{
    CRb87SASApp app;
    return app.Run(argc, argv);
}

/*









// Philipp Neufeld, 2021

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>

#include <QSim/Python/Plotting.h>
#include <QSim/Util/Argparse.h>
#include <QSim/StaticQSys.h>
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
        constexpr double dip = 4.227 * QSim::ElementaryCharge_v * QSim::BohrRadius_v;
        constexpr double intProbe = QSim::GetIntensityFromRabiFrequency(dip, 3.5e6);
        constexpr double intPump = QSim::GetIntensityFromRabiFrequency(dip, 10e6);
        constexpr double decay = 6.065e6;
        constexpr double mass = 1.44316060e-25;

        // Create system
        QSim::TStaticQSys<2> system({"S1_2", "P3_2"}, {0, QSim::SpeedOfLight_v / 780.241e-9});
        system.SetDipoleElementByName("S1_2", "P3_2", dip);
        system.AddLaserByName("Probe", "S1_2", "P3_2", intProbe, false);
        system.AddLaserByName("Pump", "S1_2", "P3_2", intPump, true);
        system.SetDecayByName("P3_2", "S1_2", decay);
        system.SetMass(mass);

        // dt << Rabi^-1, Doppler^-1, detuning^-1
        double dt = 1e-10;
        constexpr double tint = 15 / decay;

        auto rho0 = system.CreateGroundState();

        auto laserDetunings = QSim::CreateLinspaceRow(-1e9, 1e9, 1001);
        QSim::TDynamicMatrix<double> detunings(2, laserDetunings.Size());
        QSim::SetRow(detunings, laserDetunings, 0);
        QSim::SetRow(detunings, laserDetunings, 1);

        auto start_ts = std::chrono::high_resolution_clock::now();

        auto func = [&](auto dets)
        { 
            auto rho = system.GetDensityMatrixAv(
                dets, rho0, 0.0, tint, 0.25*tint, dt);
            return rho.GetPopulation("P3_2");
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
        auto figure = matplotlib.CreateFigure();
        auto ax = figure.AddSubplot();
        ax.Plot(x_axis.data(), y_axis.data(), x_axis.size());
        matplotlib.RunGUILoop();
    }

    return 0;
}
*/