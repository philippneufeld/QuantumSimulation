// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>
#include <QSim/Util/CLIProgressBar.h>

#include <chrono>

class CTestApp : public QSim::CalcApp
{
public:

    virtual void DoCalculation() override
    {
        
    }

    virtual void Plot() override
    {
        
    }
};

int main(int argc, const char* argv[])
{
    QSim::CLIProgBar progress;
    progress.Start();
        
    for (std::size_t i = 0; i <= 1000; i++)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        progress.SetProgress(i * 0.001);
    }
    
    progress.WaitUntilFinished();

    // CTestApp app;
    // return app.Run(argc, argv);
}
