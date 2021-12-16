// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>

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
    CTestApp app;
    return app.Run(argc, argv);
}
