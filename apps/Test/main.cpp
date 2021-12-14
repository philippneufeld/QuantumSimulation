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
        auto test = QSim::CreateLinspaceRow(0.0, 1.0, 101);
        this->StoreMatrix("Linspace", test);
    }

    virtual void Plot() override
    {
        auto data = LoadMatrix("Linspace");
        
        QSim::PythonMatplotlib matplotlib;
        auto fig = matplotlib.CreateFigure();
        auto ax = fig.AddSubplot();
        ax.Plot(data.Data(), data.Data(), data.Size());

        matplotlib.RunGUILoop();
    }
};

int main(int argc, const char* argv[])
{
    CTestApp app;
    return app.Run(argc, argv);
}
