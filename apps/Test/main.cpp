// Philipp Neufeld, 2021

#include <QSim/Util/CalcApp.h>
#include <QSim/NLevel/Laser.h>
#include <QSim/NLevel/NLevelSystem.h>
#include <QSim/Util/ThreadPool.h>
#include <QSim/Python/Plotting.h>

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


class CLIProgressBar
{
public:
    CLIProgressBar() : m_started(0), 
        m_progress(0), m_width(80) { }
    
    virtual void Start() 
    { 
        m_started = true;
        SetProgress(m_progress); 
    }
    
    void SetProgress(double prog) 
    { 
        m_progress = std::min(1.0, std::max(0.0, prog));

        if (m_started)
        {
            std::string str;
            str.reserve(m_width);
            
            str.append("[");
            std::size_t barWidth = m_width - 7;
            std::size_t iProg = static_cast<std::size_t>(barWidth * m_progress);
            str.append(iProg, '#');
            str.append(barWidth - iProg, ' ');
            str.append("] ");

            std::string strPerc = std::to_string(static_cast<std::size_t>(m_progress * 100));
            str.append(4 - strPerc.size(), ' ');
            str.append(strPerc);
            str.push_back('%');

            std::cout << "\u001b[1000D" << str << std::flush;

            if (m_progress == 1.0)
            {
                m_started = false;
                std::cout << std::endl;
            }
        }    
    }

private:
    std::string m_title = "";
    bool m_started;
    double m_progress;
    unsigned int m_width;
};


int main(int argc, const char* argv[])
{
    CLIProgressBar progress;
    progress.Start();
    for (std::size_t i = 0; i <= 100; i++)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        progress.SetProgress(i / 100.0);
    }

    // CTestApp app;
    // return app.Run(argc, argv);
}
