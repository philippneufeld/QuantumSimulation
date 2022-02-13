// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_CLIProgressBar_H_
#define QSim_Util_CLIProgressBar_H_

#include "../Execution/Progress.h"

#include <cstdint>
#include <string>
#include <thread>
#include <atomic>
#include <memory>

namespace QSim
{

    std::size_t GetCLIWidth();

    class CLIProgBar
    {
        using Timestamp_t = std::chrono::high_resolution_clock::time_point;
    public:
        CLIProgBar(std::size_t total, const std::string& title="");
        CLIProgBar(std::shared_ptr<Progress> pHandle, const std::string& title="");
        ~CLIProgBar();

        operator std::shared_ptr<Progress>() { return m_pHandle; }
        operator std::shared_ptr<const Progress>() const { return m_pHandle; } 

        std::shared_ptr<Progress> GetHandle() { return m_pHandle; }
        std::shared_ptr<const Progress> GetHandle() const { return m_pHandle; } 

        void WaitOnDestruction(bool wait) { m_waitOnDtor = wait; }

        void WaitUntilFinished();
        void Update();

        void SetTitle();
        void IncrementCount(std::size_t inc = 1);
        
        void Start() {}
        void Stop() {}

    private:
        // Helper function fo the generation of the bar string
        static std::string TimeToString(double secs, bool millis);
        
        static std::string GetTimeText(double elapsedTime, double totalTimeEst, bool finished);
        
        static std::string GetProgressText(std::size_t cnt, std::size_t tot);
        static std::string GetProgressTextFrac(std::size_t cnt, std::size_t tot);
        
        std::string GenerateBar(std::size_t cnt, std::size_t tot, 
            double elapsedTime, double totalTimeEst) const;

        void WorkerThread();

    private:      
        // Threading
        std::thread m_thread;
        std::atomic<bool> m_stopThread;

        std::atomic<bool> m_waitOnDtor;
        std::shared_ptr<Progress> m_pHandle;

        // Constant attributes
        const std::string m_title;
        const std::size_t m_width; 
        const char m_progressChar;   
    };

}

#endif
