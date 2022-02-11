// Philipp Neufeld, 2021-2022

#ifndef QSim_CLIProgressBar_H_
#define QSim_CLIProgressBar_H_

#include "../Platform.h"

#include <string>
#include <iostream>
#include <atomic>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <chrono>

namespace QSim
{

    std::size_t GetCLIWidth();

    class CLIProgBar
    {
        using Timestamp_t = std::chrono::high_resolution_clock::time_point;
    public:
        CLIProgBar(std::size_t total, const std::string& title="");
        virtual ~CLIProgBar();

        void SetTitle();
        void SetCount(std::size_t cnt);
        void IncrementCount();
        
        void Start();
        void WaitUntilFinished();
        void Stop();
        void Update();

    private:
        // Helper function fo the generation of the bar string
        static std::string TimeToString(double secs, bool millis);
        std::string GetTimeText(std::size_t cnt, Timestamp_t cntTs) const;
        std::string GetProgressText(std::size_t cnt, bool percentage) const;
        std::string GenerateBar() const;

        void WorkerThread();

    private:      
        // Threading
        mutable std::mutex m_mutex;
        std::thread m_thread;
        std::condition_variable m_wakeUp;
        bool m_stopThread;

        // Progress bar properties
        bool m_started;
        std::size_t m_cnt;

        // Timing
        Timestamp_t m_startTs;
        Timestamp_t m_prevTs;

        // Constants (no need for mutex)
        const std::size_t m_total;
        const std::string m_title;
        const std::size_t m_width; 
        const char m_progressChar;   
    };

}

#endif
