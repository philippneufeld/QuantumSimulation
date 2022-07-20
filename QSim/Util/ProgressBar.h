// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_ProgressBar_H_
#define QSim_Util_ProgressBar_H_

#include <cstdint>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>

namespace QSim
{

    std::size_t GetCLIWidth();

    class ProgressBar
    {
        using Timestamp_t = std::chrono::high_resolution_clock::time_point;
    public:
        ProgressBar(std::size_t total, const std::string& title="");
        ~ProgressBar();

        void WaitUntilFinished();
        void Update();

        void IncrementCount(std::size_t inc = 1);
        
    private:
        // Helper function fo the generation of the bar string
        static std::string TimeToString(double secs, bool millis);
        static std::string GetTimeText(double elapsedTime, double totalTimeEst, bool finished);
        static std::string GetProgressText(std::size_t cnt, std::size_t tot);
        static std::string GetProgressTextFrac(std::size_t cnt, std::size_t tot);
        std::string GenerateBar(std::size_t width, std::size_t cnt, std::size_t tot, 
            double elapsedTime, double totalTimeEst, bool printEst) const; 

        void WorkerThread();

    private:      
        // Threading
        std::thread m_thread;
        mutable std::mutex m_mutex;
        std::condition_variable m_wakeUp;
        bool m_stopThread;

        // Progress bar properties
        std::size_t m_cnt;
        Timestamp_t m_currTs;

        // Constants (no need for mutex)
        const Timestamp_t m_startTs;
        const std::size_t m_total; 
        const std::string m_title;
        const char m_progressChar;
    };

}

#endif
