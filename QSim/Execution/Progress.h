// Philipp Neufeld, 2021-2022

#ifndef QSim_CLIProgress_H_
#define QSim_CLIProgress_H_

#include "../Platform.h"

#include <cstdint>
#include <mutex>
#include <condition_variable>

namespace QSim
{

    class Progress
    {
        using Timestamp_t = std::chrono::high_resolution_clock::time_point;
    public:

        // (cnt, total, elapsedTime, totalTimeEst)
        using Data_t = std::tuple<std::size_t, std::size_t, double, double>;

        Progress(std::size_t total);

        void SendSignal();
        void IncrementCount(std::size_t inc = 1);
        
        Data_t GetData() const;
        Data_t WaitForSignal() const;
        Data_t WaitForSignal(double timeout) const;
        Data_t WaitUntilFinished() const;

        std::pair<std::size_t, std::size_t> GetProgress() const;
        double GetTotalTimeEst() const;
        double GetElapsedTime() const;

    private:
        static Data_t CreateDataTuple(
            std::size_t cnt, std::size_t tot, 
            Timestamp_t startTs, Timestamp_t currTs);

    private:      
        // Threading
        mutable std::mutex m_mutex;
        mutable std::condition_variable m_signal;

        // Constants (no need for mutex)
        const Timestamp_t m_startTs;
        const std::size_t m_total; 

        // Progress bar properties
        std::size_t m_cnt;
        Timestamp_t m_currTs;
    };

}

#endif
