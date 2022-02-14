// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_Progress_H_
#define QSim_Util_Progress_H_

#include "../Platform.h"

#include <cstdint>
#include <mutex>
#include <condition_variable>

namespace QSim
{

    // Threadsafe progress monitoring class
    class Progress
    {
    public:

        Progress(std::size_t total = 1);

        void Reset();
        void Reset(std::size_t total);

        void IncrementCount(std::size_t inc = 1);
        
        std::pair<std::size_t, std::size_t> WaitForUpdate() const;
        void WaitUntilFinished() const;

        std::pair<std::size_t, std::size_t> GetProgress() const;

    private:      
        // Threading
        mutable std::mutex m_mutex;
        mutable std::condition_variable m_signal;

        // Progress bar properties
        std::size_t m_total; 
        std::size_t m_cnt;
    };

}

#endif
