// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_Progress_H_
#define QSim_Execution_Progress_H_

#include <cstdint>
#include <mutex>
#include <condition_variable>

namespace QSim
{

    // Threadsafe progress tracking class
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
        // Synchronization
        mutable std::mutex m_mutex;
        mutable std::condition_variable m_signal;

        std::size_t m_total; 
        std::size_t m_cnt;
    };

    // Non thread safe progress tracker
    // Only use in single threaded contexts!
    class ProgressST
    {
    public:
        ProgressST(std::size_t total = 1) : m_total(total), m_cnt(0) {}

        void Reset() { m_cnt = 0; }
        void Reset(std::size_t total) { m_total = total; m_cnt = 0; }

        void IncrementCount(std::size_t inc = 1) { m_cnt = std::min(m_cnt + inc, m_total); }
        
        // These methods don't block (makes no sense for single threaded execution)
        std::pair<std::size_t, std::size_t> WaitForUpdate() const { return GetProgress(); }
        void WaitUntilFinished() const {}

        std::pair<std::size_t, std::size_t> GetProgress() const { return std::make_pair(m_cnt, m_total); }

    private:      
        std::size_t m_total; 
        std::size_t m_cnt;
    };

}

#endif
