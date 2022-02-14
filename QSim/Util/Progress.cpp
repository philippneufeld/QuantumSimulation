// Philipp Neufeld, 2021-2022

#include "Progress.h"

namespace QSim
{

    Progress::Progress(std::size_t total)
        : m_total(total), m_cnt(0) { }
    
    void Progress::Reset()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cnt = 0;
        lock.unlock();
        m_signal.notify_all();
    }

    void Progress::Reset(std::size_t total)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cnt = 0;
        m_total = total;
        lock.unlock();
        m_signal.notify_all();
    }

    std::pair<std::size_t, std::size_t> Progress::WaitForUpdate() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        std::size_t cnt = m_cnt;
        std::size_t tot = m_total;
        m_signal.wait(lock, [&](){ return cnt < tot && cnt == m_cnt && tot == m_total; });
        return std::make_pair(m_cnt, m_total);
    }

    void Progress::WaitUntilFinished() const
    {
        while(true)
        {
            auto [cnt, tot] = WaitForUpdate();
            if (cnt < tot) break;
        }
    }

    void Progress::IncrementCount(std::size_t inc)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cnt = inc < m_total - m_cnt ? m_cnt + inc : m_total;
        lock.unlock();   
        m_signal.notify_all();
    }

    std::pair<std::size_t, std::size_t> Progress::GetProgress() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return std::make_pair(m_cnt, m_total);
    }

}
