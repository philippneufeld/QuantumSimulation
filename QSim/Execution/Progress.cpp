// Philipp Neufeld, 2021-2022

#include "Progress.h"

namespace QSim
{

    Progress::Progress(std::size_t total)
        : m_startTs(std::chrono::high_resolution_clock::now()), 
        m_total(total), m_cnt(0), m_currTs(m_startTs) { }
    
    void Progress::SendSignal()
    {
        m_signal.notify_all();
    }

    void Progress::IncrementCount(std::size_t inc)
    {
        auto ts = std::chrono::high_resolution_clock::now();
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cnt = inc < m_total - m_cnt ? m_cnt + inc : m_total;
            m_currTs = ts;
        }
        m_signal.notify_all();
    }

    typename Progress::Data_t Progress::WaitForSignal() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_signal.wait(lock);
        return CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs);
    }

    typename Progress::Data_t Progress::WaitForSignal(double timeout) const
    {
        std::unique_lock<std::mutex> lock(m_mutex);  
        m_signal.wait_for(lock, std::chrono::duration<double>(timeout));
        return CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs);
    }

    typename Progress::Data_t Progress::WaitUntilFinished() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);  
        m_signal.wait(lock, [&](){ return m_cnt == m_total; });
        return CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs);
    }

    typename Progress::Data_t Progress::CreateDataTuple(
        std::size_t cnt, std::size_t tot, Timestamp_t startTs, Timestamp_t currTs)
    {
        auto ts = std::chrono::high_resolution_clock::now();
        double currSpan = (currTs - startTs).count() / 1e9;
        double elapsedTime = (cnt == tot ? currSpan : (ts - startTs).count() / 1e9);
        double totTimeEst = (cnt > 0 ? currSpan * tot / cnt : 2 * currSpan * tot);
        return std::make_tuple(cnt, tot, elapsedTime, totTimeEst);
    }

    typename Progress::Data_t Progress::GetData() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs);
    }

    std::pair<std::size_t, std::size_t> Progress::GetProgress() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return std::make_pair(m_cnt, m_total);
    }

    double Progress::GetTotalTimeEst() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return std::get<2>(CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs));
    }

    double Progress::GetElapsedTime() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return std::get<3>(CreateDataTuple(m_cnt, m_total, m_startTs, m_currTs));
    }

}
