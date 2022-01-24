// Philipp Neufeld, 2021-2022

#ifndef QSim_Executor_ThreadPoolExecutor_H_
#define QSim_Executor_ThreadPoolExecutor_H_

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <atomic>

#include "../Platform.h"
#include "Executor.h"

namespace QSim
{

    using namespace std::chrono_literals;

    class ThreadPoolExecutor : public TExecutor<ThreadPoolExecutor>
    {
    public:
        ThreadPoolExecutor();
        ThreadPoolExecutor(std::size_t threadCnt);
        ~ThreadPoolExecutor();

        void AddTask(const std::function<void(void)>& task);
        void WaitUntilFinished();

    private:
        void ThreadFunc();

    private:
        std::vector<std::thread> m_threads; 
        std::queue<std::function<void(void)>> m_queue;
        
        mutable std::mutex m_mutex;
        std::condition_variable m_taskAdded;
        std::condition_variable m_taskFinished;
        std::atomic<unsigned int> m_ongoingTasks;
        bool m_stopThreads;
    };

    ThreadPoolExecutor::ThreadPoolExecutor()
        : ThreadPoolExecutor(std::thread::hardware_concurrency()) { }
    
    ThreadPoolExecutor::ThreadPoolExecutor(std::size_t threadCnt)
        : m_stopThreads(false), m_ongoingTasks(0)
    {
        threadCnt = threadCnt > 0 ? threadCnt : 1;
        for (unsigned int i = 0; i < threadCnt; i++)
            m_threads.push_back(std::thread([this](){ ThreadFunc(); }));
    }

    ThreadPoolExecutor::~ThreadPoolExecutor()
    {
        WaitUntilFinished();
        m_stopThreads = true;
        m_taskAdded.notify_all(); // wake up all threads
        for (auto& t: m_threads)
            t.join();
    }

    void ThreadPoolExecutor::AddTask(const std::function<void(void)>& task)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_queue.push(task);
        m_taskAdded.notify_one();
    }

    void ThreadPoolExecutor::WaitUntilFinished()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        while (m_ongoingTasks > 0 || !m_queue.empty())
            m_taskFinished.wait_for(lock, 10ms); // wake up every 10ms if nothing happened
    }

    void ThreadPoolExecutor::ThreadFunc()
    {
        while (true)
        {
            std::function<void(void)> func;
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                while (!m_stopThreads && m_queue.empty())
                    m_taskAdded.wait_for(lock, 10ms); // wake up every 10ms if nothing happened
                
                if (m_stopThreads)
                    break;

                m_ongoingTasks++;
                func = m_queue.front();
                m_queue.pop();
            }
            
            func(); // execute function

            m_ongoingTasks--;
            m_taskFinished.notify_all();
        }
    }

}

#endif