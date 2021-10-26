// Philipp Neufeld, 2021

#ifndef QSIM_ThreadPool_H_
#define QSIM_ThreadPool_H_

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <atomic>
#include <iostream>

namespace QSim
{


    class ThreadPool
    {
    public:
        ThreadPool();
        ~ThreadPool();

        void AddTask(const std::function<void(void)>& task);
        void WaitUntilFinnished();

    private:
        void ThreadFunc();

    private:
        std::vector<std::thread> m_threads; 
        std::queue<std::function<void(void)>> m_queue;
        
        mutable std::mutex m_mutex;
        std::condition_variable m_taskAdded;
        std::condition_variable m_taskFinnished;
        std::atomic<unsigned int> m_ongoingTasks;
        bool m_stopThreads;
    };

    ThreadPool::ThreadPool()
        : m_stopThreads(false), m_ongoingTasks(0)
    {
        for (unsigned int i = 0; i < std::thread::hardware_concurrency(); i++)
            m_threads.push_back(std::thread([this]{ ThreadFunc(); }));
    }

    ThreadPool::~ThreadPool()
    {
        WaitUntilFinnished();
        m_stopThreads = true;
        m_taskAdded.notify_all(); // wake up all threads
        for (auto& t: m_threads)
            t.join();
    }

    void ThreadPool::AddTask(const std::function<void(void)>& task)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_queue.push(task);
        m_taskAdded.notify_one();
    }

    void ThreadPool::WaitUntilFinnished()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        while (m_ongoingTasks > 0 || !m_queue.empty())
            m_taskFinnished.wait(lock);
    }

    void ThreadPool::ThreadFunc()
    {
        while (true)
        {
            std::function<void(void)> func;
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                while (!m_stopThreads && m_queue.empty())
                    m_taskAdded.wait(lock);
                
                if (m_stopThreads)
                    break;

                m_ongoingTasks++;
                func = m_queue.front();
                m_queue.pop();
            }
            
            func(); // execute function

            m_ongoingTasks--;
            m_taskFinnished.notify_all();
        }
    }

}

#endif