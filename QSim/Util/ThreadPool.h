// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_ThreadPool_H_
#define QSim_Util_ThreadPool_H_

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <atomic>
#include <iostream>

#include "../Platform.h"

namespace QSim
{

    class ThreadPool
    {
    public:
        ThreadPool();
        ThreadPool(std::size_t threadCnt);
        ~ThreadPool();

        void AddTask(const std::function<void(void)>& task);
        void WaitUntilFinnished();

        template<typename Lambda, typename InputIt>
        std::vector<decltype(std::declval<Lambda>()(*std::declval<InputIt>()))> 
            Map(Lambda func, InputIt param_begin, InputIt param_end);

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
        : ThreadPool(std::thread::hardware_concurrency()) { }
    
    ThreadPool::ThreadPool(std::size_t threadCnt)
        : m_stopThreads(false), m_ongoingTasks(0)
    {
        threadCnt = threadCnt > 0 ? threadCnt : 1;
        for (unsigned int i = 0; i < threadCnt; i++)
            m_threads.push_back(std::thread([this](){ ThreadFunc(); }));
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

    template<typename Lambda, typename InputIt>
    std::vector<decltype(std::declval<Lambda>()(*std::declval<InputIt>()))> 
        ThreadPool::Map(Lambda func, InputIt param_begin, InputIt param_end)
    {
        using Ty = decltype(std::declval<Lambda>()(*std::declval<InputIt>()));
        auto diff = param_end - param_begin;
        if (diff <= 0)
            return std::vector<Ty>();
        
        std::vector<Ty> results(diff);
        auto it = param_begin;
        for (std::size_t i = 0; it < param_end; i++, it++)
        {
            auto param = *it;
            AddTask([&, i, param]() { results[i] = func(param); });
        }
        WaitUntilFinnished();

        return results;
    }

}

#endif