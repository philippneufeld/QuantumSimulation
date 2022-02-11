// Philipp Neufeld, 2021-2022

#include "ThreadPool2.h"

namespace QSim
{
    using namespace std::chrono_literals;

    ThreadPool::ThreadPool()
        : ThreadPool(std::thread::hardware_concurrency()) { }
    
    ThreadPool::ThreadPool(std::size_t threadCnt)
        : m_stopThreads(false), m_ongoingTasks(0)
    {
        threadCnt = threadCnt > 0 ? threadCnt : 1;
        for (unsigned int i = 0; i < threadCnt; i++)
            m_threads.push_back(std::thread([this](){ WorkerThread(); }));
    }

    ThreadPool::~ThreadPool()
    {   
        // wait until all tasks are done
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            while (m_ongoingTasks > 0 || !m_queue.empty())
                m_taskFinished.wait_for(lock, 10ms); // wake up every 10ms if nothing happened
        }

        // stop all threads and wait until they are finished
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_stopThreads = true;
        }
        m_taskAdded.notify_all(); // wake up all threads

        // wait until all threads have terminated
        for (auto& t: m_threads)
            t.join();
    }

    void ThreadPool::WorkerThread()
    {
        while (true)
        {
            // create empty task object
            std::function<void(void)> task;

            // extract task from queue whenever one is available
            {
                // make sure the queue is not accessed by two threads simultaneously
                std::unique_lock<std::mutex> lock(m_mutex);

                // wait until a new task is added (signalled by the m_taskAdded event)
                // spuriously wake up every 10ms
                while (!m_stopThreads && m_queue.empty())
                    m_taskAdded.wait_for(lock, 10ms); // wake up every 10ms if nothing happened
                
                // check thread exit condition
                if (m_stopThreads)
                    break;

                // a new task is awailable
                // pop it from the queue and increment the ongoing task counter
                m_ongoingTasks++;
                task = m_queue.front();
                m_queue.pop();
            }
            
            task(); // execute task

            // decrement task counter and notify all threads waiting
            // on the m_taskFinished event
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_ongoingTasks--;
            }
            m_taskFinished.notify_all();
        }
    }

}
