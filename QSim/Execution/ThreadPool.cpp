// Philipp Neufeld, 2021-2022

#include "ThreadPool.h"

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
        std::unique_lock<std::mutex> lock(m_mutex);
        m_taskFinished.wait(lock, [&](){ return m_ongoingTasks == 0 && m_queue.empty(); }); 

        // stop all threads and wait until they are finished
        m_stopThreads = true;
        lock.unlock();

        m_taskAdded.notify_all(); // wake up all threads

        // wait until all threads have terminated
        for (auto& t: m_threads)
            t.join();
    }

    void ThreadPool::WaitUntilReadyForTask()
    {
        // wait until no new tasks are in the queue and not
        // all threads are occupied
        std::unique_lock<std::mutex> lock(m_mutex);
        m_taskFinished.wait(lock, [&](){ return m_ongoingTasks < m_threads.size() && m_queue.empty(); });
    }
    
    
    void ThreadPool::WaitUntilFinished()
    {
        // wait until all tasks are done
        std::unique_lock<std::mutex> lock(m_mutex);
        m_taskFinished.wait(lock, [&](){ return m_ongoingTasks == 0 && m_queue.empty(); });
    }

    void ThreadPool::WorkerThread()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        while (!m_stopThreads)
        {
            // wait until a new task is added (signalled by the m_taskAdded event)
            m_taskAdded.wait(lock, [&](){ return m_stopThreads || !m_queue.empty(); });
            
            // check thread exit condition
            if (m_stopThreads)
                break;

            // a new task is awailable
            // pop it from the queue and increment the ongoing task counter
            m_ongoingTasks++;
            std::function<void(void)> task = m_queue.front();
            m_queue.pop();

            lock.unlock();
            std::invoke(task); // execute task
            lock.lock();

            // decrement task counter and notify all threads waiting
            // on the m_taskFinished event
            m_ongoingTasks--;
            
            // no need to hold the lock for the signalling of the condition variable
            lock.unlock();
            m_taskFinished.notify_all();
            lock.lock();
        }
    }

    unsigned int ThreadPool::GetAvailableThreads() const 
    { 
        std::unique_lock<std::mutex> lock(m_mutex);
        std::size_t queued = m_queue.size();
        
        std::size_t idleThreads = m_threads.size() - m_ongoingTasks;
        
        if (queued >= idleThreads)
            return 0;
        else
            return idleThreads - queued;
    }

}
