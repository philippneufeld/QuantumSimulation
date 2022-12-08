// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ThreadPool_H_
#define QSim_Execution_ThreadPool_H_

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <memory>
#include <type_traits>
#include <future>

#include "Progress.h"

namespace QSim
{

    class ThreadPool
    {
    public:
        ThreadPool();
        ThreadPool(std::size_t threadCnt);
        ~ThreadPool();

        static Progress CreateProgressTracker(std::size_t cnt) { return Progress(cnt); }

        template<typename Task>
        void Submit(Task&& task);

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::future<RType> SubmitWithFuture(Task&& task);

        void WaitUntilReadyForTask();
        void WaitUntilFinished();

        unsigned int GetAvailableThreads() const;

    private:
        void WorkerThread();

    private:
        std::vector<std::thread> m_threads; 
        std::queue<std::function<void(void)>> m_queue;
        
        // thread syncronization variables
        mutable std::mutex m_mutex;
        std::condition_variable m_taskAdded;
        std::condition_variable m_taskFinished;
        
        // auxilliary state variables
        unsigned int m_ongoingTasks;
        bool m_stopThreads;
    };

    //
    // Template function definitions
    //

    template<typename Task>
    void ThreadPool::Submit(Task&& task)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_queue.push(std::forward<Task>(task));
        lock.unlock();
        m_taskAdded.notify_one();
    }

    template<typename Task, typename RType>
    inline std::future<RType> ThreadPool::SubmitWithFuture(Task&& task)
    {
        // create task and retrieve the future object
        std::packaged_task<RType(void)> packed(std::forward<Task>(task));
        auto future = packed.get_future();

        // std::packaged_task is only movable (not copyable) and
        // std::function requires a copyable functor
        // workaround by using a shared pointer
        auto pTask = std::make_shared<decltype(packed)>(std::move(packed));
        auto func = [=](){ std::invoke(*pTask); };

        Submit(std::move(func));

        return std::move(future);
    }

}

#endif