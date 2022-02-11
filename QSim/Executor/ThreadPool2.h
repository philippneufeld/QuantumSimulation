// Philipp Neufeld, 2021-2022

#ifndef QSim_Executor_ThreadPool_H_
#define QSim_Executor_ThreadPool_H_

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include <memory>

#include <type_traits>
#include <future>

#include "../Platform.h"
#include "Executor.h"

namespace QSim
{

    class ThreadPool
    {
    public:
        ThreadPool();
        ThreadPool(std::size_t threadCnt);
        ~ThreadPool();

        template<typename Func, typename RType=std::invoke_result_t<Func>>
        std::future<RType> Submit(Func& func);

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

    template<typename Func, typename RType>
    std::future<RType> ThreadPool::Submit(Func& func)
    {
        // create task and retrieve the future object
        std::packaged_task<RType()> task(func);
        auto future = task.get_future();

        // std::packaged_task is only movable (not copyable) and
        // std::function requires a copyable functor
        // workaround by using a shared pointer
        auto pTask = std::make_shared<decltype(task)>(std::move(task));

        // add task to queue
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_queue.push([=](){ (*pTask)(); });
        }

        m_taskAdded.notify_one();
        return future;
    }


}

#endif