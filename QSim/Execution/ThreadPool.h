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

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::future<RType> Submit(Task&& task);

    private:
        template<typename Task, typename RType>
        static std::pair<std::function<void(void)>, std::future<RType>> 
            PackageTask(Task&& task);

        void EnqueueTask(const std::function<void(void)>& func);
        void EnqueueTask(std::function<void(void)>&& func);
        
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

    template<typename Task, typename RType>
    inline std::future<RType> ThreadPool::Submit(Task&& task)
    {
        auto [func, future] = PackageTask<Task, RType>(std::forward<Task&&>(task));
        EnqueueTask(func);
        return std::move(future);
    }

    template<typename Task, typename RType>
    inline std::pair<std::function<void(void)>, std::future<RType>> 
        ThreadPool::PackageTask(Task&& task)
    {
        // create task and retrieve the future object
        std::packaged_task<RType(void)> packed(std::forward<Task&&>(task));
        auto future = packed.get_future();

        // std::packaged_task is only movable (not copyable) and
        // std::function requires a copyable functor
        // workaround by using a shared pointer
        auto pTask = std::make_shared<decltype(packed)>(std::move(packed));
        auto func = [=](){ std::invoke(*pTask); };

        return std::make_pair(std::move(func), std::move(future));
    }

}

#endif