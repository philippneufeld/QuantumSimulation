// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_SingleThreaded_H_
#define QSim_Execution_SingleThreaded_H_

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

    class SingleThreaded
    {
    public:
        SingleThreaded() = default;

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::enable_if_t<!std::is_same_v<RType, void>, std::future<RType>> 
            Submit(Task task, std::shared_ptr<Progress> pProgress = nullptr)
        {
            std::promise<RType> res;
            res.set_value(std::invoke(task));
            if (pProgress)
                pProgress->IncrementCount();
            return res.get_future();
        }

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::enable_if_t<std::is_same_v<RType, void>, std::future<void>> 
            Submit(Task task, std::shared_ptr<Progress> pProgress = nullptr)
        {
            std::promise<RType> res;
            std::invoke(task);
            res.set_value();
            if (pProgress)
                pProgress->IncrementCount();
            return res.get_future();
        }
    };

}

#endif