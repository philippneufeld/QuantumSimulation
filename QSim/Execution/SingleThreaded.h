// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_SingleThreaded_H_
#define QSim_Execution_SingleThreaded_H_

#include <functional>
#include <type_traits>
#include <future>

#include "Progress.h"

namespace QSim
{

    class SingleThreaded
    {
    public:
        SingleThreaded() = default;

        static ProgressST CreateProgressTracker(std::size_t cnt) { return ProgressST(cnt); }

        template<typename Task>
        void Submit(Task&& task)
        {
            std::invoke(task);
        }

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::enable_if_t<!std::is_same_v<RType, void>, std::future<RType>> 
            SubmitWithFuture(Task&& task)
        {
            std::promise<RType> res;
            res.set_value(std::invoke(task));
            return res.get_future();
        }

        template<typename Task, typename RType=std::invoke_result_t<Task>>
        std::enable_if_t<std::is_same_v<RType, void>, std::future<void>> 
            SubmitWithFuture(Task&& task)
        {
            std::promise<RType> res;
            std::invoke(task);
            res.set_value();
            return res.get_future();
        }

        void WaitUntilFinished() {}
    };

}

#endif