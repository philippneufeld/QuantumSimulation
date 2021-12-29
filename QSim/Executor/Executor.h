// Philipp Neufeld, 2021-2022

#ifndef QSim_Executor_Executor_H_
#define QSim_Executor_Executor_H_

#include <functional>

#include "../Util/CRTP.h"

namespace QSim
{
    // Base Executor CRTP class
    // 
    // Implementations (classes deriving from TExecutor)
    // must implement:
    //  1.) void AddTask(const std::function<void(void)>& task);
    //  2.) void WaitUntilFinnished();
    template<typename T>
    class TExecutor : public TCRTP<T> 
    {
    protected:
        TExecutor() = default;

    public:
        template<typename Lambda, typename InputIt, typename Container>
        void MapNonBlocking(const Lambda& func, Container& dest, InputIt param_begin, InputIt param_end);
        template<typename Lambda, typename InputIt, typename Container>
        void Map(const Lambda& func, Container& dest, InputIt param_begin, InputIt param_end);
    };
    
    template<typename T>
    template<typename Lambda, typename InputIt, typename Container>
    void TExecutor<T>::MapNonBlocking(const Lambda& func, Container& dest, InputIt paramBegin, InputIt paramEnd)
    {
        auto it = paramBegin;
        for (std::size_t i = 0; it < paramEnd; it++, i++)
        {
            auto param = *it;
            (~(*this)).AddTask([&, i, param]() { dest[i] = func(param); });
        }

    }

    template<typename T>
    template<typename Lambda, typename InputIt, typename Container>
    void TExecutor<T>::Map(const Lambda& func, Container& dest, InputIt paramBegin, InputIt paramEnd)
    {
        MapNonBlocking(func, dest, paramBegin, paramEnd);
        (~(*this)).WaitUntilFinnished();
    }

    // Default executor
    // Just calls task directly and blocks until task is finnished
    class DefaultExecutor : public TExecutor<DefaultExecutor>
    {
    public:
        DefaultExecutor() = default;
        void AddTask(const std::function<void(void)>& task) { task(); }
        void WaitUntilFinnished() {}
    };

}

// Include specific implementations
#include "ThreadPool.h"

#endif
