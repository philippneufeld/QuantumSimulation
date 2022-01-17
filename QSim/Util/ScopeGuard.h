// Copyright 2021-2022, Philipp Neufeld

#ifndef QSim_Util_ScopeGuard_H_
#define QSim_Util_ScopeGuard_H_

#include <type_traits>
#include <utility>
#include <functional>

namespace QSim
{
        
    template<typename Lambda>
    class TScopeGuard
    {
    public:
        TScopeGuard(Lambda onScopeExit)
            : m_onScopeExit(onScopeExit), m_dismissed(false) { }
        ~TScopeGuard() { if (!m_dismissed) m_onScopeExit(); }

        TScopeGuard(const TScopeGuard&) = delete;
        TScopeGuard& operator=(const TScopeGuard&) = delete;

        TScopeGuard(TScopeGuard&& other)
            : m_onScopeExit(std::move(other.m_onScopeExit)),
            m_dismissed(false) {
            other.Dismiss();
        }
        TScopeGuard& operator=(TScopeGuard&& other) = delete;

        void Dismiss() { m_dismissed = true; }

    private:
        bool m_dismissed;
        Lambda m_onScopeExit;
    };

    // Caution: When using this function ALWAYS save the return value into 
    // a variable. Otherwise, the scope of the ScopeGuard ends just after the 
    // function call because the ScopeGuard object is discarded!
    // Proper usage:
    // auto scopeGuard = CreateScopeGuard(cleanup_func);
    template<typename Lambda>
    auto CreateScopeGuard(Lambda onScopeExit)
    {
        return TScopeGuard<Lambda>{onScopeExit};
    }

}

#endif // !QUANTUM_Util_ScopeGuard_H_
