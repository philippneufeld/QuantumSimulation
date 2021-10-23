// Philipp Neufeld, 2021

#ifndef QSim_NLevelSystem_H_
#define QSim_NLevelSystem_H_

#include <cstdint>
#include <algorithm>
#include <vector>
#include <cassert>

#include "Matrix.h"

namespace QSim
{

    template<std::size_t N, typename Ty = double>
    class TNLevelSystem
    {
        struct Transition 
        {
            std::size_t lvl1;
            std::size_t lvl2;
            Ty rabi;
        };

    public:
        // constructors
        TNLevelSystem(const Ty (&levels)[N]) : TNLevelSystem(static_cast<const Ty*>(levels)) { } 
        TNLevelSystem(const std::initializer_list<Ty> levels) : TNLevelSystem(levels.begin()) { assert(levels.size() == N); }
        TNLevelSystem(const Ty* levels); 
        ~TNLevelSystem();

        // copy operations
        TNLevelSystem(const TNLevelSystem& rhs) = default;
        TNLevelSystem& operator=(const TNLevelSystem& rhs) = default;

    private:
        Ty m_levels[N];
    };

    template<std::size_t N, typename Ty>
    TNLevelSystem<N, Ty>::TNLevelSystem(const Ty* levels)
        : m_levels{}
    {
        std::copy(levels, levels + N, m_levels);
    }

    template<std::size_t N, typename Ty>
    TNLevelSystem<N, Ty>::~TNLevelSystem()
    {

    }

}

#endif
