// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_ConstexprFor_H_
#define QSim_Util_ConstexprFor_H_

#include "../Platform.h"
#include <type_traits>

namespace QSim
{
    // constexpr for
    // see https://artificial-mind.net/blog/2020/10/31/constexpr-for
    namespace Internal
    {
        template <typename ItType, ItType Start, ItType End, ItType Inc, class F, bool end>
        struct TConstexprForHelper
        {
            QSim_ALWAYS_INLINE constexpr static void Exec(F &&f)
            {
                f(std::integral_constant<ItType, Start>());
                TConstexprForHelper<ItType, Start + Inc, End, Inc, F, std::less<ItType>{}(Start + Inc, End)>::Exec(std::forward<F>(f));
            }
        };

        template <typename ItType, ItType Start, ItType End, ItType Inc, class F>
        struct TConstexprForHelper<ItType, Start, End, Inc, F, false>
        {
            QSim_ALWAYS_INLINE constexpr static void Exec(F &&f) {}
        };
    }

    template <typename ItType, ItType Start, ItType End, ItType Inc, class F>
    QSim_ALWAYS_INLINE constexpr void ConstexprFor(F &&f)
    {
        Internal::TConstexprForHelper<ItType, Start, End, Inc, F, std::less<ItType>{}(Start, End)>::Exec(std::forward<F>(f));
    }
}

#endif
