// Philipp Neufeld, 2021

#ifndef QSim_Util_ConstList_H_
#define QSim_Util_ConstList_H_

#include "../Platform.h"

namespace QSim
{
    template<typename Ty, Ty... Cs>
    struct TConstList;

    template<typename CL>
    struct TConstListSizeof;
    template<typename Ty, Ty... Cs>
    struct TConstListSizeof<TConstList<Ty, Cs...>>
        : public std::integral_constant<std::size_t, sizeof...(Cs)> {};
    template<typename CL>
    constexpr auto TConstListSizeof_v = TConstListSizeof<CL>::value;

    template<std::size_t idx, typename CL>
    struct TConstListGet;
    template<std::size_t idx, typename Ty, Ty C, Ty... Cs>
    struct TConstListGet<idx, TConstList<Ty, C, Cs...>>
        : public TConstListGet<idx-1, TConstList<Ty, Cs...>> {};
    template<typename Ty, Ty C, Ty... Cs>
    struct TConstListGet<0, TConstList<Ty, C, Cs...>>
        : public std::integral_constant<Ty, C> {};
    template<std::size_t idx, typename CL>
    constexpr auto TConstListGet_v = TConstListGet<idx, CL>::value;
    
    template<typename CL>
    struct TConstListFront : public TConstListGet<0, CL> {};
    template<typename CL>
    constexpr auto TConstListFront_v = TConstListFront<CL>::value;

    template<typename CL>
    struct TConstListBack : public TConstListGet<TConstListSizeof_v<CL> - 1, CL> {};
    template<typename CL>
    constexpr auto TConstListBack_v = TConstListBack<CL>::value;

    template<typename CL1, typename CL2>
    struct TConstListConcatenate;
    template<typename Ty, Ty... Cs1, Ty... Cs2>
    struct TConstListConcatenate<TConstList<Ty, Cs1...>, TConstList<Ty, Cs2...>>
    {
        using type = TConstList<Ty, Cs1..., Cs2...>;
    };
    template<typename CL1, typename CL2>
    using TConstListConcatenate_t = typename TConstListConcatenate<CL1, CL2>::type;

    template<typename Ty, Ty C, typename CL>
    struct TConstListPrepend : public TConstListConcatenate<TConstList<Ty, C>, CL> {};
    template<typename Ty, Ty C, typename CL>
    using TConstListPrepend_t = typename TConstListPrepend<Ty, C, CL>::type;

    template<typename Ty, Ty C, typename CL>
    struct TConstListAppend : public TConstListConcatenate<CL, TConstList<Ty, C>> {};
    template<typename Ty, Ty C, typename CL>
    using TConstListAppend_t = typename TConstListAppend<Ty, C, CL>::type;

    template<std::size_t idx, typename CL>
    struct TConstListErase;
    template<std::size_t idx, typename Ty, Ty C, Ty... Cs>
    struct TConstListErase<idx, TConstList<Ty, C, Cs...>>
        : public TConstListPrepend<Ty, C, typename TConstListErase<idx-1, TConstList<Ty, Cs...>>::type> {};
    template<typename Ty, Ty C, Ty... Cs>
    struct TConstListErase<0, TConstList<Ty, C, Cs...>>
    {
        using type = TConstList<Ty, Cs...>;
    };
    template<std::size_t idx, typename CL>
    using TConstListErase_t = typename TConstListErase<idx, CL>::type;

    template<typename CL>
    struct TConstListPopFront : public TConstListErase<0, CL> {};
    template<typename CL>
    using TConstListPopFront_t = typename TConstListPopFront<CL>::type;

    template<typename CL>
    struct TConstListPopBack : public TConstListErase<TConstListSizeof_v<CL>-1, CL> {};
    template<typename CL>
    using TConstListPopBack_t = typename TConstListPopBack<CL>::type;

    template<std::size_t idx, typename Ty, Ty C, typename CL>
    struct TConstListInsert 
        : public TConstListPrepend<Ty, TConstListFront_v<CL>, typename TConstListInsert<idx-1, Ty, C, TConstListPopFront_t<CL>>::type> {};
    template<typename Ty, Ty C, typename CL>
    struct TConstListInsert<0, Ty, C, CL>
        : public TConstListPrepend<Ty, C, CL> {};
    template<std::size_t idx, typename Ty, Ty C, typename CL>
    using TConstListInsert_t = typename TConstListInsert<idx, Ty, C, CL>::type;
}

#endif
