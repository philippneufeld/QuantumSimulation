// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_Functor_H_
#define QSim_Math_Functor_H_

#include <type_traits>

namespace QSim
{

#ifndef QSIM_DEFINE_HAS_MEMBER_FUNCTION
#define QSIM_DEFINE_HAS_MEMBER_FUNCTION(funcname)\
    namespace Internal\
    {\
        template<typename... Args>\
        class TIsCallableHelper_##funcname\
        {\
        public:\
            template<typename Func, typename=std::void_t<decltype(std::declval<Func>().funcname (std::declval<Args>()...))>>\
            static std::true_type Test(int);\
            template<typename Func, typename=void>\
            static std::false_type Test(long);\
        };\
    }\
    template<typename Ty, typename... Args>\
    class TIsCallable_##funcname : \
        public decltype(Internal::TIsCallableHelper_##funcname <Args...>::template Test<Ty>(0)) {};\
    template<typename Ty, typename... Args> \
    constexpr bool TIsCallable_##funcname##_v = TIsCallable_##funcname <Ty, Args...>::value;
#endif

    template<typename Func, typename... Args>
    class Functor
    {
    public:

    };
}

#endif
