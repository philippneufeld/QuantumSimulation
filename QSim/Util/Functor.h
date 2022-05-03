// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_Functor_H_
#define QSim_Util_Functor_H_

#include <tuple>
#include <type_traits>

namespace QSim
{

    // makro that defines a type trait that can be used to
    // check the existance of a given member functions
#ifndef QSIM_DEFINE_HAS_MEMBER_FUNCTION
#define QSIM_DEFINE_HAS_MEMBER_FUNCTION(traitname, funcname)\
    namespace Internal\
    {\
        template<typename... Args>\
        class traitname##Helper\
        {\
        public:\
            template<typename Func, typename=std::void_t<decltype(std::declval<Func>().funcname (std::declval<Args>()...))>>\
            static std::true_type Test(int);\
            template<typename Func, typename=void>\
            static std::false_type Test(long);\
        };\
    }\
    template<typename Ty, typename... Args>\
    class traitname : \
        public decltype(Internal:: traitname##Helper <Args...>::template Test<Ty>(0)) {};\
    template<typename Ty, typename... Args> \
    constexpr bool traitname##_v = traitname <Ty, Args...>::value;
#endif

    QSIM_DEFINE_HAS_MEMBER_FUNCTION(THasInvocationOperator, operator())

    template<typename R, typename... Args>
    class TFunctionPrototype
    {
    public:
        using RType = R;
        using ArgTypes = std::tuple<Args...>;
    };

    namespace Internal
    {
        template<typename Func>
        class TFunctionTraitsHelper;
        template<typename R, typename... Args>
        class TFunctionTraitsHelper<R(Args...)>
        {
        public:
            using Prototype = TFunctionPrototype<R, Args...>;
            using RType = typename Prototype::RType;
            using ArgTypes = typename Prototype::ArgTypes;
        };
        template<typename R, typename... Args>
        class TFunctionTraitsHelper<R(*)(Args...)> : public TFunctionTraitsHelper<R(Args...)> {};
        template<typename R, typename... Args>
        class TFunctionTraitsHelper<R(&)(Args...)> : public TFunctionTraitsHelper<R(Args...)> {};
        template<typename C, typename R, typename... Args>
        class TFunctionTraitsHelper<R(C::*)(Args...)> : public TFunctionTraitsHelper<R(Args...)> {};
        template<typename C, typename R, typename... Args>
        class TFunctionTraitsHelper<R(C::*)(Args...) const> : public TFunctionTraitsHelper<R(Args...)> {};
    }

    template<typename Func, typename=void>
    class TFunctionTraits : public Internal::TFunctionTraitsHelper<Func> {};
    template<typename Func>
    class TFunctionTraits<Func, std::enable_if_t<THasInvocationOperator_v<Func>, void>> 
        : public Internal::TFunctionTraitsHelper<decltype(&Func::operator())> {};

    //
    // Functor
    //

    template<typename Func, typename R, typename ArgTuple>
    class TFunctor;
    template<typename Func, typename R, typename... Args>
    class TFunctor<Func, R, std::tuple<Args...>>
    {
        static_assert(std::is_invocable_r_v<R, Func, Args...>, "Func must be invokable");
    public:
        using Prototype = TFunctionPrototype<R, Args...>;
        
        TFunctor(Func& func) : m_func(func) {}

        // copy operations
        TFunctor(const TFunctor& rhs) = default;
        TFunctor(TFunctor&& rhs) = default; 
        TFunctor& operator=(const TFunctor& rhs) = default;
        TFunctor& operator=(TFunctor&& rhs) = default;

        auto Invoke(Args... args) { return static_cast<R>(std::invoke(m_func, std::forward<Args>(args)...)); }
        auto operator()(Args... args) { return Invoke(std::forward<Args>(args)...); }

    protected:
        Func& m_func;
    };

    template<typename Func>
    auto CreateFunctor(Func& func)
    {
        using RType = typename TFunctionTraits<Func>::RType;
        using ArgTypes = typename TFunctionTraits<Func>::ArgTypes;
        return TFunctor<Func, RType, ArgTypes>(func);
    }
    
}

#endif
