// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_CRTP_H_
#define QSim_Util_CRTP_H_

#include <type_traits>

namespace QSim
{
    // CRTP
    template <typename T>
    class TCRTP
    {
    protected:
        ~TCRTP() = default;
        
    public:
        // CRTP conversion operators to the actual vector type (type-safe downcast)
        constexpr T &operator~() noexcept { return static_cast<T &>(*this); }
        constexpr const T &operator~() const noexcept { return static_cast<const T &>(*this); }
    };

    // Trait that checks if type is a crtp type
    template <typename T, template <typename> class CRTP = TCRTP, typename = void>
    struct TIsCRTP : std::false_type {};
    template <typename T, template <typename> class CRTP>
    struct TIsCRTP<T, CRTP, std::enable_if_t<
        std::is_base_of<TCRTP<std::decay_t<decltype(~std::declval<T>())>>, std::decay_t<T>>::value && 
        std::is_base_of<CRTP<std::decay_t<decltype(~std::declval<T>())>>, std::decay_t<T>>::value>> 
        : std::true_type {};

    template <typename MT>
    constexpr bool TIsCRTP_v = TIsCRTP<MT>::value;

    // Get the inner CRTP type
    template <typename T, typename = void>
    struct TDecayCRTP;
    template <typename T>
    struct TDecayCRTP<T, std::enable_if_t<TIsCRTP_v<T>>>
    {
        using type = std::decay_t<decltype(~std::declval<T>())>;
    };

    template <typename T>
    using TDecayCRTP_t = typename TDecayCRTP<T>::type;

    // Get the inner CRTP type (also works for non-CRTP types and uses std::decay on them)
    template <typename T, typename = void>
    struct TDecayCRTPGeneral 
    {
        using type = std::decay_t<T>;
    };
    template <typename T>
    struct TDecayCRTPGeneral<T, std::enable_if_t<TIsCRTP_v<T>>>
    {
        using type = std::decay_t<decltype(~std::declval<T>())>;
    };

    template <typename T>
    using TDecayCRTPGeneral_t = typename TDecayCRTPGeneral<T>::type;
}

#endif
