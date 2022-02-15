// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_MatrixTraits_H
#define QSim_Math_MatrixTraits_H

#include <type_traits>
#include <Eigen/Dense>

namespace QSim
{
    

    namespace Internal
    {
        template<typename Ty>
        struct TMatrixEvalTypeHelper
        {
        private:
            template<typename Ty2, typename=std::void_t<std::decay_t<decltype(std::declval<std::decay_t<Ty2>>().eval())>>>
            static auto Detect(int) -> std::decay_t<decltype(std::declval<std::decay_t<Ty2>>().eval())>;
            template<typename Ty2, typename=void>
            static auto Detect(long) -> Ty2;

        public:
            using type = decltype(Detect<Ty>(0));
        };
    }

    template<typename Ty>
    struct TIsMatrix : std::integral_constant<bool, std::is_base_of_v<
        Eigen::EigenBase<typename Internal::TMatrixEvalTypeHelper<Ty>::type>, 
        typename Internal::TMatrixEvalTypeHelper<Ty>::type>> {};
    template<typename Ty>
    constexpr bool TIsMatrix_v = TIsMatrix<Ty>::value;

    template<typename Ty>
    struct TMatrixEvalType
    {
        using type = std::conditional_t<TIsMatrix_v<Ty>, 
            typename Internal::TMatrixEvalTypeHelper<Ty>::type, 
            std::decay_t<Ty>>;
    };
    template<typename Ty>
    using TMatrixEvalType_t = typename TMatrixEvalType<Ty>::type;

    // Element type
    namespace Internal
    {
        template<typename Ty>
        struct TMatrixElementTypeHelper
        {
            using type = std::decay_t<Ty>;
        };

        template<typename Ty, int N, int M>
        struct TMatrixElementTypeHelper<Eigen::Matrix<Ty, N, M>>
        {
            using type = std::decay_t<Ty>;
        };
    }
    template<typename Ty>
    struct TMatrixElementType 
        : Internal::TMatrixElementTypeHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    using TMatrixElementType_t = typename TMatrixElementType<Ty>::type;

    // Class that makes the type and value of the 
    // length of a dx vector (or scalar) easily accessible
    namespace Internal
    {
        template<typename Ty>
        struct TDxLengthHelper
        {
            using type = Ty;
            static Ty Get(Ty dx) { return dx; }
        };

        template<typename Ty, int N, int M>
        struct TDxLengthHelper<Eigen::Matrix<Ty, N, M>>
        {
            using type = Ty;
            static Ty Get(const Eigen::Matrix<Ty, N, M>& dx) { return dx.norm(); }
        };
    }

    template<typename Ty>
    struct TDxLength : public Internal::TDxLengthHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    using TDxLength_t = typename TDxLength<Ty>::type;

    // Matrix rows and cols at compile time
    namespace Internal
    {
        template<typename Ty>
        struct TMatrixRowsAtCompileTimeHelper 
            : public std::integral_constant<int, 1> {};
        template<typename Ty, int N, int M>
        struct TMatrixRowsAtCompileTimeHelper<Eigen::Matrix<Ty, N, M>> 
            : public std::integral_constant<int, N> {};
    
        template<typename Ty>
        struct TMatrixColsAtCompileTimeHelper 
            : public std::integral_constant<int, 1> {};
        template<typename Ty, int N, int M>
        struct TMatrixColsAtCompileTimeHelper<Eigen::Matrix<Ty, N, M>> 
            : public std::integral_constant<int, M> {};

        template<int N, int M>
        struct TMatrixSizeAtCompileTimeHelper
            : public std::integral_constant<int, (N < 0 || M < 0) ? Eigen::Dynamic : N*M> {};
    }

    template<typename Ty>
    struct TMatrixRowsAtCompileTime 
        : public Internal::TMatrixRowsAtCompileTimeHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixRowsAtCompileTime_v = TMatrixRowsAtCompileTime<Ty>::value;

    template<typename Ty>
    struct TMatrixColsAtCompileTime 
        : public Internal::TMatrixColsAtCompileTimeHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixColsAtCompileTime_v = TMatrixColsAtCompileTime<Ty>::value;

    template<typename Ty>
    struct TMatrixSizeAtCompileTime : public Internal::TMatrixSizeAtCompileTimeHelper<
        TMatrixRowsAtCompileTime_v<Ty>, TMatrixColsAtCompileTime_v<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixSizeAtCompileTime_v = TMatrixSizeAtCompileTime<Ty>::value;

}

#endif
