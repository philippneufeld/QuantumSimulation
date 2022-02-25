// Philipp Neufeld, 2021-2022

#ifndef QSim_Math_MatrixTraits_H
#define QSim_Math_MatrixTraits_H

#include <utility>
#include <type_traits>
#include <Eigen/Dense>

namespace QSim
{
    

    //
    // TMatrixEvalType<Ty> trait
    // Queryies the return type of .eval() if Ty has this method implemented 
    // and otherwise returns std::decay_t<Ty>
    //
    // TIsMatrix<Ty> trait
    // Checks if TMatrixEvalType<Ty> is derived from Eigen::EigenBase
    //
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

    //
    // MatrixEval()
    // Calls the eval() method of the argument if such a method exists
    //
    namespace Internal
    {
        template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
        struct TMatrixEvalHelper
        {
            static_assert(TIsMatrix_v<Ty>, "Ty must be a matrix type");
            static TMatrixEvalType_t<Ty> Eval(Ty&& val) { return val.eval(); }
        };

        template<typename Ty>
        struct TMatrixEvalHelper<Ty, false>
        {
            static_assert(!TIsMatrix_v<Ty>, "Ty mut not be a matrix type");
            static TMatrixEvalType_t<Ty> Eval(Ty&& val) { return val; }
        };
    }
    template<typename Ty>
    auto MatrixEval(Ty&& val) 
    { 
        return Internal::TMatrixEvalHelper<Ty>::Eval(std::forward<Ty>(val)); 
    }
    

    //
    // TMatrixElementType<Ty> trait
    // Queries the type of a matrix element of Ty (if Ty is a matrix type)
    // or Ty itself
    //
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

    //
    // TMatrixCompileTimeRows, TMatrixCompileTimeCols, TMatrixCompileTimeSize
    // Traits to query the compile-time dimensions of a matrix (<0 if dynamic)
    //
    namespace Internal
    {
        template<typename Ty>
        struct TMatrixCompileTimeRowsHelper 
            : public std::integral_constant<int, 1> {};
        template<typename Ty, int N, int M>
        struct TMatrixCompileTimeRowsHelper<Eigen::Matrix<Ty, N, M>> 
            : public std::integral_constant<int, N> {};
    
        template<typename Ty>
        struct TMatrixCompileTimeColsHelper 
            : public std::integral_constant<int, 1> {};
        template<typename Ty, int N, int M>
        struct TMatrixCompileTimeColsHelper<Eigen::Matrix<Ty, N, M>> 
            : public std::integral_constant<int, M> {};

        template<int N, int M>
        struct TMatrixCompileTimeSizeHelper
            : public std::integral_constant<int, (N < 0 || M < 0) ? Eigen::Dynamic : N*M> {};
    }

    template<typename Ty>
    struct TMatrixCompileTimeRows 
        : public Internal::TMatrixCompileTimeRowsHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixCompileTimeRows_v = TMatrixCompileTimeRows<Ty>::value;

    template<typename Ty>
    struct TMatrixCompileTimeCols 
        : public Internal::TMatrixCompileTimeColsHelper<TMatrixEvalType_t<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixCompileTimeCols_v = TMatrixCompileTimeCols<Ty>::value;

    template<typename Ty>
    struct TMatrixCompileTimeSize : public Internal::TMatrixCompileTimeSizeHelper<
        TMatrixCompileTimeRows_v<Ty>, TMatrixCompileTimeCols_v<Ty>> {};
    template<typename Ty>
    constexpr int TMatrixCompileTimeSize_v = TMatrixCompileTimeSize<Ty>::value;

    //
    // TMatrixAddResult trait
    // gets the result type of an addition expression
    //
    template<typename XTy1, typename... XTys>
    struct TMatrixAddResult
    {
        using type = TMatrixEvalType_t<decltype(std::declval<XTy1>() + std::declval<typename TMatrixAddResult<XTys...>::type>())>;
    };
    template<typename XTy1>
    struct TMatrixAddResult<XTy1> : public TMatrixEvalType<XTy1> {};
    template<typename XTyA, typename... XTys>
    using TMatrixAddResult_t = typename TMatrixAddResult<XTyA, XTys...>::type;

    template<typename XTy1, typename... XTys>
    struct TMatrixAddResultFP
    {
        using type = std::conditional_t<std::is_integral_v<TMatrixAddResult_t<XTy1, XTys...>>, 
            double, TMatrixAddResult_t<XTy1, XTys...>>;
    };
    template<typename XTy1, typename... XTys>
    using TMatrixAddResultFP_t = typename TMatrixAddResultFP<XTy1, XTys...>::type;

    //
    // TMatrixMulResult trait
    // gets the result type of a multiplication expression
    //
    template<typename XTy1, typename... XTys>
    struct TMatrixMulResult
    {
        using type = TMatrixEvalType_t<decltype(std::declval<XTy1>() * 
            std::declval<typename TMatrixMulResult<XTys...>::type>())>;
    };
    template<typename XTy1>
    struct TMatrixMulResult<XTy1> : public TMatrixEvalType<XTy1> {};
    template<typename XTyA, typename... XTys>
    using TMatrixMulResult_t = typename TMatrixMulResult<XTyA, XTys...>::type;

    template<typename XTy1, typename... XTys>
    struct TMatrixMulResultFP
    {
        using type = std::conditional_t<std::is_integral_v<TMatrixMulResult_t<XTy1, XTys...>>, 
            double, TMatrixMulResult_t<XTy1, XTys...>>;
    };
    template<typename XTy1, typename... XTys>
    using TMatrixMulResultFP_t = typename TMatrixMulResultFP<XTy1, XTys...>::type;

    //
    // TMatrixNorm<Ty>
    // Provides a trait used to call the .norm() method of a matrix
    // (or std::abs value in case of scalars)
    //
    template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
    struct TMatrixNorm
    {
        static_assert(TIsMatrix_v<Ty>);
        using type = TMatrixElementType_t<Ty>;
        static auto Get(const Ty& dx) { return dx.norm(); } 
    };
    template<typename Ty>
    struct TMatrixNorm<Ty, false>
    {
        static_assert(!TIsMatrix_v<Ty>);
        using type = TMatrixEvalType_t<Ty>;
        static auto Get(Ty dx) { return std::abs(dx); }
    };
    template<typename Ty>
    using TMatrixNorm_t = typename TMatrixNorm<Ty>::type;

    //
    // TMatrixOnesLike<Ty>
    // Implements a unifying interface for creating scalar constants 
    // and matrix constants in which all elements have the value 1
    //
    template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
    struct TMatrixOnesLike
    {
        static_assert(TIsMatrix_v<Ty>);
        static auto Get(const Ty& like) 
        { 
            return TMatrixEvalType_t<Ty>::Ones(like.rows(), like.cols()); 
        }
    };
    template<typename Ty>
    struct TMatrixOnesLike<Ty, false>
    {
        static_assert(!TIsMatrix_v<Ty>);
        static auto Get(const Ty&) { return Ty{1}; }
    };

    //
    // TMatrixCwiseAbs
    // Implements a unifying interface for calculating the 
    // componentwise absolute value (works for matrices and scalars)
    //
    template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
    struct TMatrixCwiseAbs
    {
        static_assert(TIsMatrix_v<Ty>);
        static auto Get(const Ty& mat) { return mat.cwiseAbs().eval(); }
    };
    template<typename Ty>
    struct TMatrixCwiseAbs<Ty, false>
    {
        static_assert(!TIsMatrix_v<Ty>);
        static auto Get(const Ty& val) { return std::abs(val); }
    };

    //
    // TMatrixCwiseAbsMax
    // Implements a unifying interface for calculating the 
    // componentwise absolute maximum value (works for matrices and scalars)
    //
    template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
    struct TMatrixCwiseAbsMax
    {
        static_assert(TIsMatrix_v<Ty>);
        static auto Get(const Ty& mat1, const Ty& mat2) 
        { 
            return (mat1.cwiseAbs()).cwiseMax((mat2.cwiseAbs())).eval(); 
        }
    };
    template<typename Ty>
    struct TMatrixCwiseAbsMax<Ty, false>
    {
        static_assert(!TIsMatrix_v<Ty>);
        static auto Get(const Ty& val1, const Ty& val2) 
        { 
            return std::max(std::abs(val1), std::abs(val2)); 
        }
    };


    //
    // TMatrixAnyCwiseLess
    // Implements a unifying interface for checking if any component 
    // of the first argument is less than the right argument
    //
    template<typename Ty, bool isMat=TIsMatrix_v<Ty>>
    struct TMatrixAnyCwiseLess
    {
        static bool Get(const Ty& mat1, const Ty& mat2) 
        { 
            return (mat1.array() < mat2.array()).any(); 
        }
    };
    template<typename Ty>
    struct TMatrixAnyCwiseLess<Ty, false>
    {
        static bool Get(const Ty& val1, const Ty& val2) 
        { 
            return val1 < val2; 
        }
    };

}

#endif
