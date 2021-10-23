// Philipp Neufeld, 2021

#ifndef QSim_NLevel_Matrix_H_
#define QSim_NLevel_Matrix_H_

#include <cstdint>
#include <type_traits>
#include <algorithm>
#include <complex>

namespace QSim
{

    namespace Internal
    {
        template<typename T>
        struct TIsComplex_impl : std::false_type {};
        template<typename T>
        struct TIsComplex_impl<std::complex<T>> {};
    }

    template<typename T>
    struct TIsComplex : Internal::TIsComplex_impl<std::decay_t<T>> {};
    template<typename T>
    constexpr auto TIsComplex_v = TIsComplex<T>::value;


    template<typename Ty, std::size_t N, std::size_t M>
    class TMatrix
    {
        template<typename U, std::size_t L, std::size_t O>
        friend class TMatrix;

        template<typename Dummy>
        constexpr static bool TIsVector_v = ((N==1) || (M==1));

        template<typename U, typename Dummy>
        using TEnableIfComplex_t = std::enable_if_t<std::is_same<Dummy, Dummy>::value && TIsComplex_v<Ty>, U>;
        template<typename U, typename Dummy>
        using TDisableIfComplex_t = std::enable_if_t<std::is_same<Dummy, Dummy>::value && !TIsComplex_v<Ty>, U>;
    public:
        TMatrix() : m_data{} {} // m_data is initialized to its default value this way
        template<typename... Ts, typename=std::enable_if_t<sizeof...(Ts) == N*M>>
        TMatrix(const Ts&... vals) : m_data{ static_cast<Ty>(std::forward<const Ts&>(vals))... } {}
        ~TMatrix() = default;

        TMatrix(const TMatrix& rhs) { std::copy(rhs.m_data, rhs.m_data + N*M, m_data); } 
        TMatrix& operator=(const TMatrix& rhs) { std::copy(rhs.m_data, rhs.m_data + N*M, m_data); return *this; } 
       
        Ty& operator()(std::size_t i, std::size_t j) { return m_data[i*M+j]; }
        const Ty& operator()(std::size_t i, std::size_t j) const { return m_data[i*M+j]; }

        template<typename Dummy=void, typename=std::enable_if_t<TIsVector_v<Dummy>>>
        const Ty& operator()(std::size_t i) const { return m_data[i]; } 
        template<typename Dummy=void, typename=std::enable_if_t<TIsVector_v<Dummy>>>
        Ty& operator()(std::size_t i) { return m_data[i]; }

        TMatrix& operator+=(const TMatrix& rhs);
        TMatrix& operator-=(const TMatrix& rhs);
        TMatrix& operator*=(const TMatrix<Ty, M, M>& rhs) { return this->operator=((*this) * rhs); }

        friend TMatrix operator+(const TMatrix& lhs, const TMatrix& rhs) { TMatrix tmp(lhs); tmp += rhs; return tmp; }
        friend TMatrix operator-(const TMatrix& lhs, const TMatrix& rhs) { TMatrix tmp(lhs); tmp -= rhs; return tmp; }
        //template<std::size_t L>
        //friend TMatrix<Ty, N, L> operator*(const TMatrix<Ty, N, M>& lhs, const TMatrix<Ty, M, L>& rhs);

        TMatrix<Ty, M, N> transpose() const;
        template<typename U=void> 
        TEnableIfComplex_t<TMatrix<Ty, M, N>, U> adjoint() const;
        template<typename U=void> 
        TDisableIfComplex_t<TMatrix<Ty, M, N>, U> adjoint() const { return transpose(); }

        void SetZero() { std::fill(m_data, m_data + N*M, static_cast<Ty>(0)); }

    private:
        Ty m_data[N*M];
    };

    template<typename Ty, std::size_t N, std::size_t M>
    TMatrix<Ty, N, M>& TMatrix<Ty, N, M>::operator+=(const TMatrix& rhs)
    {
        for (std::size_t i = 0; i < N*M; i++)
            m_data[i] += rhs.m_data[i];
        return *this;
    }

    template<typename Ty, std::size_t N, std::size_t M>
    TMatrix<Ty, N, M>& TMatrix<Ty, N, M>::operator-=(const TMatrix& rhs)
    {
        for (std::size_t i = 0; i < N*M; i++)
            m_data[i] += rhs.m_data[i];
        return *this;
    }
        
    template<typename Ty, std::size_t N, std::size_t M, std::size_t L>
    TMatrix<Ty, N, L> operator*(const TMatrix<Ty, N, M>& lhs, const TMatrix<Ty, M, L>& rhs)
    {
        TMatrix<Ty, N, L> res;
        for (std::size_t i = 0; i < N; i++) 
        {
            for (std::size_t k = 0; k < M; k++)
            {
                for (std::size_t j = 0; j < L; j++)
                    res(i, j) += lhs(i, k) * rhs(k, j); 
            }
        }
        return res;
    }

    
    template<typename Ty, std::size_t N, std::size_t M>
    TMatrix<Ty, M, N> TMatrix<Ty, N, M>::transpose() const
    {
        TMatrix<Ty, M, N> res;
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < M; j++)
                res(j, i) = (*this)(i, j);
        }
        return res;
    }
    
    template<typename Ty, std::size_t N, std::size_t M>
    template<typename U>
    typename TMatrix<Ty, N, M>::TEnableIfComplex_t<TMatrix<Ty, M, N>, U> TMatrix<Ty, N, M>::adjoint() const
    {
        TMatrix<Ty, M, N> res;
        for (std::size_t i = 0; i < N; i++)
        {
            for (std::size_t j = 0; j < M; j++)
                res(j, i) = std::conj((*this)(i, j));
        }
        return res;
    }

    template<typename Ty, std::size_t N>
    TMatrix<Ty, N, 1> LinearSolve(TMatrix<Ty, N, N> A, TMatrix<Ty, N, 1> b)
    {
        // Do LU decomposition and only keep the U matrix (here A is transformed to U).
        // The transformation b' = (L^-1)*b is done on the fly
        for (std::size_t k = 0; k < N - 1; k++)
        {
            // find pivot
            std::size_t pivot_row = k;
            for (std::size_t i = k + 1; i < N; i++)
            {
                if (std::abs(A(i, k)) > std::abs(A(pivot_row, k)))
                    pivot_row = i;
            }

            // swap k-th row with pivoted row
            if (k != pivot_row)
            {
                for (std::size_t i = k; i < N; i++)
                    std::swap(A(k, i), A(pivot_row, i));
            }

            // eliminate
            for (std::size_t i = k + 1; i < N; i++)
            {
                Ty lambda = A(i, k) / A(k, k);
                for (std::size_t j = k; j < N; j++)
                    A(i, j) -= lambda * A(k, j);
                b(i) -= lambda * b(k);
            }
        }

        // A is in upper triagonal form (lower triganonal is not set to zero expilicitly for efficiency)
        // b is the original b transformed by the L matrix from the L-U decomposition
        // Now one can solve for x
        TMatrix<Ty, N, 1> x = b;
        for (std::size_t l = 0; l < N; l++)
        {
            std::size_t i = N - l - 1;
            for (std::size_t j = i + 1; j < N; j++)
                x(i) -= A(i, j) * x(j);
            x(i) /= A(i, i);
        }

        return x;
    }
}

#endif
