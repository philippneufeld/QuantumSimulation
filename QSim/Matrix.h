// Philipp Neufeld, 2021

#ifndef QSim_Matrix_H_
#define QSim_Matrix_H_

#include <cstdint>
#include <type_traits>
#include <algorithm>
#include <complex>
#include <cassert>

namespace QSim
{
    // forward declare CRTP class
    template<typename MT>
    class TMatrix;
    template<typename MT>
    class TVector;

    namespace Internal
    {
        template<typename MT, typename = void>
        struct TIsMatrix : std::integral_constant<bool, false> {};
        template<typename MT>
        struct TIsMatrix<MT, std::enable_if_t<std::is_base_of<TMatrix<std::decay_t<decltype(~(std::declval<MT>()))>>, std::decay_t<MT>>::value>> 
            : std::integral_constant<bool, true> {};
        template<typename MT>
        constexpr bool TIsMatrix_v = TIsMatrix<MT>::value;

        template<typename MT, typename=std::enable_if_t<TIsMatrix_v<MT>>>
        struct TMatrixDecay
        {
            using type = std::decay_t<decltype(~(std::declval<MT>()))>;
        };
        template<typename MT>
        using TMatrixDecay_t = typename TMatrixDecay<MT>::type;

        template<typename MT1, typename MT2>
        struct TMatrixAdditionResultImpl;
        template<typename MT1, typename MT2>
        struct TMatrixAdditionResult 
            : TMatrixAdditionResultImpl<TMatrixDecay_t<MT1>, TMatrixDecay_t<MT2>> {};
        template<typename MT1, typename MT2>
        using TMatrixAdditionResult_t = typename TMatrixAdditionResult<MT1, MT2>::type;
        
        template<typename MT1, typename MT2>
        struct TMatrixSubtractionResultImpl;
        template<typename MT1, typename MT2>
        struct TMatrixSubtractionResult 
            : TMatrixSubtractionResultImpl<TMatrixDecay_t<MT1>, TMatrixDecay_t<MT2>> {};
        template<typename MT1, typename MT2>
        using TMatrixSubtractionResult_t = typename TMatrixSubtractionResult<MT1, MT2>::type;
        
        template<typename MT1, typename MT2>
        struct TMatrixMultiplicationResultImpl;
        template<typename MT1, typename MT2>
        struct TMatrixMultiplicationResult 
            : TMatrixMultiplicationResultImpl<TMatrixDecay_t<MT1>, TMatrixDecay_t<MT2>> {};
        template<typename MT1, typename MT2>
        using TMatrixMultiplicationResult_t = typename TMatrixMultiplicationResult<MT1, MT2>::type;

        template<typename MT>
        struct TMatrixTranspositionResultImpl;
        template<typename MT>
        struct TMatrixTranspositionResult 
            : TMatrixTranspositionResultImpl<TMatrixDecay_t<MT>> {};
        template<typename MT>
        using TMatrixTranspositionResult_t = typename TMatrixTranspositionResult<MT>::type;

        template<typename MT>
        struct TMatrixRowTypeImpl;
        template<typename MT>
        struct TMatrixRowType
            : TMatrixRowTypeImpl<TMatrixDecay_t<MT>> {};
        template<typename MT>
        using TMatrixRowType_t = typename TMatrixRowType<MT>::type;

        template<typename MT>
        struct TMatrixColTypeImpl;
        template<typename MT>
        struct TMatrixColType
            : TMatrixColTypeImpl<TMatrixDecay_t<MT>> {};
        template<typename MT>
        using TMatrixColType_t = typename TMatrixColType<MT>::type;
    }


    template<typename MT>
    class TMatrix
    {
    public:  
        MT& operator~() { return static_cast<MT&>(*this); }
        const MT& operator~() const { return static_cast<const MT&>(*this); }

        void SetZero();

        template<typename MT2> 
        MT& operator+=(const TMatrix<MT2>& rhs);
        template<typename MT2> 
        MT& operator-=(const TMatrix<MT2>& rhs);
        template<typename MT2> 
        MT& operator*=(const TMatrix<MT2>& rhs);
        template<typename Ty, typename=std::enable_if_t<!Internal::TIsMatrix_v<Ty>>> 
        MT& operator*=(Ty s);
    };

    template<typename VT>
    class TVector : public TMatrix<VT> {};
    template<typename VT>
    class TRowVector : public TVector<VT> {};
    template<typename VT>
    class TColVector : public TVector<VT> {};

    // 1x1 matrices are considered column vectors
    template<typename MT>
    class TSingleElementMatrix : public TColVector<MT> {};
    
    template<typename MT>
    void TMatrix<MT>::SetZero()
    {
        for (std::size_t i = 0; i < (~(*this)).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~(*this)).Cols(); j++)
                (~(*this))(i, j) = 0;
        }
    }

    template<typename MT1, typename MT2, typename MT3,
        typename=std::enable_if_t<Internal::TIsMatrix_v<MT1> && Internal::TIsMatrix_v<MT2> && Internal::TIsMatrix_v<MT3>>>
    void MatrixAdd(TMatrix<MT1>& c, const TMatrix<MT2>& a, const TMatrix<MT3>& b)
    {
        assert((~a).Rows() == (~b).Rows());
        assert((~a).Cols() == (~b).Cols());
        assert((~a).Rows() == (~c).Rows());
        assert((~a).Cols() == (~c).Cols());

        for (std::size_t i = 0; i < (~c).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~c).Cols(); j++)
                (~c)(i, j) = (~a)(i, j) + (~b)(i, j);            
        }
    }
    
    template<typename MT1, typename MT2, typename MT3>
    void MatrixSub(TMatrix<MT1>& c, const TMatrix<MT2>& a, const TMatrix<MT3>& b)
    {
        assert((~a).Rows() == (~b).Rows());
        assert((~a).Cols() == (~b).Cols());
        assert((~a).Rows() == (~c).Rows());
        assert((~a).Cols() == (~c).Cols());

        for (std::size_t i = 0; i < (~c).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~c).Cols(); j++)
                (~c)(i, j) = (~a)(i, j) - (~b)(i, j);            
        }
    }

    template<typename MT1, typename MT2, typename MT3>
    void MatrixMul(TMatrix<MT1>& c, const TMatrix<MT2>& a, const TMatrix<MT3>& b)
    {
        assert((~a).Cols() == (~b).Rows());
        assert((~c).Rows() == (~a).Rows());
        assert((~c).Cols() == (~b).Cols());   

        for (std::size_t i = 0; i < (~c).Rows(); i++)
        {
            for (std::size_t k = 0; k < (~a).Cols(); k++)
            {
                for (std::size_t j = 0; j < (~c).Cols(); j++)
                    (~c)(i, j) += (~a)(i, k) * (~b)(k, j);
            }                
        }
    }

    template<typename MT1, typename MT2, typename Ty>
    void MatrixScalarMul(TMatrix<MT1>& c, const TMatrix<MT2>& a, Ty s)
    {
        assert((~c).Rows() == (~a).Rows());
        assert((~c).Cols() == (~a).Cols());   

        for (std::size_t i = 0; i < (~c).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~c).Cols(); j++)
                (~c)(i, j) = (~a)(i, j) * s;                
        }
    }
    
    template<typename MT>
    Internal::TMatrixTranspositionResult_t<MT> Transpose(const TMatrix<MT>& mat)
    {
        Internal::TMatrixTranspositionResult_t<MT> res((~mat).Cols(), (~mat).Rows());
        for (std::size_t i = 0; i < (~res).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~res).Cols(); j++)
                (~res)(i, j) = (~mat)(j, i);
        }
        return res;
    }

    template<typename MT>
    Internal::TMatrixTranspositionResult_t<MT> Adjoint(const TMatrix<MT>& mat)
    {
        Internal::TMatrixTranspositionResult_t<MT> res((~mat).Cols(), (~mat).Rows());
        for (std::size_t i = 0; i < (~res).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~res).Cols(); j++)
                (~res)(i, j) = std::conj((~mat)(j, i));
        }
        return res;
    }

    template<typename MT>
    Internal::TMatrixRowType_t<MT> GetRow(const TMatrix<MT>& mat, std::size_t i)
    {
        Internal::TMatrixRowType_t<MT> row((~mat).Cols());
        for (std::size_t j = 0; j < (~mat).Cols(); j++)
            row(j) = (~mat)(i, j);
        return row;
    }

    template<typename MT>
    Internal::TMatrixColType_t<MT> GetCol(const TMatrix<MT>& mat, std::size_t i)
    {
        Internal::TMatrixColType_t<MT> col((~mat).Rows());
        for (std::size_t j = 0; j < (~mat).Rows(); j++)
            col(j) = (~mat)(j, i);
        return col;
    }

    template<typename MT, typename VT>
    void SetRow(TMatrix<MT>& mat, const TVector<VT>& v, std::size_t i)
    {
        assert((~mat).Cols() == (~v).Size());
        for (std::size_t j = 0; j < (~mat).Cols(); j++)
            (~mat)(i, j) = (~v)(j);
    }

    template<typename MT, typename VT>
    void SetCol(TMatrix<MT>& mat, const TVector<VT>& v, std::size_t i)
    {
        assert((~mat).Rows() == (~v).Size());
        for (std::size_t j = 0; j < (~mat).Rows(); j++)
            (~mat)(j, i) = (~v)(j);
    }

    template<typename MT1, typename MT2, typename=void>
    Internal::TMatrixAdditionResult_t<MT1, MT2> operator+(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    { 
        Internal::TMatrixAdditionResult_t<MT1, MT2> res((~rhs).Rows(), (~rhs).Cols());
        MatrixAdd(res, lhs, rhs);   
        return res;
    }

    template<typename MT1, typename MT2, typename=void>
    Internal::TMatrixSubtractionResult_t<MT1, MT2> operator-(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    {
        Internal::TMatrixAdditionResult_t<MT1, MT2> res((~rhs).Rows(), (~rhs).Cols());
        MatrixSub(res, lhs, rhs);   
        return res;
    }

    template<typename MT1, typename MT2, typename=void>
    Internal::TMatrixMultiplicationResult_t<MT1, MT2> operator*(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    { 
        Internal::TMatrixMultiplicationResult_t<MT1, MT2> res((~lhs).Rows(), (~rhs).Cols());
        MatrixMul(res, lhs, rhs);
        return res;
    }

    template<typename MT1, typename Ty, typename=std::enable_if_t<!Internal::TIsMatrix_v<Ty>>>
    MT1 operator*(const TMatrix<MT1>& lhs, Ty s) 
    { 
        MT1 res((~lhs).Rows(), (~lhs).Cols());
        MatrixScalarMul(res, lhs, s);
        return res;
    }

    template<typename MT1, typename Ty, typename=std::enable_if_t<!Internal::TIsMatrix_v<Ty>>>
    MT1 operator*(Ty s, const TMatrix<MT1>& rhs) 
    { 
        return rhs * s;
    }

    template<typename MT>
    template<typename MT2> 
    MT& TMatrix<MT>::operator+=(const TMatrix<MT2>& rhs) 
    {
        MatrixAdd(*this, *this, rhs);
        return ~(*this);
    }

    template<typename MT>
    template<typename MT2> 
    MT& TMatrix<MT>::operator-=(const TMatrix<MT2>& rhs) 
    {
        MatrixSub(*this, *this, rhs);
        return ~(*this);
    }

    template<typename MT>
    template<typename MT2> 
    MT& TMatrix<MT>::operator*=(const TMatrix<MT2>& rhs) 
    {
        MatrixMul(*this, *this, rhs);
        return ~(*this);
    }

    template<typename MT>
    template<typename Ty, typename> 
    MT& TMatrix<MT>::operator*=(Ty s) 
    {
        MatrixScalarMul(*this, *this, s);
        return ~(*this);
    }

    namespace Internal
    {
        template<typename MT, typename IT>
        class TMatrixRCIterator_base
        {
        public:
            TMatrixRCIterator_base(const TMatrix<MT>& mat, std::size_t idx) 
                : m_mat(mat), m_idx(idx) { }
            TMatrixRCIterator_base(const TMatrixRCIterator_base&) = default;
            TMatrixRCIterator_base& operator=(TMatrixRCIterator_base&) = default;

            IT& operator++() { m_idx++; return *this; }
            IT& operator--() { m_idx--; return *this; }
            IT operator++(int) { return IT(m_mat, m_idx++); }
            IT operator--(int) { return IT(m_mat, m_idx--); }

            bool operator==(const IT& rhs) { return (m_idx == rhs.m_idx) && (&m_mat) == (&rhs.m_mat); }
            bool operator!=(const IT& rhs) { return !((*this) == rhs); }
            bool operator<(const IT& rhs) { return ((&m_mat) == (&rhs.m_mat) && (m_idx < rhs.m_idx)) || (&m_mat) < (&rhs.m_mat); }
            bool operator<=(const IT& rhs) { return ((*this) == rhs|| (*this) < rhs); }
            bool operator>=(const IT& rhs) { return !((*this) < rhs); }
            bool operator>(const IT& rhs) { return !((*this) <= rhs); }

            std::ptrdiff_t operator-(const IT& rhs) 
            { 
                if ((&m_mat) != (&rhs.m_mat)) return -1;
                return static_cast<std::ptrdiff_t>(m_idx) - static_cast<std::ptrdiff_t>(rhs.m_idx); 
            }
            IT operator+(std::ptrdiff_t off) const { return IT(this->m_mat, this->m_idx + off); }
            IT operator-(std::ptrdiff_t off) const { return IT(this->m_mat, this->m_idx - off); }

        protected:
            const TMatrix<MT>& m_mat;
            std::size_t m_idx;
        };
    }

    template<typename MT>
    class TMatrixRowIterator 
        : public Internal::TMatrixRCIterator_base<MT, TMatrixRowIterator<MT>>
    {
        using MyT = TMatrixRowIterator<MT>;
    public:
        TMatrixRowIterator(const TMatrix<MT>& mat, std::size_t idx) 
            : Internal::TMatrixRCIterator_base<MT, MyT>(mat, idx) { }
        Internal::TMatrixRowType_t<MT> operator*() const { return GetRow(this->m_mat, this->m_idx); }
    };

    template<typename MT>
    class TMatrixColIterator
        : public Internal::TMatrixRCIterator_base<MT, TMatrixColIterator<MT>>
    {
        using MyT = TMatrixColIterator<MT>;
    public:
        TMatrixColIterator(const TMatrix<MT>& mat, std::size_t idx) 
            : Internal::TMatrixRCIterator_base<MT, MyT>(mat, idx) { }
        Internal::TMatrixColType_t<MT> operator*() const { return GetCol(this->m_mat, this->m_idx); }
    };

    template<typename MT>
    TMatrixRowIterator<MT> GetRowIteratorBegin(const TMatrix<MT>& mat)
    {
        return TMatrixRowIterator<MT>(mat, 0);
    }

    template<typename MT>
    TMatrixRowIterator<MT> GetRowIteratorEnd(const TMatrix<MT>& mat)
    {
        return TMatrixRowIterator<MT>(mat, (~mat).Rows());
    }

    template<typename MT>
    TMatrixColIterator<MT> GetColIteratorBegin(const TMatrix<MT>& mat)
    {
        return TMatrixColIterator<MT>(mat, 0);
    }

    template<typename MT>
    TMatrixColIterator<MT> GetColIteratorEnd(const TMatrix<MT>& mat)
    {
        return TMatrixColIterator<MT>(mat, (~mat).Cols());
    }

    //
    // Statically sized matrix
    //

    namespace Internal
    {
        template<typename Ty, std::size_t N, std::size_t M, typename MyT>
        struct TStaticMatrixHelper
        {
            using type = TMatrix<MyT>;
        };

        template<typename Ty, std::size_t N, typename MyT>
        struct TStaticMatrixHelper<Ty, N, 1, MyT>
        {
            using type = TColVector<MyT>;
        };

        template<typename Ty, std::size_t M, typename MyT>
        struct TStaticMatrixHelper<Ty, 1, M, MyT>
        {
            using type = TRowVector<MyT>;
        };

        template<typename Ty, typename MyT>
        struct TStaticMatrixHelper<Ty, 1, 1, MyT>
        {
            using type = TSingleElementMatrix<MyT>;
        };


        template<typename Ty, std::size_t N, std::size_t M, typename MyT>
        using TStaticMatrixHelper_t = typename TStaticMatrixHelper<Ty, N, M, MyT>::type;
    }

    template<typename Ty, std::size_t N, std::size_t M>
    class TStaticMatrix : public Internal::TStaticMatrixHelper_t<Ty, N, M, TStaticMatrix<Ty, N, M>>
    {
    public:
        TStaticMatrix() : m_data{} {} // m_data is initialized to its default value this way
        TStaticMatrix(const Ty(&data)[N*M]);
        template<std::size_t Dummy = 1, typename = std::enable_if_t<Dummy == 1 && (N == 1 || M == 1)>>
        TStaticMatrix(std::size_t) : TStaticMatrix() {}
        TStaticMatrix(std::size_t, std::size_t) : TStaticMatrix() {}
        TStaticMatrix(std::size_t, std::size_t, const Ty(&data)[N*M]) : TStaticMatrix(data) {}
        ~TStaticMatrix() = default;

        template<typename MT>
        TStaticMatrix(const TMatrix<MT>& rhs);
        template<typename MT>
        TStaticMatrix& operator=(const TMatrix<MT>& rhs);
       
        Ty& operator()(std::size_t i, std::size_t j) { return m_data[i*M+j]; }
        const Ty& operator()(std::size_t i, std::size_t j) const { return m_data[i*M+j]; }

        template<std::size_t K=N, std::size_t L=M, typename=std::enable_if_t<K==1 || L==1>>
        Ty& operator()(std::size_t i) { return m_data[i]; }
        template<std::size_t K=N, std::size_t L=M, typename=std::enable_if_t<K==1 || L==1>>
        const Ty& operator()(std::size_t i) const { return m_data[i]; }

        Ty& operator[](std::size_t i) { return m_data[i]; }
        const Ty& operator[](std::size_t i) const { return m_data[i]; }

        constexpr std::size_t Rows() const { return N; }
        constexpr std::size_t Cols() const { return M; }
        constexpr std::size_t Size() const { return Rows()*Cols(); }

    private:
        Ty m_data[N*M];
    };

    template<typename Ty, std::size_t N>
    using TStaticRowVector = TStaticMatrix<Ty, 1, N>;
    template<typename Ty, std::size_t N>
    using TStaticColVector = TStaticMatrix<Ty, N, 1>;

    template<typename Ty, std::size_t N, std::size_t M>
    TStaticMatrix<Ty, N, M>::TStaticMatrix(const Ty(&data)[N*M])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = data[i*M + j];
        }
    }

    template<typename Ty, std::size_t N, std::size_t M>
    template<typename MT>
    TStaticMatrix<Ty, N, M>::TStaticMatrix(const TMatrix<MT>& rhs)
    {
        assert((~rhs).Rows() == N);
        assert((~rhs).Cols() == M);

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);
        }
    }
    
    template<typename Ty, std::size_t N, std::size_t M>
    template<typename MT>
    TStaticMatrix<Ty, N, M>& TStaticMatrix<Ty, N, M>::operator=(const TMatrix<MT>& rhs)
    {
        assert((~rhs).Rows() == N);
        assert((~rhs).Cols() == M);

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);
        }

        return *this;
    }   

    namespace Internal
    {
        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<TStaticMatrix<Ty, N, M>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<TStaticMatrix<Ty, N, M>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M, std::size_t L>
        struct TMatrixMultiplicationResultImpl<TStaticMatrix<Ty, N, L>, TStaticMatrix<Ty, L, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixTranspositionResultImpl<TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, M, N>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixRowTypeImpl<TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticRowVector<Ty, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixColTypeImpl<TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticColVector<Ty, N>;
        };
    }

    //
    // Hybridly sized matrix
    //

    template<typename Ty, std::size_t N, bool colDyn>
    class THybridMatrix : public std::conditional_t<N==1, 
        std::conditional_t<colDyn, TRowVector<THybridMatrix<Ty, N, colDyn>>, TColVector<THybridMatrix<Ty, N, colDyn>>>, TMatrix<THybridMatrix<Ty, N, colDyn>>>
    {
    public:
        THybridMatrix() : m_dynDim(0), m_data(nullptr) {}
        template<std::size_t Dummy = 1, typename = std::enable_if_t<Dummy == 1 && (N == 1)>>
        THybridMatrix(std::size_t size) : THybridMatrix(colDyn ? N : size, colDyn ? size : N) {}
        THybridMatrix(std::size_t rows, std::size_t cols);
        THybridMatrix(std::size_t rows, std::size_t cols, const Ty* data);
        THybridMatrix(std::size_t rows, std::size_t cols, const std::initializer_list<double>& lst)
            : THybridMatrix(rows, cols, lst.begin()) {}
        ~THybridMatrix() { if (m_data) delete[] m_data; m_data = nullptr; }

        template<typename MT>
        THybridMatrix(const TMatrix<MT>& rhs);
        template<typename MT>
        THybridMatrix& operator=(const TMatrix<MT>& rhs);

        THybridMatrix(const THybridMatrix& rhs);
        THybridMatrix& operator=(const THybridMatrix& rhs);

        THybridMatrix(THybridMatrix&& rhs) 
            : m_dynDim(rhs.m_dynDim), m_data(rhs.m_data) { rhs.m_data = nullptr; }
        THybridMatrix& operator=(THybridMatrix&& rhs);
       
        Ty& operator()(std::size_t i, std::size_t j) { return m_data[i*Cols()+j]; }
        const Ty& operator()(std::size_t i, std::size_t j) const { return m_data[i*Cols()+j]; }

        template<std::size_t K=N, typename=std::enable_if_t<K==1>>
        Ty& operator()(std::size_t i) { return m_data[i]; }
        template<std::size_t K=N, typename=std::enable_if_t<K==1>>
        const Ty& operator()(std::size_t i) const { return m_data[i]; }

        Ty& operator[](std::size_t i) { return m_data[i]; }
        const Ty& operator[](std::size_t i) const { return m_data[i]; }

        std::size_t Rows() const { return colDyn ? N : m_dynDim; }
        std::size_t Cols() const { return colDyn ? m_dynDim : N; }
        constexpr std::size_t Size() const { return Rows()*Cols(); }

        void Resize(std::size_t dynDim);

    private:
        std::size_t m_dynDim;
        Ty* m_data;
    };

    template<typename Ty>
    using TDynamicRowVector = THybridMatrix<Ty, 1, true>;
    template<typename Ty>
    using TDynamicColVector = THybridMatrix<Ty, 1, false>;


    template<typename Ty, std::size_t N, bool colDyn>
    THybridMatrix<Ty, N, colDyn>::THybridMatrix(std::size_t rows, std::size_t cols)
        : m_dynDim(colDyn ? cols : rows), m_data(new Ty[N*m_dynDim])
    {
        assert((colDyn ? rows : cols) == N);
        this->SetZero();
    }

    template<typename Ty, std::size_t N, bool colDyn>
    THybridMatrix<Ty, N, colDyn>::THybridMatrix(std::size_t rows, std::size_t cols, const Ty* data)
        : m_dynDim(colDyn ? cols : rows), m_data(new Ty[N*m_dynDim])
    {
        assert((colDyn ? rows : cols) == N);
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = data[i*Cols()+j];
        } 
    }

    template<typename Ty, std::size_t N, bool colDyn>
    template<typename MT>
    THybridMatrix<Ty, N, colDyn>::THybridMatrix(const TMatrix<MT>& rhs)
        : m_dynDim(colDyn ? (~rhs).Cols() : (~rhs).Rows()), m_data(new Ty[N*m_dynDim])
    {
        assert((colDyn ? (~rhs).Rows() : (~rhs).Cols()) == N);
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }    
    }

    template<typename Ty, std::size_t N, bool colDyn>
    template<typename MT>
    THybridMatrix<Ty, N, colDyn>& THybridMatrix<Ty, N, colDyn>::operator=(const TMatrix<MT>& rhs)
    {
        assert((colDyn ? (~rhs).Rows() : (~rhs).Cols()) == N);
        this->Resize(colDyn ? (~rhs).Cols() : (~rhs).Rows());

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }

        return *this;
    }

    template<typename Ty, std::size_t N, bool colDyn>
    THybridMatrix<Ty, N, colDyn>::THybridMatrix(const THybridMatrix<Ty, N, colDyn>& rhs)
        : m_dynDim(rhs.m_dynDim), m_data(new Ty[N*m_dynDim])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }
    }
    
    template<typename Ty, std::size_t N, bool colDyn>
    THybridMatrix<Ty, N, colDyn>& THybridMatrix<Ty, N, colDyn>::operator=(const THybridMatrix<Ty, N, colDyn>& rhs)
    {
        this->Resize(rhs.m_dynDim);

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }

        return *this;
    }

    template<typename Ty, std::size_t N, bool colDyn>
    THybridMatrix<Ty, N, colDyn>& THybridMatrix<Ty, N, colDyn>::operator=(THybridMatrix<Ty, N, colDyn>&& rhs)
    {
        std::swap(m_dynDim, rhs.m_dynDim);
        std::swap(m_data, rhs.m_data);
        return *this;
    }

    template<typename Ty, std::size_t N, bool colDyn>
    void THybridMatrix<Ty, N, colDyn>::Resize(std::size_t dynDim)
    {
        if (m_dynDim != dynDim)
        {
            this->~THybridMatrix();
            if (dynDim > 0)
                m_data = new Ty[N*dynDim];
        }
        m_dynDim = dynDim;
    }

    namespace Internal
    {
        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, N, true>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, M, false>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<THybridMatrix<Ty, N, true>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<THybridMatrix<Ty, M, false>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, bool colDyn>
        struct TMatrixAdditionResultImpl<THybridMatrix<Ty, N, colDyn>, THybridMatrix<Ty, N, colDyn>>
        {
            using type = THybridMatrix<Ty, N, colDyn>;
        };

        template<typename Ty, std::size_t N, std::size_t M, bool colDyn>
        struct TMatrixAdditionResultImpl<THybridMatrix<Ty, N, colDyn>, THybridMatrix<Ty, M, !colDyn>>
        {
            using type = TStaticMatrix<Ty, colDyn ? N : M, colDyn ? M : N>;
        };


        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, N, true>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, M, false>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<THybridMatrix<Ty, N, true>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<THybridMatrix<Ty, M, false>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, bool colDyn>
        struct TMatrixSubtractionResultImpl<THybridMatrix<Ty, N, colDyn>, THybridMatrix<Ty, N, colDyn>>
        {
            using type = THybridMatrix<Ty, N, colDyn>;
        };

        template<typename Ty, std::size_t N, std::size_t M, bool colDyn>
        struct TMatrixSubtractionResultImpl<THybridMatrix<Ty, N, colDyn>, THybridMatrix<Ty, M, !colDyn>>
        {
            using type = TStaticMatrix<Ty, colDyn ? N : M, colDyn ? M : N>;
        };



        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, M, true>>
        {
            using type = THybridMatrix<Ty, N, true>;
        };

        template<typename Ty, std::size_t N, std::size_t M, std::size_t K>
        struct TMatrixMultiplicationResultImpl<TStaticMatrix<Ty, N, M>, THybridMatrix<Ty, K, false>>
        {
            using type = TStaticMatrix<Ty, N, K>;
        };

        template<typename Ty, std::size_t N, std::size_t M, std::size_t K>
        struct TMatrixMultiplicationResultImpl<THybridMatrix<Ty, K, true>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, K, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResultImpl<THybridMatrix<Ty, N, false>, TStaticMatrix<Ty, N, M>>
        {
            using type = THybridMatrix<Ty, M, false>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResultImpl<THybridMatrix<Ty, N, true>, THybridMatrix<Ty, M, false>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M, bool colDyn>
        struct TMatrixMultiplicationResultImpl<THybridMatrix<Ty, N, colDyn>, THybridMatrix<Ty, M, colDyn>>
        {
            using type = THybridMatrix<Ty, colDyn ? N : M, colDyn>;
        };

        template<typename Ty, std::size_t N, bool colDyn>
        struct TMatrixTranspositionResultImpl<THybridMatrix<Ty, N, colDyn>>
        {
            using type = THybridMatrix<Ty, N, !colDyn>;
        };  

        template<typename Ty, std::size_t N, bool colDyn>
        struct TMatrixRowTypeImpl<THybridMatrix<Ty, N, colDyn>>
        {
            using type = std::conditional_t<colDyn, TStaticRowVector<Ty, N>, TDynamicRowVector<Ty>>;
        };

        template<typename Ty, std::size_t N, bool colDyn>
        struct TMatrixColTypeImpl<THybridMatrix<Ty, N, colDyn>>
        {
            using type = std::conditional_t<!colDyn, TStaticColVector<Ty, N>, TDynamicColVector<Ty>>;
        };
    }

    //
    // Dynamically sized matrix
    //

    template<typename Ty>
    class TDynamicMatrix : public TMatrix<TDynamicMatrix<Ty>>
    {
    public:
        TDynamicMatrix() : m_rows(0), m_cols(0), m_data(nullptr) {}
        TDynamicMatrix(std::size_t rows, std::size_t cols);
        TDynamicMatrix(std::size_t rows, std::size_t cols, const Ty* data);
        TDynamicMatrix(std::size_t rows, std::size_t cols, const std::initializer_list<double>& lst)
            : TDynamicMatrix(rows, cols, lst.begin()) {}
        ~TDynamicMatrix() { if (m_data) delete[] m_data; m_data = nullptr; }

        template<typename MT>
        TDynamicMatrix(const TMatrix<MT>& rhs);
        template<typename MT>
        TDynamicMatrix& operator=(const TMatrix<MT>& rhs);

        TDynamicMatrix(const TDynamicMatrix& rhs);
        TDynamicMatrix& operator=(const TDynamicMatrix& rhs);

        TDynamicMatrix(TDynamicMatrix&& rhs) 
            : m_rows(rhs.m_rows), m_cols(rhs.m_cols), m_data(rhs.m_data) { rhs.m_data = nullptr; }
        TDynamicMatrix& operator=(TDynamicMatrix&& rhs);
       
        Ty& operator()(std::size_t i, std::size_t j) { return m_data[i*Cols()+j]; }
        const Ty& operator()(std::size_t i, std::size_t j) const { return m_data[i*Cols()+j]; }

        Ty& operator[](std::size_t i) { return m_data[i]; }
        const Ty& operator[](std::size_t i) const { return m_data[i]; }

        std::size_t Rows() const { return m_rows; }
        std::size_t Cols() const { return m_cols; }
        constexpr std::size_t Size() const { return Rows()*Cols(); }

        void Resize(std::size_t rows, std::size_t cols);

    private:
        std::size_t m_rows;
        std::size_t m_cols;
        Ty* m_data;
    };


    template<typename Ty>
    TDynamicMatrix<Ty>::TDynamicMatrix(std::size_t rows, std::size_t cols)
        : m_rows(rows), m_cols(cols), m_data(new Ty[m_rows*m_cols])
    {
        this->SetZero();
    }

    template<typename Ty>
    TDynamicMatrix<Ty>::TDynamicMatrix(std::size_t rows, std::size_t cols, const Ty* data)
        : m_rows(rows), m_cols(cols), m_data(new Ty[m_rows*m_cols])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = data[i*m_cols+j];
        } 
    }

    template<typename Ty>
    template<typename MT>
    TDynamicMatrix<Ty>::TDynamicMatrix(const TMatrix<MT>& rhs)
        : m_rows((~rhs).Rows()), m_cols((~rhs).Cols()), m_data(new Ty[m_rows*m_cols])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }    
    }

    template<typename Ty>
    template<typename MT>
    TDynamicMatrix<Ty>& TDynamicMatrix<Ty>::operator=(const TMatrix<MT>& rhs)
    {
        Resize((~rhs).Rows(), (~rhs).Cols());

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }

        return *this;
    }

    template<typename Ty>
    TDynamicMatrix<Ty>::TDynamicMatrix(const TDynamicMatrix<Ty>& rhs)
        : m_rows(rhs.m_rows), m_cols(rhs.m_cols), m_data(new Ty[m_rows*m_cols])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }
    }
    
    template<typename Ty>
    TDynamicMatrix<Ty>& TDynamicMatrix<Ty>::operator=(const TDynamicMatrix<Ty>& rhs)
    {
        this->Resize(rhs.Rows(), rhs.Cols());

        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = (~rhs)(i, j);            
        }

        return *this;
    }

    template<typename Ty>
    TDynamicMatrix<Ty>& TDynamicMatrix<Ty>::operator=(TDynamicMatrix<Ty>&& rhs)
    {
        std::swap(m_rows, rhs.m_rows);
        std::swap(m_cols, rhs.m_cols);
        std::swap(m_data, rhs.m_data);
        return *this;
    }

    template<typename Ty>
    void TDynamicMatrix<Ty>::Resize(std::size_t rows, std::size_t cols)
    {
        if (m_rows*m_cols != rows*cols)
        {
            this->~TDynamicMatrix();
            if (rows*cols > 0)
                m_data = new Ty[rows*cols];
        }
        m_rows = rows;
        m_cols = cols;
    }

    namespace Internal
    {
        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResultImpl<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty>
        struct TMatrixAdditionResultImpl<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };


        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResultImpl<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty>
        struct TMatrixSubtractionResultImpl<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };


        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResultImpl<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = THybridMatrix<Ty, N, true>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResultImpl<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = THybridMatrix<Ty, M, false>;
        };

        template<typename Ty>
        struct TMatrixMultiplicationResultImpl<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty, std::size_t N>
        struct TMatrixMultiplicationResultImpl<THybridMatrix<Ty, N, false>, THybridMatrix<Ty, N, true>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty>
        struct TMatrixTranspositionResultImpl<TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        }; 

        template<typename Ty>
        struct TMatrixRowTypeImpl<TDynamicMatrix<Ty>>
        {
            using type = TDynamicRowVector<Ty>;
        };

        template<typename Ty>
        struct TMatrixColTypeImpl<TDynamicMatrix<Ty>>
        {
            using type = TDynamicColVector<Ty>;
        };
    }


    //
    // Algorithms
    //

    template<typename Ty = double, bool colVec = true>
    auto CreateZeros(std::size_t steps)
    {
        std::conditional_t<colVec, TDynamicColVector<Ty>, TDynamicRowVector<Ty>> vec;
        vec.Resize(steps);
        vec.SetZero();
        return vec;
    }

    template<typename MT>
    auto CreateZerosLike(const TMatrix<MT>& like)
    {
        Internal::TMatrixDecay_t<MT> vec((~like).Rows(), (~like).Cols());
        vec.SetZero();
        return vec;
    }

    template<typename Ty, bool colVec = true>
    auto CreateLinspace(Ty start, Ty stop, std::size_t steps)
    {
        std::conditional_t<colVec, TDynamicColVector<Ty>, TDynamicRowVector<Ty>> vec;
        vec.Resize(steps);

        if (steps == 1)
            vec(0) = (start + stop) / 2;
        else if(steps > 1)
        {
            Ty tyStep = (stop - start) / (steps - 1);
            for (std::size_t i = 0; i < steps; i++)
                vec(i) = start + i * tyStep;
        }

        return vec;
    }

    template<typename Ty>
    auto CreateLinspaceRow(Ty start, Ty stop, std::size_t steps)
    {
        return CreateLinspace<Ty, false>(start, stop, steps);
    }

    template<typename Ty>
    auto CreateLinspaceCol(Ty start, Ty stop, std::size_t steps)
    {
        return CreateLinspace<Ty, true>(start, stop, steps);
    }

    template<typename MT, typename VT>
    Internal::TMatrixDecay_t<VT> LinearSolve(const TMatrix<MT>& A, const TVector<VT>& b)
    {
        assert((~A).Rows() == (~A).Cols());
        assert((~A).Cols() == (~b).Rows());
        assert((~b).Cols() == 1);
        
        Internal::TMatrixDecay_t<MT> U(A);
        Internal::TMatrixDecay_t<VT> y(b);

        // Do LU decomposition and only keep the U matrix (here A is transformed to U).
        // The transformation y = (L^-1)*b is done on the fly
        for (std::size_t k = 0; k < (~U).Rows() - 1; k++)
        {
            // find pivot
            std::size_t pivot_row = k;
            for (std::size_t i = k + 1; i < (~U).Rows(); i++)
            {
                if (std::abs((~U)(i, k)) > std::abs((~U)(pivot_row, k)))
                    pivot_row = i;
            }

            // swap k-th row with pivoted row
            if (k != pivot_row)
            {
                for (std::size_t i = k; i < (~U).Rows(); i++)
                    std::swap((~U)(k, i), (~U)(pivot_row, i));
                std::swap((~y)(k), (~y)(pivot_row));
            }

            // eliminate
            auto factor = 1.0 / (~U)(k, k);
            for (std::size_t i = k; i < (~U).Rows(); i++)
                (~U)(k, i) *= factor;
            (~U)(k, k) = 1;
            (~y)(k) *= factor;
            

            for (std::size_t i = k + 1; i < (~U).Rows(); i++)
            {
                auto lambda = (~U)(i, k); // / (~U)(k, k);
                for (std::size_t j = k; j < (~U).Rows(); j++)
                    (~U)(i, j) -= lambda * (~U)(k, j);
                (~y)(i) -= lambda * (~y)(k);
            }
        }

        // U is in upper triagonal form (lower triganonal is not set to zero expilicitly for efficiency)
        // y is the original y transformed by the L matrix from the L-U decomposition
        // Now one can solve for x
        Internal::TMatrixDecay_t<VT> x = (~y);
        for (std::size_t l = 0; l < (~U).Rows(); l++)
        {
            std::size_t i = (~U).Rows() - l - 1;
            for (std::size_t j = i + 1; j < (~U).Rows(); j++)
                (~x)(i) -= (~U)(i, j) * (~x)(j);
            (~x)(i) /= (~U)(i, i);
        }

        return (~x);
    }

}

#endif
