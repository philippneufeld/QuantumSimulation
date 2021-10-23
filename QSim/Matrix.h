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
        template<typename MT1, typename MT2>
        struct TMatrixAdditionResult;
        template<typename MT1, typename MT2>
        using TMatrixAdditionResult_t = typename TMatrixAdditionResult<MT1, MT2>::type;
        
        template<typename MT1, typename MT2>
        struct TMatrixSubtractionResult;
        template<typename MT1, typename MT2>
        using TMatrixSubtractionResult_t = typename TMatrixSubtractionResult<MT1, MT2>::type;
        
        template<typename MT1, typename MT2>
        struct TMatrixMultiplicationResult;
        template<typename MT1, typename MT2>
        using TMatrixMultiplicationResult_t = typename TMatrixMultiplicationResult<MT1, MT2>::type;

        template<typename MT>
        struct TMatrixTranspositionResult;
        template<typename MT>
        using TMatrixTranspositionResult_t = typename TMatrixTranspositionResult<MT>::type;
    }


    template<typename MT>
    class TMatrix
    {
    public:  
        MT& operator~() { return static_cast<MT&>(*this); }
        const MT& operator~() const { return static_cast<const MT&>(*this); }
    };

    template<typename MT1, typename MT2>
    Internal::TMatrixAdditionResult_t<MT1, MT2> operator+(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    { 
        assert((~lhs).Rows() == (~rhs).Rows());
        assert((~lhs).Cols() == (~rhs).Cols());

        Internal::TMatrixAdditionResult_t<MT1, MT2> res((~rhs).Rows(), (~rhs).Cols());
        for (std::size_t i = 0; i < (~res).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~res).Cols(); j++)
                (~res)(i, j) = (~lhs)(i, j) + (~rhs)(i, j);            
        }
        
        return res;
    }

    template<typename MT1, typename MT2>
    Internal::TMatrixSubtractionResult_t<MT1, MT2> operator-(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    { 
        assert((~lhs).Rows() == (~rhs).Rows());
        assert((~lhs).Cols() == (~rhs).Cols());

        Internal::TMatrixSubtractionResult_t<MT1, MT2> res((~rhs).Rows(), (~rhs).Cols());
        for (std::size_t i = 0; i < (~res).Rows(); i++)
        {
            for (std::size_t j = 0; j < (~res).Cols(); j++)
                (~res)(i, j) = (~lhs)(i, j) - (~rhs)(i, j);            
        }

        return res;
    }
    
    template<typename MT1, typename MT2>
    Internal::TMatrixMultiplicationResult_t<MT1, MT2> operator*(const TMatrix<MT1>& lhs, const TMatrix<MT2>& rhs) 
    { 
        assert((~rhs).Cols() == (~rhs).Rows());

        Internal::TMatrixMultiplicationResult_t<MT1, MT2> res((~lhs).Rows(), (~rhs).Cols());
        for (std::size_t i = 0; i < (~res).Rows(); i++)
        {
            for (std::size_t k = 0; k < (~lhs).Cols(); k++)
            {
                for (std::size_t j = 0; j < (~res).Cols(); j++)
                    (~res)(i, j) += (~lhs)(i, k) * (~rhs)(k, j);
            }                
        }

        return res;
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

    
    //
    // Statically sized matrix
    //

    template<typename Ty, std::size_t N, std::size_t M>
    class TStaticMatrix : public TMatrix<TStaticMatrix<Ty, N, M>>
    {
    public:
        TStaticMatrix() : m_data{} {} // m_data is initialized to its default value this way
        TStaticMatrix(const Ty(&data)[N*M]);
        TStaticMatrix(std::size_t, std::size_t) : TStaticMatrix() {}
        TStaticMatrix(std::size_t, std::size_t, const Ty(&data)[N*M]) : TStaticMatrix(data) {}
        ~TStaticMatrix() = default;

        template<typename MT>
        TStaticMatrix(const TMatrix<MT>& rhs);
        template<typename MT>
        TStaticMatrix& operator=(const TMatrix<MT>& rhs);
       
        Ty& operator()(std::size_t i, std::size_t j) { return m_data[i*M+j]; }
        const Ty& operator()(std::size_t i, std::size_t j) const { return m_data[i*M+j]; }

        constexpr std::size_t Rows() const { return N; }
        constexpr std::size_t Cols() const { return M; }

    private:
        Ty m_data[N*M];
    };

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
        struct TMatrixAdditionResult<TStaticMatrix<Ty, N, M>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResult<TStaticMatrix<Ty, N, M>, TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M, std::size_t L>
        struct TMatrixMultiplicationResult<TStaticMatrix<Ty, N, L>, TStaticMatrix<Ty, L, M>>
        {
            using type = TStaticMatrix<Ty, N, M>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixTranspositionResult<TStaticMatrix<Ty, N, M>>
        {
            using type = TStaticMatrix<Ty, M, N>;
        };
    }

    //
    // Dynamically sized matrix
    //

    template<typename Ty>
    class TDynamicMatrix : public TMatrix<TDynamicMatrix<Ty>>
    {
    public:
        TDynamicMatrix(std::size_t rows, std::size_t cols);
        TDynamicMatrix(std::size_t rows, std::size_t cols, const Ty* data);
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

        std::size_t Rows() const { return m_rows; }
        std::size_t Cols() const { return m_cols; }

    private:
        std::size_t m_rows;
        std::size_t m_cols;
        Ty* m_data;
    };


    template<typename Ty>
    TDynamicMatrix<Ty>::TDynamicMatrix(std::size_t rows, std::size_t cols)
        : m_rows(rows), m_cols(cols), m_data(new Ty[m_rows*m_cols])
    {
        for (std::size_t i = 0; i < Rows(); i++)
        {
            for (std::size_t j = 0; j < Cols(); j++)
                (*this)(i, j) = Ty();
        } 
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
        if (m_rows*m_cols != (~rhs).Rows()*(~rhs).Cols())
        {
            this->~TDynamicMatrix();
            m_data = new Ty[(~rhs).Rows()*(~rhs).Cols()];
        }
        m_rows = (~rhs).Rows();
        m_cols = (~rhs).Cols(); 

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
        if (m_rows*m_cols != (~rhs).Rows()*(~rhs).Cols())
        {
            this->~TDynamicMatrix();
            m_data = new Ty[(~rhs).Rows()*(~rhs).Cols()];
        }
        m_rows = (~rhs).Rows();
        m_cols = (~rhs).Cols();

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

    namespace Internal
    {
        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResult<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixAdditionResult<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty>
        struct TMatrixAdditionResult<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };


        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResult<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixSubtractionResult<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty>
        struct TMatrixSubtractionResult<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };


        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResult<TStaticMatrix<Ty, N, M>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty, std::size_t N, std::size_t M>
        struct TMatrixMultiplicationResult<TDynamicMatrix<Ty>, TStaticMatrix<Ty, N, M>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty>
        struct TMatrixMultiplicationResult<TDynamicMatrix<Ty>, TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        };

        template<typename Ty>
        struct TMatrixTranspositionResult<TDynamicMatrix<Ty>>
        {
            using type = TDynamicMatrix<Ty>;
        }; 
    }


    template<typename MT, typename VT>
    std::decay_t<decltype(~(std::declval<VT>()))> LinearSolve(const TMatrix<MT>& A, const TMatrix<VT>& b)
    {
        assert((~A).Rows() == (~A).Cols());
        assert((~A).Cols() == (~b).Rows());
        assert((~b).Cols() == 1);
        
        std::decay_t<decltype(~(std::declval<MT>()))> U(A);
        std::decay_t<decltype(~(std::declval<VT>()))> y(b);

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
            }

            // eliminate
            for (std::size_t i = k + 1; i < (~U).Rows(); i++)
            {
                auto lambda = (~U)(i, k) / (~U)(k, k);
                for (std::size_t j = k; j < (~U).Rows(); j++)
                    (~U)(i, j) -= lambda * (~U)(k, j);
                (~y)(i, 0) -= lambda * (~y)(k, 0);
            }
        }

        // U is in upper triagonal form (lower triganonal is not set to zero expilicitly for efficiency)
        // y is the original y transformed by the L matrix from the L-U decomposition
        // Now one can solve for x
        std::decay_t<decltype(~(std::declval<VT>()))> x = (~y);
        for (std::size_t l = 0; l < (~U).Rows(); l++)
        {
            std::size_t i = (~U).Rows() - l - 1;
            for (std::size_t j = i + 1; j < (~U).Rows(); j++)
                (~x)(i, 0) -= (~U)(i, j) * (~x)(j, 0);
            (~x)(i, 0) /= (~U)(i, i);
        }

        return (~x);
    }

}

#endif
