// Philipp Neufeld, 2021

#ifndef QSim_DensityMatrix_H_
#define QSim_DensityMatrix_H_

#include <cstdint>
#include <string>
#include <map>

#include "../Math/Matrix.h"

namespace QSim
{
    //
    // Density matrix class
    //

    template<std::size_t N>
    class TStaticDensityMatrix : public TStaticMatrix<std::complex<double>, N, N>
    {
        using MyParent = TStaticMatrix<std::complex<double>, N, N>;
    public:
        TStaticDensityMatrix(std::array<std::string, N> levelNames)
            : MyParent(N, N), m_levelNames(levelNames) {}
        TStaticDensityMatrix(std::array<std::string, N> levelNames, 
            const TStaticMatrix<std::complex<double>, N, N>& densityMatrix)
            : MyParent(densityMatrix), m_levelNames(levelNames) {}

        TStaticDensityMatrix(const TStaticDensityMatrix&) = default;
        TStaticDensityMatrix& operator=(const TStaticDensityMatrix&) = default;
        
        double GetPopulation(const std::string& lvl) const;
        double GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const;

        const std::array<std::string, N>& GetLevelNames() const { return m_levelNames; }
        const MyParent& GetMatrix() const { return static_cast<const MyParent&>(*this); }

    private:
        std::array<std::string, N> m_levelNames;
    };

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetPopulation(const std::string& lvl) const 
    { 
        auto it = std::find(m_levelNames.begin(), m_levelNames.end(), lvl);
        if (it == m_levelNames.end())
            return 0.0;
        std::size_t idx = it - m_levelNames.begin();
        return std::real((*this)(idx, idx)); 
    }

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const 
    { 
        auto it1 = std::find(m_levelNames.begin(), m_levelNames.end(), lvl1);
        if (it1 == m_levelNames.end())
            return 0.0;
        std::size_t idx1 = it1 - m_levelNames.begin();

        auto it2 = std::find(m_levelNames.begin(), m_levelNames.end(), lvl2);
        if (it2 == m_levelNames.end())
            return 0.0;
        std::size_t idx2 = it2 - m_levelNames.begin();

        return std::imag((*this)(idx1, idx2)); 
    }

}

#endif
