// Philipp Neufeld, 2021

#ifndef QSim_DensityMatrix_H_
#define QSim_DensityMatrix_H_

#include <cstdint>
#include <string>
#include <map>

#include "Math/Matrix.h"

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
        TStaticDensityMatrix(std::map<std::string, std::size_t> levelNames)
            : MyParent(N, N), m_levelNames(levelNames) {}
        TStaticDensityMatrix(std::map<std::string, std::size_t> levelNames, 
            const TStaticMatrix<std::complex<double>, N, N>& densityMatrix)
            : MyParent(densityMatrix), m_levelNames(levelNames) {}

        TStaticDensityMatrix(const TStaticDensityMatrix&) = default;
        TStaticDensityMatrix& operator=(const TStaticDensityMatrix&) = default;
        
        double GetPopulation(const std::string& lvl) const;
        double GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const;

        const MyParent& GetMatrix() const { return static_cast<const MyParent&>(*this); }

    private:
        std::map<std::string, std::size_t> m_levelNames;
    };

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetPopulation(const std::string& lvl) const 
    { 
        std::size_t idx = m_levelNames.at(lvl);
        return std::real((*this)(idx, idx)); 
    }

    template<std::size_t N>
    double TStaticDensityMatrix<N>::GetAbsCoeff(const std::string& lvl1, const std::string& lvl2) const 
    { 
        std::size_t idx1 = m_levelNames.at(lvl1);
        std::size_t idx2 = m_levelNames.at(lvl2);
        return std::imag((*this)(idx1, idx2)); 
    }

}

#endif
