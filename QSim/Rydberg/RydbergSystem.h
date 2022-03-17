// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_RydbergSystem_H_
#define QSim_Rydberg_RydbergSystem_H_

#include <cstdint>

#include <Eigen/Dense>

#include "../Constants.h"
#include "../Math/Quad.h"
#include "../Math/Numerov.h"
#include "../Math/Wigner.h"

namespace QSim
{

    class RydbergSystem
    {
    public:

        RydbergSystem(double mass);
        virtual ~RydbergSystem() = default;

        static double GetQuantumDefectFromCoeffs(int n, const double* coeffs, int coeffCnt);
        static double GetExtrapolatedQuantumDefect(double defect, int l, int lnew);
        virtual double GetQuantumDefect(int n, int l, double j) const = 0;

        double GetEnergy(int n, int l, double j) const;

        double CorePotential(double r, int n, int l) const;
        double FSPotential(double r, int n, int l, double j) const;

        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWF(int n, int l, double j,
            double rInner, double rOuter, std::size_t steps) const;
        
        double GetDipoleMERad(int n1, int l1, double j1, int n2, int l2, double j2) const;
        double GetDipoleMEAng(int l1, double j1, double mj1, int l2, double j2, double mj2) const;
        double GetDipoleME(int n1, int l1, double j1, double m1, int n2, int l2, double j2, double m2) const;
        
    private:
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFTransformed(int n, int l, double j,
            double xInner, double xOuter, std::size_t steps) const;

    private:
        double m_reducedMass;
        std::vector<double> m_quantumDefectCoeffs;

    };

}

#endif
