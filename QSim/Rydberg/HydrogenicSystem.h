// Philipp Neufeld, 2021-2022

#ifndef QSim_Rydberg_HydrogenicSystem_H_
#define QSim_Rydberg_HydrogenicSystem_H_

#include <cstdint>
#include <utility>
#include <cmath>
#include <complex>
#include <algorithm>

#include <Eigen/Dense>

#include <iostream>
#include "../Constants.h"
#include "../Math/Quad.h"
#include "../Math/Numerov.h"
#include "../Math/Wigner.h"

namespace QSim
{

    class HydrogenicSystem
    {
    public:
        double GetEnergy(int n) const;

        double CorePotential(double r, int n, int l) const;
        double FSPotential(double r, int n, int l, double j) const;

        /*std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFLinear(int n, int l, 
            double rInner, double rOuter, std::size_t steps) const;*/
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWF(int n, int l, double j,
            double rInner, double rOuter, std::size_t steps) const;
        
        double GetDipoleMERad(int n1, int l1, double j1, int n2, int l2, double j2) const;
        double GetDipoleMEAng(int l1, int m1, int l2, int m2) const;
        double GetDipoleMEAng2(int l1, double j1, double mj1, int l2, double j2, double mj2) const;
        //double GetDipoleME(int n1, int l1, int m1, int n2, int l2, int m2) const;
        double GetDipoleME2(int n1, int l1, double j1, double m1, int n2, int l2, double j2, double m2) const;
        
    public: // TODO: make private
        std::pair<Eigen::VectorXd, Eigen::VectorXd> GetRadialWFTransformed(int n, int l, double j,
            double xInner, double xOuter, std::size_t steps) const;
        

    };

    namespace Internal
    {
        class HydSysNumerovMonitor
        {
        public:
            HydSysNumerovMonitor(int distThreshold) 
                : m_distThreshold(distThreshold), 
                m_maxIdx(0), m_maxVal(0.0) {};

            template<typename Xs, typename Ys>
            bool operator()(int i, Xs& xs, Ys& ys)
            {
                double absYi = std::abs(ys[i]);
                if (i - m_maxIdx < m_distThreshold)
                {
                    // get height of first peak
                    if (absYi > m_maxVal)
                    {
                        m_maxVal = absYi;
                        m_maxIdx = i;
                    }
                }
                else
                {
                    // subsequent peaks may not be higher than the first peak
                    if (absYi > m_maxVal)
                    {
                        // clean up (set to zero until node where divergence started)
                        for(;i > 1 && std::abs(ys[i]) > std::abs(ys[i-1]); i--);
                        while(i < ys.size()) ys[i++] = 0;
                        return false;
                    }
                }

                return true;
            }

        private:
            const int m_distThreshold;
            int m_maxIdx;
            double m_maxVal;
        };
    }
    
}

#endif
