// Philipp Neufeld, 2021-2022

//! \file Wigner.h
//! This file contains functions that help with the calculation of the
//! Wigner-3nj symbols.
//!
//! The Wigner symbols can be expanded as fractions of products of factorials.
//! The quotients of these fracctions can get quite large such that one has to
//! handle fractions of two large numbers that result in numbers that are small
//! Calculating these fractions directly makes the procedure very inprecise.
//! This problem can be solved by utilizing the following identities:
//!   1. \f$ \Gamma(n) = (n-1)! \f$
//!   2. \f$ a \cdot b = \exp(\log(a) + \log(b)) \f$
//!   3. \f$ a^b = \exp(b \cdot \log(a)) \f$
//!
//! It is possible to directly calculate the logarithm of the gamma function
//! and thus the logarithm of the factorials. By using identity 2 one can
//! do the whole calculation logarithmically and exponentiating it in the end.
//! 

#ifndef QSim_Math_Wigner_H
#define QSim_Math_Wigner_H

#include <cmath>
#include "Gamma.h"

namespace QSim
{

    //! This function calculates the Wigner-3j symbol
    //! \f[
    //!     \left(\matrix{j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3} \right)
    //! \f]
    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);
    
    //! This function calculates the Wigner-6j symbol
    //! \f[
    //!     \left\{\matrix{j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6} \right\}
    //! \f]
    double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6);

    //! This function calculates the clebsh-gordan coefficients \f$ C(j_1, m_{j_1}, j_2, m_{j_2}; j_3, m_{j_3}) \f$
    //! \f[
    //!     \left.\vert j_1, j_2, j_3, m_{j_3} \right\rangle = \sum_{j_3, m_{j_3}} \left.\vert j_1, m_{j_3}, j_2, m_{j_2} \right\rangle
    //!     \left\langle j_1, m_{j_3}, j_2, m_{j_2} \right.\vert \left. j_1, j_2, j_3, m_{j_3} \right\rangle
    //!     = \sum_{j_3, m_{j_3}} C(j_1, m_{j_1}, j_2, m_{j_2}; j_3, m_{j_3}) \left.\vert j_1, m_{j_3}, j_2, m_{j_2} \right\rangle.
    //! \f]
    //! The clebsh-gordan coefficients are calculated by using the Wigner-3j symbols internally,
    //! \f[
    //!     C(j_1, m_{j_1}, j_2, m_{j_2}; j_3, m_{j_3}) = (-1)^{j_1 - j_2 + m_{j_3}} \sqrt{2j + 1} 
    //!     \left(\matrix{j_1 & j_2 & j_3 \\ m_1 & m_2 & -m_3} \right)
    //! \f]    
    //! \sa Wigner3j
    double ClebshGordan(double j1, double j2, double j3, 
        double m1, double m2, double m3);
}

#endif 
