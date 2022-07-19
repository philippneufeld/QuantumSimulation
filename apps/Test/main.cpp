// Philipp Neufeld, 2021-2022

#include <iostream>
#include <Eigen/Dense>

#include <QSim/Rydberg/RydbergDiatomic.h>

using namespace QSim;
using namespace Eigen;

constexpr std::array<double, 81> g0 = { 0.99999, 0.99710, 0.98845, 0.97415, 0.95436, 0.92933, 0.89933, 0.86472, 0.82590, 0.78331, 0.73743, 0.68878, 0.63789, 0.58533, 0.53164, 0.47741, 0.42318, 0.36949, 0.31687, 0.26581, 0.21676, 0.17013, 0.12631, 0.08561, 0.04830, 0.01459, -0.01536, -0.04145, -0.06363, -0.08192, -0.09638, -0.10712, -0.11430, -0.11814, -0.11886, -0.11674, -0.11209, -0.10521, -0.09644, -0.08614, -0.07463, -0.06227, -0.04939, -0.03630, -0.02331, -0.01071, 0.00127, 0.01240, 0.02249, 0.03139, 0.03896, 0.04513, 0.04985, 0.05309, 0.05488, 0.05526, 0.05430, 0.05209, 0.04877, 0.04447, 0.03934, 0.03354, 0.02725, 0.02062, 0.01384, 0.00707, 0.00046, -0.00585, -0.01171, -0.01701, -0.02166, -0.02558, -0.02871, -0.03101, -0.03246, -0.03307, -0.03287, -0.03187, -0.03016, -0.02779, -0.02485 };
constexpr std::array<double, 81> g1 = { 0, -0.04157, -0.08258, -0.12249, -0.16076, -0.19688, -0.23040, -0.26088, -0.28795, -0.31128, -0.33061, -0.34574, -0.35655, -0.36296, -0.36498, -0.36269, -0.35621, -0.34576, -0.33158, -0.31399, -0.29334, -0.27004, -0.24450, -0.21719, -0.18857, -0.15913, -0.12934, -0.09967, -0.07058, -0.04250, -0.01582, 0.00909, 0.03191, 0.05237, 0.07025, 0.08537, 0.09763, 0.10696, 0.11335, 0.11686, 0.11759, 0.11567, 0.11129, 0.10469, 0.09612, 0.08587, 0.07423, 0.06153, 0.04810, 0.03426, 0.02032, 0.00660, -0.00661, -0.01905, -0.03047, -0.04068, -0.04949, -0.05677, -0.06242, -0.06639, -0.06865, -0.06922, -0.06816, -0.06555, -0.06151, -0.05619, -0.04976, -0.04241, -0.03433, -0.02574, -0.01686, -0.00790, 0.00094, 0.00945, 0.01744, 0.02476, 0.03124, 0.03678, 0.04126, 0.04462, 0.04682 };
constexpr std::array<double, 81> g2 = { 0.0005, 0.0017, 0.0053, 0.0112, 0.0193, 0.0294, 0.0414, 0.0548, 0.0695, 0.0850, 0.1012, 0.1175, 0.1337, 0.1493, 0.1641, 0.1776, 0.1896, 0.1998, 0.2078, 0.2137, 0.2170, 0.2178, 0.2160, 0.2116, 0.2045, 0.1950, 0.1830, 0.1690, 0.1529, 0.1352, 0.1162, 0.0960, 0.0752, 0.0541, 0.0329, 0.0122, -0.0079, -0.0269, -0.0445, -0.0604, -0.0745, -0.0865, -0.0962, -0.1035, -0.1084, -0.1109, -0.1109, -0.1085, -0.1039, -0.0973, -0.0887, -0.0786, -0.0671, -0.0544, -0.0410, -0.0272, -0.0131, 0.0008, 0.0143, 0.0270, 0.0388, 0.0493, 0.0584, 0.0660, 0.0718, 0.0759, 0.0781, 0.0785, 0.0770, 0.0739, 0.0691, 0.0628, 0.0553, 0.0466, 0.0371, 0.0270, 0.0164, 0.0057, -0.0048, -0.0151, -0.0247 };
constexpr std::array<double, 81> g3 = { 0, -0.017, -0.033, -0.050, -0.066, -0.081, -0.096, -0.110, -0.123, -0.136, -0.147, -0.157, -0.166, -0.173, -0.180, -0.184, -0.188, -0.189, -0.190, -0.188, -0.186, -0.181, -0.176, -0.168, -0.160, -0.150, -0.140, -0.128, -0.115, -0.102, -0.088, -0.074, -0.060, -0.045, -0.031, -0.017, -0.003, 0.010, 0.022, 0.033, 0.043, 0.052, 0.060, 0.066, 0.071, 0.074, 0.076, 0.076, 0.075, 0.073, 0.069, 0.064, 0.058, 0.051, 0.044, 0.035, 0.026, 0.017, 0.008, -0.001, -0.010, -0.018, -0.026, -0.033, -0.039, -0.045, -0.049, -0.052, -0.054, -0.055, -0.055, -0.053, -0.051, -0.047, -0.042, -0.037, -0.031, -0.024, -0.017, -0.010, -0.003};

double GetMatrixElement(double n1, double l1, double n2, double l2)
{
    if (n2 > n1)
        return GetMatrixElement(n2, l2, n1, l1);

    double lc = std::max(l1, l2);
    double nc = 2.0/(1.0/n1 + 1.0/n2);
    double dl = l2 - l1;
    double gamma = dl * lc / nc;
    double s = n1 - n2;

    int idx = static_cast<int>(std::floor(s / 0.05));
    double t = s - idx*0.05;

    double result = 0.0;
    std::array<double, 81> gcoeffs[] = {g0, g1, g2, g3};
    for (int i=0; i<3; i++)
        result += std::pow(gamma, i) * ((1-t)*gcoeffs[i][idx] + t*gcoeffs[i][idx+1]);
    
    result *= 1.5*nc*nc*std::sqrt(1 - (lc/nc)*(lc/nc));

    return result * BohrRadius_v;
}

int main(int argc, const char *argv[])
{
    NitricOxide molecule;

    int n1 = 33;
    int n2 = 33;
    int l1 = 1;
    int l2 = 2;

    auto st1 = std::make_tuple(n1, l1, 0, 1, 0);
    auto st2 = std::make_tuple(n2, l2, 0, 1, 0);
    double rad = molecule.GetDipMeRadHelper(
        st1, st2, 1,
        molecule.GetIntegrationRange(n1), 
        molecule.GetIntegrationRange(n2), 100);
    
    std::cout << rad << std::endl;

    double n1Eff = n1 - molecule.GetQuantumDefect(st1);
    double n2Eff = n2 - molecule.GetQuantumDefect(st2);
    double rad2 = GetMatrixElement(n1Eff, l1, n2Eff, l2);
    std::cout << rad2 << std::endl;


    int i=0, j=0;
    for (int R1=0; R1<=2; R1+=2)
    {
        for (int l1=0; l1<=2; l1++)
        {
            for (int N1=std::abs(R1-l1); N1<=R1+l1; N1++, i++)
            {
               for (int R2=0; R2<=2; R2+=2)
                {
                    for (int l2=0; l2<=2; l2++)
                    {
                        for (int N2=std::abs(R2-l2); N2<=R2+l2; N2++, j++)
                        {
                            auto st1 = std::make_tuple(n1, l1, R1, N1, 0);
                            auto st2 = std::make_tuple(n1, l2, R2, N2, 0);
    
                            double h = molecule.GetCoreInteractionME(st1, st2);
                            h += molecule.GetDipoleME(st1, st2);
                            std::cout << (h != 0) << " " << std::flush;
                        }
                    }
                    std::cout << "| ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << "______|___________________|" << std::endl;
    }
    std::cout << std::endl;

    // H state -> 3d state
    // C. Jungen, Rydberg Series in the NO Spectrum: An Interpretation of Quantum Defects and Intensities in the sand d Series
    constexpr double hc = SpeedOfLight_v * PlanckConstant_v;
    double threshold = 74720;
    double B = 1e-2 * molecule.GetRotationalConstant();
    
    double Jx = 5.5;
    double Ja = Jx - 1;
    double Jh = Ja + 1;

    double xstate = (-threshold + B*Jx*(Jx+1) - B*0.5*1.5) * EnergyInverseCm_v;
    double astate = (-threshold + 44200 + B*Ja*(Ja+1)) * EnergyInverseCm_v;
    double hstate = (-threshold + 62705 + B*Jh*(Jh+1)) * EnergyInverseCm_v;

    astate = xstate + hc / (907.88e-9/4);
    hstate = astate + hc / (1080.92e-9/2);

    std::cout << "X(J=" << Jx << "): " << xstate / EnergyInverseCm_v << " cm^-1" << std::endl;
    std::cout << "A(J=" << Ja << "): " << astate / EnergyInverseCm_v << " cm^-1" << std::endl;
    std::cout << "H(J=" << Jh << "): " << hstate / EnergyInverseCm_v << " cm^-1" << std::endl;

    std::cout << "A<-X: " << hc / (astate - xstate) * 1e9 << " nm" << std::endl;
    std::cout << "H<-A: " << hc / (hstate - astate) * 1e9 << " nm" << std::endl;
    std::cout << std::endl;

    std::vector<std::tuple<double, int, int>> quantumNumbers;

    for (int n = 4; n<=100; n++)
    {
        for (int R=0; R<=8; R++)
        {
            double rydberg = molecule.GetEnergy(std::make_tuple(n, 3, R, 0, 0));
            double lambda = hc / (rydberg - hstate) * 1e9;
            quantumNumbers.emplace_back(lambda, n, R);

        }
    }

    std::sort(quantumNumbers.begin(), quantumNumbers.end(), [=](auto x1, auto x2) 
    {
        constexpr double redLaserWL = 834.92;
        auto [l1, n1, R1] = x1;
        auto [l2, n2, R2] = x2;
        return std::abs(l1-redLaserWL) < std::abs(l2-redLaserWL);
    });

    std::cout << "Possible rydberg states:" << std::endl;
    for (int i=0; i < 15 && i < quantumNumbers.size(); i++)
    {
        auto [lambda, n, R] = quantumNumbers[i];
        std::cout << "n=" << n << " R=" << R << " (" << lambda << "nm)" << std::endl;
    }

    return 0;
}
