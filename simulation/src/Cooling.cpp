#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>
#include <iomanip>

void Cooling::coolingRoutine(Particle* particle)
{
    double density = particle->rho;
    double T = particle->T;

    double nH_m3  = Constants::X_H * (density / Constants::prtn);  // [m^-3]
    double nH_cm3 = nH_m3 * 1e-6;                                  // [cm^-3]
    double nHe_m3  = (Constants::X_He / 4.0) * (density / Constants::prtn); // [m^-3]
    double nHe_cm3 = nHe_m3 * 1e-6; 

    CIEFractions fr = getCIEFractions(T, nH_cm3, nHe_cm3);

    double Lambda_cgs = coolingRatePrimordialCgs(
        T, fr.n_e,
        fr.n_H0, fr.n_Hplus,
        fr.n_He0, fr.n_Heplus, fr.n_Heplusplus, particle
    );

    double netCoolingRateSI = Lambda_cgs * 1e-7 * 1e6; // [erg cm^-3 s^-1] -> [J m^-3 s^-1]
    double dudt = - (netCoolingRateSI / density);
    //std::cout << std::scientific << std::setprecision(6);
    //std::cout << "lambda: " << particle->T << std::endl;
    particle->dUdt += dudt;
}


double Cooling::coolingRatePrimordialCgs(
    double T,         
    double n_e,        
    double n_H0,       
    double n_Hplus,    
    double n_He0,      
    double n_Heplus,   
    double n_Heplusplus,
    Particle* particle
)
{
    // collisional excitation
    double lambda_ce_H0     = lambdaExcitationH0(T)     * n_e * n_H0;
    double lambda_ce_Heplus = lambdaExcitationHeplus(T) * n_e * n_Heplus;
    double Lambda_exc = lambda_ce_H0 + lambda_ce_Heplus;

    // collisional ionisation
    const double eV_in_erg = 1.602176634e-12; // 1 eV ~ 1.60e-12 erg
    double E_H0     = 13.6  * eV_in_erg;
    double E_He0    = 24.6  * eV_in_erg;
    double E_Heplus = 54.4  * eV_in_erg;

    double ci_H0       = kIonH0(T)      * n_e * n_H0;      // [1/s cm^-3]
    double ci_He0      = kIonHe0(T)     * n_e * n_He0;
    double ci_Heplus   = kIonHeplus(T)  * n_e * n_Heplus;

    double Lambda_ci = ci_H0      * E_H0
                     + ci_He0     * E_He0
                     + ci_Heplus  * E_Heplus;
    double rec_Hplus      = alphaHplus(T)         * n_e * n_Hplus;   // [1/s cm^-3]
    double rec_Heplus     = alphaHeplus(T)        * n_e * n_Heplus;
    double rec_Heplusplus = alphaHeplusplus(T)    * n_e * n_Heplusplus;
    double dielec_Heplus  = alphaDielecHeplus(T)  * n_e * n_Heplus;
    double epsRec_H   = 0.75 * E_H0;
    double epsRec_He0 = 0.75 * E_He0;
    double epsRec_He1 = 0.75 * E_Heplus;

    double Lambda_rec = rec_Hplus      * epsRec_H
                      + rec_Heplus     * epsRec_He0
                      + rec_Heplusplus * epsRec_He1
                      + dielec_Heplus  * epsRec_He0;

    //free free cooling
    double gff  = gauntFactorFF(T);
    double ff   = 1.42e-27 * gff * std::sqrt(T)
                 * (n_Hplus + n_Heplus + 4.0*n_Heplusplus)
                 * n_e;


    double Lambda_total = Lambda_exc + Lambda_ci + Lambda_rec + ff;
    return Lambda_total;
}

inline double Cooling::kIonH0(double T)  // Coll. Ionisation H^0 -> H^+
{
    // 1.27e-21 sqrt(T) e^(-157809/T) / (1 + sqrt(T))
    // cgs: returns [cm^3 / s]
    double sT = std::sqrt(T);
    return 1.27e-21 * sT * std::exp(-157809.1 / T) / (1.0 + sT);
}


inline double Cooling::kIonHe0(double T) // Coll. Ionisation He^0 -> He^+
{
    // 9.38e-22 sqrt(T) e^(-285335.4 / T) / (1 + sqrt(T))
    return 9.38e-22 * std::sqrt(T)
                   * std::exp(-285335.4 / T)
                   / (1.0 + std::sqrt(T));
}
CIEFractions Cooling::getCIEFractions(
        double T,
        double n_H,
        double n_He
)
{
    double xHplus        = 0.99;
    double xHeplus       = 0.01; 
    double xHeplusplus   = 0.88;
    const int    maxIter = 50;
    const double tol     = 1e-6; 

    for(int iter=0; iter<maxIter; ++iter)
    {
        double old_xHplus      = xHplus;
        double old_xHeplus     = xHeplus;
        double old_xHeplusplus = xHeplusplus;
        //    kIonH0(T)*n_e*n_H0 = alphaHplus(T)*n_e*n_Hplus
        //    => xHplus / (1 - xHplus) = kIonH0(T)/alphaHplus(T)
        //    => xHplus = ratio / (1 + ratio).
        double ratioH = kIonH0(T) / std::max(alphaHplus(T), 1e-30);
        double new_xHplus = ratioH / (1.0 + ratioH);

        //    1: He^0 <-> He^+
        //         kIonHe0(T)*n_e*n_He0 = [alphaHeplus + alphaDielecHeplus]*n_e*n_Heplus
        //         => xHeplus / xHe0 = kIonHe0 / (alphaHeplus + alphaDielecHeplus)
        //
        //    2: He^+ <-> He^{++}
        //         kIonHeplus(T)*n_e*n_Heplus = alphaHeplusplus(T)*n_e*n_He^{++}
        //         => xHeplusplus / xHeplus = kIonHeplus / alphaHeplusplus
        //
        //    mit xHe0 + xHeplus + xHeplusplus = 1
        //    => xHe0 = 1 - xHeplus - xHeplusplus.
        double kI0   = kIonHe0(T);
        double kI1   = kIonHeplus(T);
        double aHe1  = alphaHeplus(T) + alphaDielecHeplus(T);
        double aHe2  = alphaHeplusplus(T);

        // ratio1 = xHeplus / xHe0
        double ratio1 = 0.0;
        if(aHe1 > 0.0) ratio1 = kI0 / aHe1;
        // ratio2 = xHeplusplus / xHeplus
        double ratio2 = 0.0;
        if(aHe2 > 0.0) ratio2 = kI1 / aHe2;

        // xHeplus + xHe0 = xHeplus + (1 - xHeplus - xHeplusplus) = 1 - xHeplusplus
        // => xHeplus / (1 - xHeplusplus - xHeplus) = ratio1 => ...
        //    xHe0 = 1 - (xHeplus + xHeplusplus).
        //    xHeplus / xHe0 = ratio1 => xHeplus = ratio1 * xHe0
        //    => xHeplus = ratio1 * [1 - (xHeplus + xHeplusplus)]
        // xHeplusplus / xHeplus = ratio2 => xHeplusplus = ratio2 * xHeplus
        double denom = 1.0 + ratio1*(1.0 + ratio2);
        double new_xHeplus = (denom > 0.0) ? (ratio1 / denom) : 0.0;
        double new_xHeplusplus = ratio2 * new_xHeplus;
        double new_xHe0 = 1.0 - new_xHeplus - new_xHeplusplus;
        if(new_xHe0 < 0.0) new_xHe0 = 0.0;

        xHplus      = new_xHplus;
        xHeplus     = new_xHeplus;
        xHeplusplus = new_xHeplusplus;

        double err = std::fabs(xHplus - old_xHplus)
                   + std::fabs(xHeplus - old_xHeplus)
                   + std::fabs(xHeplusplus - old_xHeplusplus);
        if(err < tol) break;
    }
    CIEFractions out;
    out.n_Hplus       = xHplus      * n_H;
    out.n_H0          = (1.0 - xHplus) * n_H;
    out.n_Heplus      = xHeplus     * n_He;
    out.n_Heplusplus  = xHeplusplus * n_He;
    out.n_He0         = (1.0 - xHeplus - xHeplusplus) * n_He;
    out.n_e           = out.n_Hplus + out.n_Heplus + 2.0*out.n_Heplusplus;
    return out;
}

inline double Cooling::kIonHeplus(double T) // Coll. Ion. He^+ -> He^{++}
{
    // 4.95e-22 sqrt(T) e^(-631515.0 / T) / (1 + sqrt(T))
    return 4.95e-22 * std::sqrt(T)
                    * std::exp(-631515.0 / T)
                    / (1.0 + std::sqrt(T));
}
inline double Cooling::alphaHplus(double T) // Recomb H^+ -> H^0
{
    double sT = std::sqrt(T);
    return 8.70e-27 * sT * std::pow(T, -0.2) 
                    / (1.0 + std::pow(T/1.0e5, 0.7));
}

inline double Cooling::alphaHeplus(double T) // Recomb He^+ -> He^0
{
    return 1.55e-26 * std::pow(T, 0.3647);
}

inline double Cooling::alphaHeplusplus(double T) // Recomb He^{++} -> He^+
{
    // ~ 3.48e-26 T^0.5 ...
    return 3.48e-26 * std::pow(T, 0.5) * std::pow(T, -0.2); 
}
inline double Cooling::alphaDielecHeplus(double T)
{
    // 1.24e-13 T^-1.5 exp(-470000 / T) (1 + 0.3 e^-94000/T)
    return 1.24e-13 * std::pow(T, -1.5)
                     * std::exp(-470000.0 / T)
                     * (1.0 + 0.3 * std::exp(-94000.0 / T));
}

inline double Cooling::lambdaExcitationH0(double T)
{
    // 7.50e-19 T^-0.118 exp(-118348/T) / (1 + sqrt(T))
    // cgs => [erg cm^3 / s], multiplied by n_e n_H0 => [erg cm^-3 s^-1]
    double sT = std::sqrt(T);
    return 7.50e-19 * std::pow(T, -0.118)
                    * std::exp(-118348.0 / T)
                    / (1.0 + sT);
}

inline double Cooling::lambdaExcitationHeplus(double T)
{
    // 5.54e-17 T^0.397 exp(-4736380/T)/(1 + sqrt(T))
    return 5.54e-17 * std::pow(T, 0.397)
                    * std::exp(-4736380.0 / T)
                    / (1.0 + std::sqrt(T));
}
inline double Cooling::gauntFactorFF(double T)
{
    double logT = std::log10(T);
    return 1.1 + 0.34 * std::exp(-std::pow(5.5 - logT, 2.0) / 3.0);
}
