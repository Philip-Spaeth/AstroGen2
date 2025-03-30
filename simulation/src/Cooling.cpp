#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Particle.h"

CoolingResult Cooling::computeCoolingRates(double T, double z, const IonizationResult &ion) {
    CoolingResult cool;
    if (T < 1e3) return cool;
    
    cool.collExc_H = coolingExcitation_H(T, ion.n_e, ion.n_H0);
    cool.collExc_HePlus = coolingExcitation_HePlus(T, ion.n_e, ion.n_Heplus);
    cool.collIon_H = coolingIonization_H(T, ion.n_e, ion.n_H0);
    cool.collIon_He0 = coolingIonization_He0(T, ion.n_e, ion.n_He0);
    cool.collIon_HePlus = coolingIonization_HePlus(T, ion.n_e, ion.n_Heplus);
    cool.recomb_H = coolingRecombination_H(T, ion.n_e, ion.n_Hplus);
    cool.recomb_HePlus = coolingRecombination_HePlus(T, ion.n_e, ion.n_Heplus);
    cool.recomb_Heplusplus = coolingRecombination_Heplusplus(T, ion.n_e, ion.n_Heplusplus);
    cool.dielRecomb_HePlus = coolingDielectricRecombination_HePlus(T, ion.n_e, ion.n_Heplus);
    cool.freeFree = coolingFreeFree(T, ion.n_e, ion.n_Hplus, ion.n_Heplus, ion.n_Heplusplus);
    cool.invCompton = coolingInverseCompton(T, ion.n_e, z);
    
    cool.total = cool.collExc_H + cool.collExc_HePlus + cool.collIon_H + cool.collIon_He0 + 
                cool.collIon_HePlus + cool.recomb_H + cool.recomb_HePlus + cool.recomb_Heplusplus + 
                cool.dielRecomb_HePlus + cool.freeFree + cool.invCompton;
    
    return cool;
}

void Cooling::coolingRoutine(Particle* particle) {
    if (particle->delayedCoolingTime > 0) {
        particle->delayedCoolingTime -= particle->timeStep;
        return;
    } else {
        particle->delayedCoolingTime = 0;
    }

    double n_H = 1.0;
    double T = particle->T;
    double z = 0.0;
    double Gamma_H0 = 1e-12;    
    double Gamma_He0 = 1e-12;
    double Gamma_Heplus = 1e-12;
    double rci_H0 = 1e-10;
    double rci_He0 = 1e-10;
    double rci_Heplus = 1e-10;
    double alpha_Hplus = 1e-13;
    double alpha_Heplus = 1e-13;
    double alpha_Heplusplus = 1e-13;

    IonizationResult ion = computeIonization(
        n_H, T, Gamma_H0, Gamma_He0, Gamma_Heplus,
        rci_H0, rci_He0, rci_Heplus,
        alpha_Hplus, alpha_Heplus, alpha_Heplusplus
    );

    CoolingResult cool = computeCoolingRates(T, z, ion);
    double lambda_cgs = cool.total;
    double lambda_SI = lambda_cgs * 1e-7;
    double dudt = - (lambda_SI / particle->rho);
    particle->dUdt += dudt;
}

IonizationResult Cooling::computeIonization(
    double n_H, double T,
    double Gamma_H0, double Gamma_He0, double Gamma_Heplus,
    double rci_H0, double rci_He0, double rci_Heplus,
    double alpha_Hplus, double alpha_Heplus, double alpha_Heplusplus
) {
    IonizationResult result;

    double ratio_H = alpha_Hplus / (rci_H0 + Gamma_H0);
    result.n_Hplus = n_H / (1.0 + ratio_H);
    result.n_H0 = n_H - result.n_Hplus;

    double n_He = (0.24 / (4.0 * 0.76)) * n_H;  // ca. 0.0789*n_H
    double factor0 = (rci_He0 + Gamma_He0) / alpha_Heplus;
    double factor1 = (rci_Heplus + Gamma_Heplus) / alpha_Heplusplus;
    double denom = 1.0 + factor0 + factor0 * factor1;
    result.n_He0 = n_He / denom;
    result.n_Heplus = result.n_He0 * factor0;
    result.n_Heplusplus = result.n_He0 * factor0 * factor1;
    result.n_e = result.n_Hplus + result.n_Heplus + 2.0 * result.n_Heplusplus;

    return result;
}

double Cooling::coolingExcitation_H(double T, double n_e, double n_H0) {
    return 7.50e-19 * pow(T, -0.5) * exp(-118348.0 / T) * n_e * n_H0;
}

double Cooling::coolingExcitation_HePlus(double T, double n_e, double n_Heplus) {
    return 5.54e-17 * pow(T, -0.397) * exp(-4736380.0 / T) / (1.0 + pow(T, 0.5)) * n_e * n_Heplus;
}

double Cooling::coolingIonization_H(double T, double n_e, double n_H0) {
    return 1.27e-21 * pow(T, 0.5) * exp(-157809.0 / T) * n_e * n_H0;
}

double Cooling::coolingIonization_He0(double T, double n_e, double n_He0) {
    return 9.38e-22 * pow(T, 0.5) * exp(-285335.0 / T) * n_e * n_He0;
}

double Cooling::coolingIonization_HePlus(double T, double n_e, double n_Heplus) {
    return 4.95e-22 * pow(T, 0.5) * exp(-285335.0 / T) * n_e * n_Heplus;
}

double Cooling::coolingRecombination_H(double T, double n_e, double n_Hplus) {
    return 8.70e-27 * pow(T, 0.5) * n_e * n_Hplus;
}

double Cooling::coolingRecombination_HePlus(double T, double n_e, double n_Heplus) {
    return 1.55e-26 * n_e * n_Heplus;
}

double Cooling::coolingRecombination_Heplusplus(double T, double n_e, double n_Heplusplus) {
    return 3.48e-26 * pow(T, 0.5) * n_e * n_Heplusplus;
}

double Cooling::coolingDielectricRecombination_HePlus(double T, double n_e, double n_Heplus) {
    return 1.24e-21 * n_e * n_Heplus;
}

double Cooling::coolingFreeFree(double T, double n_e, double n_Hplus, double n_Heplus, double n_Heplusplus) {
    double n_i = n_Hplus + n_Heplus + 4.0 * n_Heplusplus;
    double gaunt = 1.1 + 0.34 * exp(-pow(5.5 - log10(T), 2) / 3.0);
    return 1.42e-27 * gaunt * pow(T, 0.5) * n_e * n_i;
}

double Cooling::coolingInverseCompton(double T, double n_e, double z) {
    return 5.41e-36 * n_e * T * pow((1.0 + z), 4);
}