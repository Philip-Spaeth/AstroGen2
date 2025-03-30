#ifndef COOLING_H
#define COOLING_H

#include "Constants.h"
#include <cmath>
#include <iostream>
#include <iomanip>

struct Particle;

struct IonizationResult {
    double n_H0;
    double n_Hplus;
    double n_He0;
    double n_Heplus;
    double n_Heplusplus;
    double n_e;
};

struct CoolingResult {
    double collExc_H = 0;         // Collisional excitation (H°)
    double collExc_HePlus = 0;    // Collisional excitation (He+)
    double collIon_H = 0;         // Collisional ionization (H°)
    double collIon_He0 = 0;       // Collisional ionization (He0)
    double collIon_HePlus = 0;    // Collisional ionization (He+)
    double recomb_H = 0;          // Recombination (H+)
    double recomb_HePlus = 0;     // Recombination (He+)
    double recomb_Heplusplus = 0; // Recombination (He++)
    double dielRecomb_HePlus = 0; // Dielectric recombination (He+)
    double freeFree = 0;          // Free-free (bremsstrahlung)
    double invCompton = 0;        // Inverse Compton cooling
    double total = 0;             // Total cooling rate
};

class Cooling {
public:
    static IonizationResult computeIonization(
        double n_H, double T,
        double Gamma_H0, double Gamma_He0, double Gamma_Heplus,
        double rci_H0, double rci_He0, double rci_Heplus,
        double alpha_Hplus, double alpha_Heplus, double alpha_Heplusplus
    );

    static double coolingExcitation_H(double T, double n_e, double n_H0);
    static double coolingExcitation_HePlus(double T, double n_e, double n_Heplus);
    static double coolingIonization_H(double T, double n_e, double n_H0);
    static double coolingIonization_He0(double T, double n_e, double n_He0);
    static double coolingIonization_HePlus(double T, double n_e, double n_Heplus);
    static double coolingRecombination_H(double T, double n_e, double n_Hplus);
    static double coolingRecombination_HePlus(double T, double n_e, double n_Heplus);
    static double coolingRecombination_Heplusplus(double T, double n_e, double n_Heplusplus);
    static double coolingDielectricRecombination_HePlus(double T, double n_e, double n_Heplus);
    static double coolingFreeFree(double T, double n_e, double n_Hplus, double n_Heplus, double n_Heplusplus);
    static double coolingInverseCompton(double T, double n_e, double z);

    static CoolingResult computeCoolingRates(double T, double z, const IonizationResult &ion);

    static void coolingRoutine(Particle* particle);
};

#endif