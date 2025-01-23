#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

void Cooling::coolingRoutine(Particle* particle)
{
    double density = particle->rho;
    double temperature = particle->T;// * 1e5; // convert from K to 10^5 K

    double nH  = Constants::X_H * ((density) / (Constants::prtn)) * 1e-6;
    double nHe = (Constants::X_He / 4.0) * ((density) / (Constants::prtn)) * 1e-6;
    double nH0        = 0.0;
    double nHplus     = nH - nH0;
    double nHeplus    = 0.0;
    double nHeplusplus= nHe;
    double ne = nHplus + 2.0 * nHeplusplus;
    
    //all ions
    // Bremsstrahlung cooling in erg * s^-1 * cm^-3
    double gff = 1.1 + 0.34 * exp(-pow(5.5 - log10(temperature), 2.0) / 3.0);
    double bremsstrahlung = 1.42e-27 * gff * sqrt(temperature) * (nHplus + nHeplus + 4.0*nHeplusplus) * ne;
    bremsstrahlung *= 1e-7; // convert from erg to J
    bremsstrahlung *= 1e6; // convert from cm^-3 to m^-3

    double netCoolingRate = bremsstrahlung;
    double dudt = - (netCoolingRate / density);
    //std::cout << "dudt: " << dudt << std::endl;
    particle->dUdt += dudt;
}