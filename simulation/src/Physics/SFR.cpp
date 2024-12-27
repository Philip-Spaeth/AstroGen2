#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"
#include <random>

double randomUniform()
{
    static thread_local std::mt19937_64 rng(42); // fixierter Seed nur als Beispiel
    static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

void SFR::sfrRoutine(Particle* particle)
{
    if (particle->type != 2) return;

    double DensityThreshold = 1.0e-22;
    double TemperatureThreshold = 1.0e4;
    double CStar            = 0.1;

    double rho  = particle->rho;   // lokale Gasdichte
    double mass = particle->mass;  // Masse des Gaspartikels

    if (rho < DensityThreshold)
        return;
    
    if (particle->T > TemperatureThreshold)
        return;

    double tDyn = std::sqrt(3.0 * M_PI / (16.0 * Constants::G * particle->rho));

    double sfrMass = CStar * (mass / tDyn);
    particle->sfr = sfrMass;

    double pForm = 1.0 - std::exp(- (sfrMass * particle->timeStep) / mass);

    // Zufallszahl ziehen
    double r = randomUniform();

    if (r < pForm)
    {
        //std::cout << "Star Formation!" << std::endl;
        particle->type = 1;
        particle->U = 0.0;
    }
}