#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"
#include <random>

double randomUniform() 
{
    return static_cast<double>(rand()) / RAND_MAX;
}

void SFR::sfrRoutine(Particle* particle)
{
    if(!particle) return;
    if(particle->rho <= 0) return;
    if(particle->U <= 0) return;
    if(particle->mass <= 0) return;

    //density threshold in SI units
    double densityThreshold = 1e-23;
    //temperature threshold in SI units
    double temperatureThreshold = 1e4;

    if (particle->rho < densityThreshold || particle->T > temperatureThreshold) return;

    double c = 0.01;
    double tdyn = sqrt(3 * M_PI / (32 * Constants::G * particle->rho));

    double sfr = c * particle->rho / tdyn;
    double Vpart = particle->mass / particle->rho;
    particle->sfr = ((sfr * Vpart) / 1.98847e30) * 3.154e7; // in solar masses per year
    totalSFR += particle->sfr;
    double dstarMass = sfr * particle->timeStep * Vpart;
    double P = dstarMass / particle->mass;

    if (randomUniform() < P)
    {
        particle->type = 1;
        particle->U = 0;
        particle->dUdt = 0;
        particle->sfr = 0;
        particle->T = 0;
        //std::cout << "Star formation at particle " << particle->id << std::endl;
    }
}