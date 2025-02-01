#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"
#include <random>
#include "Simulation.h"

double randomUniform() 
{
    return static_cast<double>(rand()) / RAND_MAX;
}

void SFR::sfrRoutine(Particle* particle, Simulation* sim, double& newStarMass)
{
    if(!particle) return;
    if(particle->rho <= 0) return;
    if(particle->U <= 0) return;
    if(particle->mass <= 0) return;

    double densityThreshold = 1e-23;
    double temperatureThreshold = 1e4;

    //if (particle->rho < densityThreshold || particle->T > temperatureThreshold) return;

    double c = 0.1;
    double tdyn = sqrt(3 * M_PI / (32 * Constants::G * particle->rho));
    double dt = particle->timeStep;
    double mdot = c * particle->mass / tdyn;
    double dM   = mdot * dt;

    double fraction = dM / particle->mass;
    if (fraction > 1.0) {
        fraction = 1.0;
    }

    double p = 1.0 - std::exp(-fraction);
    double r = randomUniform();
    if (r < p)
    {
        // convert gas particle to star particle
        particle->type = 1;
        newStarMass += particle->mass;

        //Supernoavae
        double snChance = 0.12; 
        double rSN = randomUniform();
        if (rSN < snChance)
        {
            //convert to hot gas
            particle->type = 2;
            particle->T = 1e7;
            particle->U = Constants::k_b * particle->T / ((Constants::GAMMA - 1.0) * Constants::prtn * particle->mu);
        }
    }
}