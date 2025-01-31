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

void SFR::sfrRoutine(Particle* particle, Simulation* sim, std::vector<Particle*>* newStars)
{
    if(!particle) return;
    if(particle->rho <= 0) return;
    if(particle->U <= 0) return;
    if(particle->mass <= 0) return;

    double densityThreshold = 1e-23;
    double temperatureThreshold = 1e4;

    //if (particle->rho < densityThreshold || particle->T > temperatureThreshold) return;

    double c = 1;
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
        double massToConvert = fraction * particle->mass;
        particle->accumulatedStarMass += massToConvert;
        particle->mass -= massToConvert;

        if (particle->mass < 1e-10) {
            particle->mass = 0.0;
        }
    }

    double minMass = 1e35;

    if (particle->accumulatedStarMass >= minMass)
    {
        Particle* newStar = new Particle();
        newStar->type     = 1;
        newStar->mass     = particle->accumulatedStarMass;
        newStar->position = particle->position;
        newStar->velocity = particle->velocity;
        
        #pragma omp critical
        {
            newStars->push_back(newStar);
        }
        particle->accumulatedStarMass = 0.0;

        particle->mass -= newStar->mass;
        if (particle->mass < 0.0) {
            particle->mass = 0.0;
        }
    }
}