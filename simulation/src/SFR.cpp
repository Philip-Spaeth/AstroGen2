#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"
#include <random>
#include "Simulation.h"
#include <algorithm> 

SFR::SFR(){rng.seed(42);}

void SFR::sfrRoutine(Particle* particle, Simulation* sim, double& newStarMass)
{
    if(!particle) return;
    if(particle->rho <= 0) return;
    if(particle->U <= 0) return;
    if(particle->mass <= 0) return;

    //if T > 10000 Kelvin
    if (particle->T > 1e4);

    //2 × 10−25 g cm−3 or 0.1 Hydrogen atoms per cm^3
    double densityThreshold = 2e-22;
    if (particle->rho < densityThreshold) return;

    //convergent flow, calculted by the Tree, see node.cpp
    if (particle->div_v > 0) return;

    //jeans unstable
    double soundspeed = sqrt(Constants::GAMMA*(particle->P / particle->rho));
    double soundcrossingTime = particle->h / soundspeed;
    double jeansCrit = 1.0 / sqrt(4*Constants::PI*Constants::G*particle->rho);
    //if (soundcrossingTime < jeansCrit) return;

    double c = sim->c_sfr;
    double tdyn = sqrt(3 * M_PI / (32 * Constants::G * particle->rho));
    double dt = particle->timeStep;
    
    double p_star = 1.0 - std::exp(-c * dt / tdyn);
    
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(rng);
    
    if (r < p_star)
    {
        particle->type = 1;
        newStarMass += particle->mass;

        if(sim->SNFeedbackEnabled)
        {
            double p_supernova = 0.12;
            double r_sn = dist(rng);
            if (r_sn < p_supernova)
            {
                particle->SN_pending = true;
                particle->type = 2;
            }
            else
            {
                particle->SN_pending = false;
            }
        }
    }
}