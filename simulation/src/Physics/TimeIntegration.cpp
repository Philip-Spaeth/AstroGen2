#include "TimeIntegration.h"

void TimeIntegration::Euler(Particle* particle, double deltaTime)
{
    // Semi implicit Euler
    particle->velocity = particle->velocity + particle->acc * deltaTime;
    particle->position = particle->position + particle->velocity * deltaTime;
}

void TimeIntegration::Kick(Particle* particle, double deltaTime)
{
    if(std::isnan(particle->acc.x) || std::isnan(particle->acc.y) || std::isnan(particle->acc.z))
    {
        std::cout << "NAN gravacc" << std::endl;
        return;
    }
    // Leapfrog Kick
    particle->velocity = particle->velocity + particle->acc * deltaTime / 2;
}

void TimeIntegration::Drift(Particle* particle, double deltaTime)
{
    // Leapfrog Drift
    particle->position = particle->position + particle->velocity * deltaTime;
}

//integrate the Internal Energy
void TimeIntegration::Ueuler(Particle* particle, double deltaTime)
{
    if (!particle)
    {
        std::cerr << "Error: Particle " << " is not initialized." << std::endl;
        return;
    }
    if(std::isnan(particle->dUdt))
    {
       return; 
    }
    if(std::isnan(particle->U))
    {
        return;
    }
    if((particle->U + particle->dUdt * deltaTime) < 0)
    {
        return;
    }
    particle->U += particle->dUdt * deltaTime;
    particle->dUdt = 0;
}