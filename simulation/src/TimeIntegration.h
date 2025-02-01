#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>


class TimeIntegration
{
public:
    TimeIntegration(){}
    ~TimeIntegration(){}

    //Semi implicit Euler
    void Euler(Particle* particle, double deltaTime);

    //Leapfrog
    void Kick(Particle* particle, double deltaTime);
    void Drift(Particle* particle, double deltaTime);

    //integrate the internal energy
    void Ueuler(Particle* particle, double deltaTime);
};