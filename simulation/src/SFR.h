#pragma once
#include <cmath>
#include <vector>
#include <memory>
#include "Particle.h"
#include "Simulation.h"

class Simulation;

class SFR
{
public:
    SFR() {}
    ~SFR() {}

    double totalSFR = 0;

    void sfrRoutine(Particle* particle, Simulation* sim, std::vector<Particle*>* newStars);
};