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
    SFR();
    ~SFR() {}

    double totalSFR = 0;
    std::mt19937 rng;

    void sfrRoutine(Particle* particle, Simulation* sim, double& newStarMass);
};