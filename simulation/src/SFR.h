#pragma once
#include <cmath>
#include <vector>
#include <memory>
#include "Particle.h"

class SFR
{
public:
    SFR() {}
    ~SFR() {}

    double totalSFR = 0;

    void sfrRoutine(Particle* particle);
};