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

    void sfrRoutine(Particle* particle);
};