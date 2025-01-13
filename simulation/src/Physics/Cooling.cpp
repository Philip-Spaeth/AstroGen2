#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

void Cooling::coolingRoutine(Particle* particle)
{
    double density = particle->rho;
    double temperature = particle->T;
    double mu = 0.6;
    double nH = density / (mu * Constants::prtn);
    double n_e = nH;

    double bremsstrahlungCooling = 1.42e-27 * std::sqrt(temperature) * nH * n_e;
    double recombinationCooling = 1.0e-26 * std::pow(temperature, -0.7) * nH * n_e;

    double lineCooling;
    if (temperature < 1e4)
    {
        lineCooling = 7.5e-19 * std::sqrt(temperature / 1e4) * (1.0 + std::sqrt(temperature / 1e4)) / (1.0 + temperature / 1e4);
    }
    else if (temperature < 1e5)
    {
        lineCooling = 1.27e-21 * std::pow(temperature / 1e4, -0.7) * nH * n_e;
    }
    else
    {
        lineCooling = 1.4e-27 * std::pow(temperature / 1e4, -1.5) * nH * n_e;
    }

    double dielectricRecombinationCooling = 1.5e-22 * std::pow(temperature, -0.5) * std::exp(-470000 / temperature) * nH * n_e;
    double netCoolingRate = bremsstrahlungCooling + recombinationCooling + lineCooling + dielectricRecombinationCooling;

    double dudt = -netCoolingRate / density;
    particle->dUdt = dudt;
}